!=====================================================================
! ef_core.F -- standalone numerical core of the VASP-EF optimizer
!
! Pure Fortran 2008, no VASP dependencies.  This file is compiled both
! inside VASP (as part of the patched source tree) and standalone by
! the test driver in tests/, so it must not USE any VASP module and
! must not contain any preprocessor directives.
!
! Contents:
!   ef_sym_eigh            symmetric eigen-decomposition (LAPACK DSYEVD)
!   ef_prfo_step           partitioned-RFO eigenvector-following step
!   ef_rfo_downhill_shift  Newton solve of the RFO secular equation
!   ef_bofill_update       Bofill Hessian update (TS default)
!   ef_powell_update       Powell Hessian update (alternative)
!   ef_model_hessian_diag  element-parametrized diagonal model Hessian
!   ef_mass_weight_freqs   mass-weighted eigenvalues -> wavenumbers
!   ef_incar_*             minimal INCAR tag readers (no VASP reader dep.)
!
! All math is done in flattened Cartesian vectors of length n = 3*NIONS,
! with the Hessian in eV/A^2 and gradients in eV/A.
!
! SPDX-License-Identifier: MIT
!=====================================================================
MODULE ef_core
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: dq = KIND(1.0D0)

  ! sqrt(eV/(A^2*amu)) -> cm^-1 :  sqrt(e/(m_u * A^2))/(2*pi*c)
  REAL(dq), PARAMETER, PUBLIC :: EF_WAVENUMBER_CONST = 521.4700_dq

  ! tolerances
  REAL(dq), PARAMETER :: TINY_LEN = 1.0D-12

  PUBLIC :: ef_sym_eigh
  PUBLIC :: ef_seed_mode
  PUBLIC :: ef_prfo_step
  PUBLIC :: ef_rfo_downhill_shift
  PUBLIC :: ef_bofill_update
  PUBLIC :: ef_powell_update
  PUBLIC :: ef_model_hessian_diag
  PUBLIC :: ef_model_hessian_full
  PUBLIC :: ef_rcov
  PUBLIC :: ef_external_penalty
  PUBLIC :: ef_mass_weight_freqs
  PUBLIC :: ef_incar_logical
  PUBLIC :: ef_incar_int
  PUBLIC :: ef_incar_int_list
  PUBLIC :: ef_incar_real
  PUBLIC :: ef_incar_str

CONTAINS

!----------------------------------------------------------------------
! Lock the curvature along unit vector v to exactly k:
!
!     H' = k v v^T + (I - v v^T) H (I - v v^T)
!
! v becomes an exact eigenvector of H' with eigenvalue k; the
! orthogonal block is the projection of H.  Used to give a user-supplied
! MODECAR direction physical meaning inside the (otherwise atom-local)
! model Hessian.
!----------------------------------------------------------------------
  SUBROUTINE ef_seed_mode(n, H, v, k)
    INTEGER,INTENT(IN)    :: n
    REAL(dq),INTENT(INOUT):: H(n,n)
    REAL(dq),INTENT(IN)   :: v(n), k

    INTEGER  :: i, j
    REAL(dq) :: hv(n), c

    hv = MATMUL(H, v)
    c  = DOT_PRODUCT(v, hv)
    DO j = 1, n
      DO i = 1, n
        H(i,j) = H(i,j) - v(i)*hv(j) - hv(i)*v(j) + (c + k)*v(i)*v(j)
      END DO
    END DO
    ! remove round-off asymmetry
    DO j = 1, n
      DO i = j+1, n
        H(i,j) = 0.5_dq*(H(i,j) + H(j,i))
        H(j,i) = H(i,j)
      END DO
    END DO
  END SUBROUTINE ef_seed_mode

!----------------------------------------------------------------------
! Symmetric eigen-decomposition  H = evec * diag(eval) * evec^T
! eval ascending.  Only the lower triangle of H is referenced.
!----------------------------------------------------------------------
  SUBROUTINE ef_sym_eigh(n, H, eval, evec, info)
    INTEGER,INTENT(IN)    :: n
    REAL(dq),INTENT(IN)   :: H(n,n)
    REAL(dq),INTENT(OUT)  :: eval(n), evec(n,n)
    INTEGER,INTENT(OUT)   :: info        ! LAPACK info; 0 on success

    CHARACTER(1)          :: jobz, uplo
    INTEGER               :: lwork, liwork, liw
    REAL(dq),ALLOCATABLE  :: work(:), a(:,:)
    INTEGER,ALLOCATABLE   :: iwork(:)

    INTERFACE
      SUBROUTINE dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, &
                        liwork, info)
        INTEGER,INTENT(IN) :: n, lda, lwork, liwork
        INTEGER,INTENT(OUT):: info
        CHARACTER(1)       :: jobz, uplo
        REAL(8)            :: a(lda,*), w(*), work(*)
        INTEGER            :: iwork(*)
      END SUBROUTINE dsyevd
    END INTERFACE

    info = 0
    IF (n <= 0) THEN
      info = -1
      RETURN
    END IF

    ALLOCATE(a(n,n))
    a = H                       ! DSYEVD destroys its input
    jobz = 'V'; uplo = 'L'

    ! workspace query
    lwork = -1; liwork = -1
    ALLOCATE(work(1), iwork(1))
    CALL dsyevd(jobz, uplo, n, a, n, eval, work, lwork, iwork, liwork, info)
    lwork  = MAX(INT(work(1)), 1)
    liwork = MAX(iwork(1), 1)
    DEALLOCATE(work, iwork)

    liw = liwork
    ALLOCATE(work(lwork), iwork(liw))
    CALL dsyevd(jobz, uplo, n, a, n, eval, work, lwork, iwork, liw, info)

    evec = a
    DEALLOCATE(a, work, iwork)
  END SUBROUTINE ef_sym_eigh

!----------------------------------------------------------------------
! Newton solve of the downhill RFO secular equation
!
!     sum_{i /= iskip} gp_i^2 / (b_i - lam)^2 = target ,  lam < min(b_i)
!
! gp are gradient projections on the eigenvectors, b the eigenvalues.
! The left-hand side is strictly increasing in lam on (-inf, bmin),
! so the root is unique.  ok=.FALSE. if the downhill block carries no
! gradient (nothing to minimize).
!----------------------------------------------------------------------
  SUBROUTINE ef_rfo_downhill_shift(n, b, gp, iskip, target, lam, ok)
    INTEGER,INTENT(IN)   :: n, iskip
    REAL(dq),INTENT(IN)  :: b(n), gp(n), target
    REAL(dq),INTENT(OUT) :: lam
    LOGICAL,INTENT(OUT)  :: ok

    INTEGER  :: i, it
    REAL(dq) :: ssum, bmin, f, df, dl, eps

    ok = .FALSE.
    ssum = 0.0_dq
    bmin = HUGE(1.0_dq)
    DO i = 1, n
      IF (i == iskip) CYCLE
      ssum = ssum + gp(i)**2
      IF (b(i) < bmin) bmin = b(i)
    END DO
    IF (ssum <= TINY_LEN) RETURN          ! no gradient outside the TS mode

    eps = MAX(1.0E-12_dq, 1.0E-12_dq*ABS(bmin))
    lam = bmin - MAX(1.0_dq, ABS(bmin))   ! safe starting point

    DO it = 1, 200
      f = -target
      df = 0.0_dq
      DO i = 1, n
        IF (i == iskip) CYCLE
        dl = b(i) - lam
        f  = f  + gp(i)**2 / dl**2
        df = df + 2.0_dq * gp(i)**2 / dl**3
      END DO
      IF (ABS(f) <= 1.0E-12_dq*MAX(1.0_dq,target)) EXIT
      IF (df <= 0.0_dq) EXIT              ! cannot happen; safety net
      lam = lam - f/df
      IF (lam >= bmin - eps) lam = 0.5_dq*(lam + (bmin - eps))  ! bisection-ish guard
    END DO
    ok = .TRUE.
  END SUBROUTINE ef_rfo_downhill_shift

!----------------------------------------------------------------------
! Partitioned-RFO eigenvector-following step (Baker 1986 / Cerjan-Miller
! style, with classical length capping):
!
!   climb=.TRUE.  : mode `imode` is followed UPHILL with a unit-length
!        component (the P-RFO uphill shift lam_p = b_m + |g_m|
!        reproduces alpha_m = +sign(g_m), i.e. motion along the gradient
!        on that mode); all other modes are minimized with the downhill
!        shift lam_n solved from the secular equation.
!   climb=.FALSE. : pure RFO-damped minimization -- the secular equation
!        runs over ALL modes including `imode`.
!
!   If the resulting step exceeds `maxstep`, it is rescaled to
!   `maxstep` (direction preserved).
!
! ierr: 0 = ok, 1 = zero gradient.
!----------------------------------------------------------------------
  SUBROUTINE ef_prfo_step(n, eval, evec, grad, imode, climb, maxstep, step, slen, ierr)
    INTEGER,INTENT(IN)   :: n, imode
    LOGICAL,INTENT(IN)   :: climb
    REAL(dq),INTENT(IN)  :: eval(n), evec(n,n), grad(n), maxstep
    REAL(dq),INTENT(OUT) :: step(n), slen
    INTEGER,INTENT(OUT)  :: ierr

    REAL(dq),ALLOCATABLE :: gp(:), alpha(:)
    REAL(dq)             :: lam, r2
    INTEGER              :: i

    ierr = 0
    slen = 0.0_dq
    step = 0.0_dq
    IF (ABS(SQRT(DOT_PRODUCT(grad,grad))) <= TINY_LEN) THEN
      ierr = 1
      RETURN
    END IF

    ALLOCATE(gp(n), alpha(n))
    gp = MATMUL(TRANSPOSE(evec), grad)
    alpha = 0.0_dq
    r2 = MIN(1.0_dq, maxstep)**2

    IF (climb) THEN
      ! single-shift eigenvector following: lam in (b_m, b_next).
      ! The TS-mode component is uphill (b_m - lam < 0), every other
      ! mode downhill, and the secular equation fixes the climb/relax
      ! balance automatically at the requested step length.
      CALL ef_shift_interval(n, eval, gp, imode, r2, lam)
    ELSE
      ! damped quasi-Newton minimization: lam below the lowest
      ! eigenvalue, solving the same equation with span doubling.
      CALL ef_shift_below_min(n, eval, gp, r2, lam)
    END IF
    DO i = 1, n
      alpha(i) = -gp(i) / (eval(i) - lam)
    END DO

    step = MATMUL(evec, alpha)
    slen = SQRT(DOT_PRODUCT(step, step))
    IF (slen > maxstep .AND. slen > 0.0_dq) THEN
      step = step * (maxstep/slen)
      slen = maxstep
    END IF
    DEALLOCATE(gp, alpha)
  END SUBROUTINE ef_prfo_step

!----------------------------------------------------------------------
! CLIMB shift: bracket lam in (b_m, b_next) and bisect the secular
! equation  len(lam)^2 = tgt .  len(lam) is strictly decreasing in lam
! on the interval, so the root is unique.
!----------------------------------------------------------------------
  SUBROUTINE ef_shift_interval(n, b, gp, imode, tgt, lam)
    INTEGER,INTENT(IN)   :: n, imode
    REAL(dq),INTENT(IN)  :: b(n), gp(n), tgt
    REAL(dq),INTENT(OUT) :: lam

    INTEGER  :: it, i
    REAL(dq) :: b_next, eps, lo, hi, hi_orig, m1, m2, l1, l2, lmin, lam_min
    REAL(dq) :: fa, fb

    b_next = HUGE(1.0_dq)
    DO i = 1, n
      IF (i == imode) CYCLE
      IF (b(i) < b_next) b_next = b(i)
    END DO
    eps = MAX(1.0E-8_dq, 1.0E-8_dq*ABS(b_next))
    lo  = b(imode) + eps
    hi  = b_next - eps
    IF (hi <= lo) hi = lo + 1.0E-8_dq
    hi_orig = hi

    ! len^2 is U-shaped on (lo,hi): locate its minimum by ternary search
    DO it = 1, 200
      IF (hi - lo < 1.0E-12_dq*MAX(1.0_dq,ABS(hi))) EXIT
      m1 = lo + (hi - lo)/3.0_dq
      m2 = hi - (hi - lo)/3.0_dq
      CALL ef_shift_len(n, b, gp, m1, l1)
      CALL ef_shift_len(n, b, gp, m2, l2)
      IF (l1 < l2) THEN
        hi = m2
      ELSE
        lo = m1
      END IF
    END DO
    lam_min = 0.5_dq*(lo + hi)
    CALL ef_shift_len(n, b, gp, lam_min, lmin)

    IF (tgt >= lmin) THEN
      ! target reachable: take the UPPER root (most damped step that
      ! still has length sqrt(tgt)); len is increasing on (lam_min, hi_orig)
      lo = lam_min
      hi = hi_orig
      fa = lmin - tgt
      fb = 0.0_dq
      DO it = 1, 200
        lam = 0.5_dq*(lo + hi)
        CALL ef_shift_len(n, b, gp, lam, fb)
        fb = fb - tgt
        IF (ABS(fb) <= 1.0E-12_dq*MAX(1.0_dq,tgt)) EXIT
        IF (fb*fa > 0.0_dq) THEN
          lo = lam
        ELSE
          hi = lam
        END IF
      END DO
    ELSE
      ! target shorter than any achievable step: use the optimal
      ! (minimum-length) shift; the outer cap rescales to maxstep
      lam = lam_min
    END IF
  END SUBROUTINE ef_shift_interval

!----------------------------------------------------------------------
! MIN shift: lam < b_min solving len(lam)^2 = tgt (span doubling plus
! bisection; len decreases as lam increases towards b_min).
!----------------------------------------------------------------------
  SUBROUTINE ef_shift_below_min(n, b, gp, tgt, lam)
    INTEGER,INTENT(IN)   :: n
    REAL(dq),INTENT(IN)  :: b(n), gp(n), tgt
    REAL(dq),INTENT(OUT) :: lam

    INTEGER  :: it
    REAL(dq) :: bmin, span, lo, hi, f

    bmin = MINVAL(b)
    hi   = bmin - MAX(1.0E-8_dq, 1.0E-8_dq*ABS(bmin))
    CALL ef_shift_len(n, b, gp, hi, f)
    IF (f <= tgt) THEN
      lam = hi
      RETURN
    END IF
    span = MAX(1.0_dq, ABS(bmin))
    lo = hi - span
    DO it = 1, 60
      CALL ef_shift_len(n, b, gp, lo, f)
      IF (f <= tgt) EXIT
      lo = hi - span
      span = span*4.0_dq
      lo = hi - span
    END DO
    DO it = 1, 200
      lam = 0.5_dq*(lo + hi)
      CALL ef_shift_len(n, b, gp, lam, f)
      IF (ABS(f - tgt) <= 1.0E-12_dq*MAX(1.0_dq,tgt)) EXIT
      IF (f > tgt) THEN
        hi = lam
      ELSE
        lo = lam
      END IF
    END DO
  END SUBROUTINE ef_shift_below_min

  ! squared length of the shifted step for a given lam
  SUBROUTINE ef_shift_len(n, b, gp, lam, len2)
    INTEGER,INTENT(IN)   :: n
    REAL(dq),INTENT(IN)  :: b(n), gp(n), lam
    REAL(dq),INTENT(OUT) :: len2
    INTEGER :: i
    REAL(dq) :: dl
    len2 = 0.0_dq
    DO i = 1, n
      dl = b(i) - lam
      len2 = len2 + gp(i)**2/dl**2
    END DO
  END SUBROUTINE ef_shift_len

!----------------------------------------------------------------------
! Bofill update (interpolation between SR1 and PSB; the standard TS
! Hessian update, cf. Bofill, J. Comput. Chem. 15, 1 (1994)):
!
!   eta = y - H s
!   H  += (eta s^T + s eta^T)/(s.s) - (s.eta) s s^T/(s.s)^2
!   H  += phi * [ eta eta^T/(s.y) - (s.eta) s s^T/(s.s)^2 ]
!   phi = (s.y)^2 / ((s.s)(y.y))
!
! applied=.FALSE. if the update is undefined (degenerate s or y).
!----------------------------------------------------------------------
  SUBROUTINE ef_bofill_update(n, H, s, y, applied)
    INTEGER,INTENT(IN)    :: n
    REAL(dq),INTENT(INOUT):: H(n,n)
    REAL(dq),INTENT(IN)   :: s(n), y(n)
    LOGICAL,INTENT(OUT)   :: applied

    REAL(dq),ALLOCATABLE  :: eta(:), ss(:,:), etas(:,:), yy(:,:)
    REAL(dq)              :: d1, d2, seta, yyv, phi
    INTEGER               :: i, j

    applied = .FALSE.
    d1 = DOT_PRODUCT(s, s)
    d2 = DOT_PRODUCT(s, y)
    IF (d1 <= TINY_LEN .OR. ABS(d2) <= TINY_LEN) RETURN
    ALLOCATE(eta(n), ss(n,n), etas(n,n), yy(n,n))
    eta = y - MATMUL(H, s)
    seta = DOT_PRODUCT(s, eta)
    yyv = DOT_PRODUCT(y, y)
    IF (yyv <= TINY_LEN) THEN
      DEALLOCATE(eta, ss, etas, yy)
      RETURN
    END IF
    phi = d2**2 / (d1*yyv)

    DO j = 1, n
      DO i = 1, n
        ss(i,j)  = s(i)*s(j)
        etas(i,j)= eta(i)*s(j) + s(i)*eta(j)
        yy(i,j)  = eta(i)*eta(j)
      END DO
    END DO
    H = H + etas/d1 - seta*ss/d1**2 &
          + phi*(yy/d2 - seta*ss/d1**2)

    ! remove round-off asymmetry
    DO j = 1, n
      DO i = j+1, n
        H(i,j) = 0.5_dq*(H(i,j) + H(j,i))
        H(j,i) = H(i,j)
      END DO
    END DO
    DEALLOCATE(eta, ss, etas, yy)
    applied = .TRUE.
  END SUBROUTINE ef_bofill_update

!----------------------------------------------------------------------
! Powell update: like PSB, but the difference vector is modified to
! keep the curvature along the step negative (required near a saddle):
!
!   if (s.y > 0)  y' = y - 2 (s.y)/(s.s) s
!   eta = y' - H s ;  H += (eta s^T + s eta^T)/(s.s) - (s.eta) s s^T/(s.s)^2
!----------------------------------------------------------------------
  SUBROUTINE ef_powell_update(n, H, s, y, applied)
    INTEGER,INTENT(IN)    :: n
    REAL(dq),INTENT(INOUT):: H(n,n)
    REAL(dq),INTENT(IN)   :: s(n)
    REAL(dq),INTENT(INOUT):: y(n)
    LOGICAL,INTENT(OUT)   :: applied

    REAL(dq),ALLOCATABLE  :: eta(:), ss(:,:), etas(:,:)
    REAL(dq)              :: d1, d2, seta
    INTEGER               :: i, j

    applied = .FALSE.
    d1 = DOT_PRODUCT(s, s)
    d2 = DOT_PRODUCT(s, y)
    IF (d1 <= TINY_LEN) RETURN
    IF (d2 > 0.0_dq) y = y - (2.0_dq*d2/d1)*s
    d2 = DOT_PRODUCT(s, y)
    IF (ABS(d2) <= TINY_LEN) RETURN

    ALLOCATE(eta(n), ss(n,n), etas(n,n))
    eta = y - MATMUL(H, s)
    seta = DOT_PRODUCT(s, eta)
    DO j = 1, n
      DO i = 1, n
        ss(i,j)  = s(i)*s(j)
        etas(i,j)= eta(i)*s(j) + s(i)*eta(j)
      END DO
    END DO
    H = H + etas/d1 - seta*ss/d1**2

    DO j = 1, n
      DO i = j+1, n
        H(i,j) = 0.5_dq*(H(i,j) + H(j,i))
        H(j,i) = H(i,j)
      END DO
    END DO
    DEALLOCATE(eta, ss, etas)
    applied = .TRUE.
  END SUBROUTINE ef_powell_update

!----------------------------------------------------------------------
! Element-parametrized diagonal model Hessian (Fischer & Almlof 1992
! spirit): one characteristic force constant per atom, in eV/A^2,
! replicated on all three Cartesian coordinates and scaled by `scale`
! (INCAR: EF_HSCALE).  Only the order of magnitude matters -- the
! Bofill/Powell updates rebuild the rest of the Hessian from forces.
!----------------------------------------------------------------------
  SUBROUTINE ef_model_hessian_diag(nat, zat, scale, kdiag)
    INTEGER,INTENT(IN)   :: nat
    INTEGER,INTENT(IN)   :: zat(nat)     ! atomic numbers
    REAL(dq),INTENT(IN)  :: scale
    REAL(dq),INTENT(OUT) :: kdiag(3,nat)

    INTEGER :: ia
    REAL(dq) :: k

    DO ia = 1, nat
      k = element_k(zat(ia)) * scale
      kdiag(1,ia) = k
      kdiag(2,ia) = k
      kdiag(3,ia) = k
    END DO
  END SUBROUTINE ef_model_hessian_diag

  ! characteristic diagonal force constant, eV/A^2, by atomic number.
  ! Scale chosen to match typical bond-stretch force constants, so the
  ! RFO-damped steps do not over-amplify stiff directions.
  FUNCTION element_k(z) RESULT(k)
    INTEGER,INTENT(IN) :: z
    REAL(dq)           :: k

    SELECT CASE (z)
      CASE (1)              ; k = 40.0_dq    ! H
      CASE (2)              ; k = 20.0_dq    ! He
      CASE (3:10)           ; k = 16.0_dq + 4.0_dq*(z-3)   ! Li .. Ne
      CASE (11:18)          ; k = 12.0_dq + 3.0_dq*(z-11)  ! Na .. Ar
      CASE (19:36)          ; k =  8.0_dq + 1.0_dq*(z-19)  ! K  .. Kr
      CASE (37:54)          ; k =  6.0_dq + 0.5_dq*(z-37)  ! Rb .. Xe
      CASE DEFAULT          ; k = 16.0_dq
    END SELECT
  END FUNCTION element_k

!----------------------------------------------------------------------
! Bond-aware model Hessian for periodic systems:
!   * element-parametrized diagonal terms (as above), plus
!   * one spring (rank-1 stiffness along the bond direction) for every
!     covalently bonded pair, detected with the covalent-radius sum
!     criterion under the minimum-image convention.
! This makes the model modes physical on surfaces: substrate atoms are
! stiffened by their bonds, adsorbate modes stay soft.
!----------------------------------------------------------------------
  SUBROUTINE ef_model_hessian_full(nat, zat, cart, a, b, bondk, btol, scale, H)
    INTEGER,INTENT(IN)   :: nat
    INTEGER,INTENT(IN)   :: zat(nat)
    REAL(dq),INTENT(IN)  :: cart(3,nat), a(3,3), b(3,3)
    REAL(dq),INTENT(IN)  :: bondk, btol, scale
    REAL(dq),INTENT(OUT) :: H(3*nat,3*nat)

    INTEGER  :: n, ia, ib, i, j, ka, kb
    REAL(dq) :: kd(3,nat), df(3), dc(3), d, rc, u(3), k

    n = 3*nat
    CALL ef_model_hessian_diag(nat, zat, scale, kd)
    H = 0.0_dq
    DO ia = 1, nat
      DO i = 1, 3
        H(3*(ia-1)+i, 3*(ia-1)+i) = kd(i,ia)
      END DO
    END DO

    DO ia = 1, nat-1
      DO ib = ia+1, nat
        dc = cart(:,ib) - cart(:,ia)
        df = MATMUL(dc, b)
        df = df - ANINT(df)                    ! minimum image (fractional)
        dc = MATMUL(df, a)
        d  = SQRT(DOT_PRODUCT(dc, dc))
        rc = ef_rcov(zat(ia)) + ef_rcov(zat(ib)) + btol
        IF (d > 1.0E-6_dq .AND. d < rc) THEN
          u = dc/d
          k  = bondk*scale
          ! spring in block form: [[+k nnT, -k nnT], [-k nnT, +k nnT]]
          DO i = 1, 3
            ka = 3*(ia-1)+i
            kb = 3*(ib-1)+i
            DO j = 1, 3
              H(ka, 3*(ia-1)+j) = H(ka, 3*(ia-1)+j) + k*u(i)*u(j)
              H(kb, 3*(ib-1)+j) = H(kb, 3*(ib-1)+j) + k*u(i)*u(j)
              H(ka, 3*(ib-1)+j) = H(ka, 3*(ib-1)+j) - k*u(i)*u(j)
              H(kb, 3*(ia-1)+j) = H(kb, 3*(ia-1)+j) - k*u(i)*u(j)
            END DO
          END DO
        END IF
      END DO
    END DO
  END SUBROUTINE ef_model_hessian_full

  ! covalent radii, A, by atomic number (standard values)
  FUNCTION ef_rcov(z) RESULT(r)
    INTEGER,INTENT(IN) :: z
    REAL(dq)           :: r
    SELECT CASE (z)
      CASE (1);  r = 0.31_dq
      CASE (2);  r = 0.28_dq
      CASE (3:10)
        r = 0.28_dq + 0.04_dq*(z-2)          ! Li .. Ne (rough)
        IF (z==3) r = 1.28_dq
        IF (z==4) r = 0.96_dq
        IF (z==5) r = 0.84_dq
        IF (z==6) r = 0.76_dq
        IF (z==7) r = 0.71_dq
        IF (z==8) r = 0.66_dq
        IF (z==9) r = 0.57_dq
      CASE (11:18)
        r = 1.66_dq - 0.06_dq*(z-11)         ! Na .. Ar
      CASE (19:36)
        r = 2.03_dq - 0.023_dq*(z-19)        ! K .. Kr
      CASE (37:54)
        r = 2.20_dq - 0.021_dq*(z-37)        ! Rb .. Xe
      CASE DEFAULT
        r = 1.20_dq
    END SELECT
  END FUNCTION ef_rcov

!----------------------------------------------------------------------
! Penalize the external (rigid) degrees of freedom of a Cartesian
! Hessian: global translations (always exact zero modes of the energy,
! in PBC too) and, for finite molecules in a box, global rotations
! about the mass centre (lrot=.TRUE.).  The rigid directions receive a
! large curvature K so they are never amplified by the RFO shift and
! never selected as the TS mode.  Applied to a WORKING copy before
! diagonalization; the physical Hessian is left untouched.
!----------------------------------------------------------------------
  SUBROUTINE ef_external_penalty(nat, cart, masses, lrot, K, H)
    INTEGER,INTENT(IN)    :: nat
    REAL(dq),INTENT(IN)   :: cart(3,nat), masses(nat)
    LOGICAL,INTENT(IN)    :: lrot
    REAL(dq),INTENT(IN)   :: K
    REAL(dq),INTENT(INOUT):: H(3*nat,3*nat)

    INTEGER  :: n, ia, i, next, kp, p, q
    REAL(dq) :: com(3), E(3*nat,6), r(3), v(3*nat), nrm

    n = 3*nat
    com = 0.0_dq
    DO ia = 1, nat
      com = com + masses(ia)*cart(:,ia)
    END DO
    com = com/MAX(SUM(masses), 1.0E-10_dq)

    E = 0.0_dq
    next = 0
    DO i = 1, 3                              ! translations
      DO ia = 1, nat
        E(3*(ia-1)+i, i) = 1.0_dq/SQRT(REAL(nat,dq))
      END DO
      next = next + 1
    END DO
    IF (lrot) THEN                           ! rotations about x, y, z
      DO i = 1, 3
        DO ia = 1, nat
          r = cart(:,ia) - com
          v(3*(ia-1)+1) = r(2)*AXK(i,3) - r(3)*AXK(i,2)
          v(3*(ia-1)+2) = r(3)*AXK(i,1) - r(1)*AXK(i,3)
          v(3*(ia-1)+3) = r(1)*AXK(i,2) - r(2)*AXK(i,1)
          E(3*(ia-1)+1:3*(ia-1)+3, 3+i) = v(3*(ia-1)+1:3*(ia-1)+3)
        END DO
      END DO
      next = 6
    END IF

    ! Gram-Schmidt orthonormalize the external directions
    DO kp = 1, next
      DO p = 1, kp-1
        nrm = DOT_PRODUCT(E(1:n,p), E(1:n,kp))
        E(1:n,kp) = E(1:n,kp) - nrm*E(1:n,p)
      END DO
      nrm = SQRT(DOT_PRODUCT(E(1:n,kp), E(1:n,kp)))
      IF (nrm > 1.0E-10_dq) E(1:n,kp) = E(1:n,kp)/nrm
    END DO

    ! H += K * E E^T
    DO kp = 1, next
      DO q = 1, n
        DO p = 1, n
          H(p,q) = H(p,q) + K*E(p,kp)*E(q,kp)
        END DO
      END DO
    END DO
  END SUBROUTINE ef_external_penalty

  ! unit axis vector: AXK(axis, comp)
  FUNCTION AXK(axis, comp) RESULT(v)
    INTEGER,INTENT(IN) :: axis, comp
    REAL(dq)           :: v
    v = 0.0_dq
    IF (axis == comp) v = 1.0_dq
  END FUNCTION AXK

!----------------------------------------------------------------------
! Mass-weighted normal-mode analysis of the Cartesian Hessian:
! builds H_mw = M^-1/2 H M^-1/2, diagonalizes it, and converts the
! eigenvalues to wavenumbers (cm^-1); negative eigenvalues are reported
! as negative (imaginary) frequencies.
!
!   masses  in amu, one per atom (VASP POMASS)
!   freqs   (3*nat) wavenumbers, ascending by eigenvalue
!   nimag   number of freqs < -zero_thresh
!----------------------------------------------------------------------
  SUBROUTINE ef_mass_weight_freqs(nat, masses, H, freqs, nimag, zero_thresh, info)
    INTEGER,INTENT(IN)   :: nat
    REAL(dq),INTENT(IN)  :: masses(nat), H(3*nat,3*nat), zero_thresh
    REAL(dq),INTENT(OUT) :: freqs(3*nat)
    INTEGER,INTENT(OUT)  :: nimag, info

    INTEGER              :: n, ia, i, j
    REAL(dq),ALLOCATABLE :: hmw(:,:), evec(:,:), dminv(:)

    n = 3*nat
    ALLOCATE(hmw(n,n), evec(n,n), dminv(n))
    DO ia = 1, nat
      DO i = 0, 2
        dminv(3*(ia-1)+i+1) = 1.0_dq/SQRT(MAX(masses(ia), 1.0E-8_dq))
      END DO
    END DO
    DO j = 1, n
      DO i = 1, n
        hmw(i,j) = dminv(i)*H(i,j)*dminv(j)
      END DO
    END DO
    CALL ef_sym_eigh(n, hmw, freqs, evec, info)
    IF (info /= 0) RETURN

    nimag = 0
    DO i = 1, n
      IF (freqs(i) >= 0.0_dq) THEN
        freqs(i) = EF_WAVENUMBER_CONST*SQRT(freqs(i))
      ELSE
        freqs(i) = -EF_WAVENUMBER_CONST*SQRT(-freqs(i))
        IF (freqs(i) < -zero_thresh) nimag = nimag + 1
      END IF
    END DO
    DEALLOCATE(hmw, evec, dminv)
  END SUBROUTINE ef_mass_weight_freqs

!----------------------------------------------------------------------
! Minimal INCAR readers.  INCAR syntax understood: one KEY = VALUE per
! line, '#' or ';' starts a trailing comment, case-insensitive keys.
! These bypass VASP's RDATAB so that the optimizer is self-contained.
!----------------------------------------------------------------------

  ! locate the right-hand side of `key` in file `fname`; internal use
  SUBROUTINE incar_lookup(fname, key, rhs, found)
    CHARACTER(*),INTENT(IN)  :: fname, key
    CHARACTER(512),INTENT(OUT):: rhs
    LOGICAL,INTENT(OUT)      :: found

    CHARACTER(512) :: line, k
    INTEGER        :: unit, ios, eq, cmt
    LOGICAL        :: exists

    found = .FALSE.; rhs = ''
    INQUIRE(FILE=fname, EXIST=exists)
    IF (.NOT. exists) RETURN
    OPEN(NEWUNIT=unit, FILE=fname, STATUS='OLD', ACTION='READ', &
         IOSTAT=ios)
    IF (ios /= 0) RETURN
    DO
      READ(unit, '(A)', IOSTAT=ios) line
      IF (ios /= 0) EXIT
      eq = SCAN(line, '=')
      IF (eq <= 1) CYCLE                     ! no '=' or empty key
      k = UPCASE(ADJUSTL(line(1:eq-1)))
      IF (TRIM(k) /= UPCASE(TRIM(key))) CYCLE
      rhs = line(eq+1:)
      cmt = SCAN(rhs, '#;!')
      IF (cmt > 0) rhs = rhs(1:cmt-1)
      found = .TRUE.
      EXIT
    END DO
    CLOSE(unit)
  END SUBROUTINE incar_lookup

  SUBROUTINE ef_incar_logical(fname, key, val, found)
    CHARACTER(*),INTENT(IN) :: fname, key
    LOGICAL,INTENT(OUT)     :: val
    LOGICAL,INTENT(OUT)     :: found
    CHARACTER(512)          :: rhs
    CALL incar_lookup(fname, key, rhs, found)
    val = .FALSE.
    IF (.NOT. found) RETURN
    SELECT CASE (UPCASE(ADJUSTL(rhs)))
      CASE ('.TRUE.','.T.','TRUE','T')
        val = .TRUE.
      CASE ('.FALSE.','.F.','FALSE','F')
        val = .FALSE.
      CASE DEFAULT
        found = .FALSE.        ! unparsable -> treat as absent
    END SELECT
  END SUBROUTINE ef_incar_logical

  SUBROUTINE ef_incar_int(fname, key, val, found)
    CHARACTER(*),INTENT(IN) :: fname, key
    INTEGER,INTENT(OUT)     :: val
    LOGICAL,INTENT(OUT)     :: found
    CHARACTER(512)          :: rhs
    INTEGER                 :: ios
    CALL incar_lookup(fname, key, rhs, found)
    val = 0
    IF (.NOT. found) RETURN
    READ(rhs, *, IOSTAT=ios) val
    IF (ios /= 0) found = .FALSE.
  END SUBROUTINE ef_incar_int

  SUBROUTINE ef_incar_real(fname, key, val, found)
    CHARACTER(*),INTENT(IN) :: fname, key
    REAL(dq),INTENT(OUT)    :: val
    LOGICAL,INTENT(OUT)     :: found
    CHARACTER(512)          :: rhs
    INTEGER                 :: ios
    CALL incar_lookup(fname, key, rhs, found)
    val = 0.0_dq
    IF (.NOT. found) RETURN
    READ(rhs, *, IOSTAT=ios) val
    IF (ios /= 0) found = .FALSE.
  END SUBROUTINE ef_incar_real

  ! first whitespace-delimited token of the value, uppercased
  SUBROUTINE ef_incar_str(fname, key, val, found)
    CHARACTER(*),INTENT(IN) :: fname, key
    CHARACTER(64),INTENT(OUT):: val
    LOGICAL,INTENT(OUT)     :: found
    CHARACTER(512)          :: rhs
    INTEGER                 :: sp
    CALL incar_lookup(fname, key, rhs, found)
    val = ''
    IF (.NOT. found) RETURN
    rhs = ADJUSTL(rhs)
    sp = SCAN(rhs, ' '//CHAR(9))
    IF (sp > 1) rhs = rhs(1:sp-1)
    val = UPCASE(TRIM(rhs))
  END SUBROUTINE ef_incar_str

  ! whitespace-separated list of integers (e.g. EF_ACTIVE = 1 2 5)
  SUBROUTINE ef_incar_int_list(fname, key, maxn, list, n, found)
    CHARACTER(*),INTENT(IN) :: fname, key
    INTEGER,INTENT(IN)      :: maxn
    INTEGER,INTENT(OUT)     :: list(maxn)
    INTEGER,INTENT(OUT)     :: n
    LOGICAL,INTENT(OUT)     :: found

    CHARACTER(512)          :: rhs, tok
    INTEGER                 :: pos, v, ios, p0, ln

    CALL incar_lookup(fname, key, rhs, found)
    list = 0
    n = 0
    IF (.NOT. found) RETURN
    found = .FALSE.
    pos = 1
    ln = LEN(rhs)
    DO
      DO WHILE (pos <= ln .AND. rhs(pos:pos) == ' ')
        pos = pos + 1
      END DO
      p0 = pos
      DO WHILE (pos <= ln .AND. rhs(pos:pos) /= ' ')
        pos = pos + 1
      END DO
      IF (pos-1 < p0) EXIT
      tok = rhs(p0:pos-1)
      READ(tok, *, IOSTAT=ios) v
      IF (ios /= 0) CYCLE
      n = n + 1
      IF (n <= maxn) list(n) = v
      found = .TRUE.
    END DO
  END SUBROUTINE ef_incar_int_list

  FUNCTION upcase(s) RESULT(u)
    CHARACTER(*),INTENT(IN) :: s
    CHARACTER(LEN(s))       :: u
    INTEGER :: i, ic
    u = s
    DO i = 1, LEN(s)
      ic = IACHAR(s(i:i))
      IF (ic >= IACHAR('a') .AND. ic <= IACHAR('z')) u(i:i) = ACHAR(ic-32)
    END DO
  END FUNCTION upcase

END MODULE ef_core
