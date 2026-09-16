!=====================================================================
! test_ef_core.F90 -- standalone test driver for ef_core.F
!
! Build & run (Linux / MSYS2 UCRT64):
!   gfortran -O2 -o test_ef_core ../src/ef_core.F test_ef_core.F90 -llapack -lblas
!   ./test_ef_core
!
! All checks are property-based (no golden files required); a nonzero
! exit status means at least one test failed.
!=====================================================================
PROGRAM test_ef_core
  USE ef_core
  IMPLICIT NONE

  INTEGER :: nfail
  nfail = 0

  CALL test_eigh(nfail)
  CALL test_seed(nfail)
  CALL test_rfo_shift(nfail)
  CALL test_prfo(nfail)
  CALL test_updates(nfail)
  CALL test_model_hessian(nfail)
  CALL test_model_bonds(nfail)
  CALL test_extpen(nfail)
  CALL test_freqs(nfail)
  CALL test_incar(nfail)

  IF (nfail == 0) THEN
    WRITE(*,'(A)') 'ALL TESTS PASSED'
  ELSE
    WRITE(*,'(A,I0,A)') 'FAILED: ', nfail, ' check(s)'
    STOP 1
  END IF

CONTAINS

  SUBROUTINE check(cond, name)
    LOGICAL,INTENT(IN) :: cond
    CHARACTER(*)       :: name
    IF (cond) THEN
      WRITE(*,'(A,A)') '  ok   ', name
    ELSE
      WRITE(*,'(A,A)') '  FAIL ', name
      nfail = nfail + 1
    END IF
  END SUBROUTINE check

  SUBROUTINE test_eigh(nfail)
    INTEGER,INTENT(INOUT) :: nfail
    INTEGER,PARAMETER :: n = 6
    REAL(dq) :: H(n,n), eval(n), evec(n,n), R(n,n)
    INTEGER  :: i, j, info
    CALL RANDOM_SEED()
    DO j = 1, n
      DO i = 1, j
        CALL RANDOM_NUMBER(H(i,j))
        H(i,j) = 2*H(i,j) - 1
        H(j,i) = H(i,j)
      END DO
    END DO
    CALL ef_sym_eigh(n, H, eval, evec, info)
    CALL check(info == 0, 'eigh: info==0')
    CALL check(ALL(eval(2:n) >= eval(1:n-1) - 1.0D-12), 'eigh: ascending')
    R = MATMUL(H, evec) - evec*SPREAD(eval, 1, n)
    CALL check(MAXVAL(ABS(R)) < 1.0D-12, 'eigh: H*v = lam*v')
    R = MATMUL(TRANSPOSE(evec), evec) - eye(n)
    CALL check(MAXVAL(ABS(R)) < 1.0D-12, 'eigh: orthonormal')
  END SUBROUTINE test_eigh

  FUNCTION eye(n) RESULT(res)
    INTEGER,INTENT(IN) :: n
    REAL(dq)           :: res(n,n)
    INTEGER :: k
    res = 0.0_dq
    DO k = 1, n
      res(k,k) = 1.0_dq
    END DO
  END FUNCTION eye

  SUBROUTINE test_seed(nfail)
    INTEGER,INTENT(INOUT) :: nfail
    INTEGER,PARAMETER :: n = 6
    REAL(dq) :: H(n,n), H0(n,n), v(n), Hv(n), w(n)
    INTEGER  :: i, j
    REAL(dq),PARAMETER :: K = -0.5_dq

    DO j = 1, n
      DO i = 1, j
        CALL RANDOM_NUMBER(H0(i,j)); H0(i,j) = 2*H0(i,j)-1
        H0(j,i) = H0(i,j)
      END DO
    END DO
    CALL RANDOM_NUMBER(v)
    v = v/SQRT(DOT_PRODUCT(v,v))
    H = H0
    CALL ef_seed_mode(n, H, v, K)

    Hv = MATMUL(H, v)                              ! exact eigenpair
    CALL check(MAXVAL(ABS(Hv - K*v)) < 1.0D-12, 'seed: H v = k v')
    CALL check(ABS(DOT_PRODUCT(v, MATMUL(H, v)) - K) < 1.0D-12, &
               'seed: v^T H v = k')

    CALL RANDOM_NUMBER(w)                          ! orthogonal test vector
    w = w - DOT_PRODUCT(v, w)*v
    w = w/SQRT(DOT_PRODUCT(w,w))
    CALL check(ABS(DOT_PRODUCT(w, MATMUL(H, w)) - &
                   DOT_PRODUCT(w, MATMUL(H0, w))) < 1.0D-12, &
               'seed: orthogonal block preserved')
    CALL check(MAXVAL(ABS(H - TRANSPOSE(H))) < 1.0D-14, 'seed: symmetric')
  END SUBROUTINE test_seed

  SUBROUTINE test_rfo_shift(nfail)
    INTEGER,INTENT(INOUT) :: nfail
    INTEGER,PARAMETER :: n = 4
    REAL(dq) :: b(n), gp(n), lam, f
    INTEGER  :: i, iskip
    LOGICAL  :: ok
    ! known eigen-spectrum; follow mode 1 uphill, minimize modes 2..4
    b  = (/ -2.0_dq, 1.0_dq, 3.0_dq, 5.0_dq /)
    gp = (/  0.7_dq, 1.0_dq, 0.5_dq, 0.2_dq /)
    iskip = 1
    CALL ef_rfo_downhill_shift(n, b, gp, iskip, 1.0_dq, lam, ok)
    CALL check(ok, 'rfo: solved')
    CALL check(lam < MINVAL(b(2:n)), 'rfo: lam below bmin')
    f = -1.0_dq
    DO i = 2, n
      f = f + gp(i)**2/(b(i)-lam)**2
    END DO
    CALL check(ABS(f) < 1.0D-10, 'rfo: secular residual')
    ! degenerate case: no gradient outside the TS mode
    gp(2:n) = 0.0_dq
    CALL ef_rfo_downhill_shift(n, b, gp, iskip, 1.0_dq, lam, ok)
    CALL check(.NOT. ok, 'rfo: rejects zero downhill gradient')
  END SUBROUTINE test_rfo_shift

  SUBROUTINE test_prfo(nfail)
    INTEGER,INTENT(INOUT) :: nfail
    INTEGER,PARAMETER :: n = 5
    REAL(dq) :: H(n,n), eval(n), evec(n,n), grad(n), step(n), slen
    REAL(dq) :: proj, gp(n)
    INTEGER  :: info, ierr, imode, i
    LOGICAL  :: ok
    REAL(dq) :: lam

    ! build H from a known spectrum
    eval = (/ -1.5_dq, 0.5_dq, 1.2_dq, 2.0_dq, 4.0_dq /)
    CALL RANDOM_NUMBER(evec)
    CALL qr_orthonormalize(n, evec)
    H = MATMUL(evec*SPREAD(eval,1,n), TRANSPOSE(evec))
    H = 0.5_dq*(H + TRANSPOSE(H))
    imode = 1
    CALL RANDOM_NUMBER(grad)
    grad = 2*grad - 1

    ! single-shift EF climb: mode m uphill, others downhill
    CALL ef_prfo_step(n, eval, evec, grad, imode, .TRUE., 1.0E6_dq, step, slen, ierr)
    CALL check(ierr == 0, 'prfo: ierr==0')
    gp = MATMUL(TRANSPOSE(evec), grad)
    proj = DOT_PRODUCT(evec(:,imode), step)
    CALL check(proj*gp(imode) > 0.0_dq, 'prfo: climb uphill along mode')
    CALL check(SQRT(DOT_PRODUCT(step,step)) > 0.0_dq, &
               'prfo: climb step has nonzero length')
    proj = DOT_PRODUCT(evec(:,3), step)
    CALL check(proj*gp(3) < 0.0_dq, 'prfo: other modes damped downhill')
    ! tiny gradient: the secular root exists, length hits R exactly
    CALL ef_prfo_step(n, eval, evec, 1.0E-2_dq*grad, imode, .TRUE., 1.0E6_dq, &
                      step, slen, ierr)
    CALL check(ABS(slen - 1.0_dq) < 1.0D-6, 'prfo: small-g secular length = R')

    ! step capping
    CALL ef_prfo_step(n, eval, evec, grad, imode, .TRUE., 0.3_dq, step, slen, ierr)
    CALL check(slen <= 0.3_dq + 1.0D-12, 'prfo: capped to maxstep')
    ! zero gradient -> error flag
    CALL ef_prfo_step(n, eval, evec, 0.0_dq*grad, imode, .TRUE., 0.3_dq, step, slen, ierr)
    CALL check(ierr == 1, 'prfo: zero-gradient flagged')

    ! minimization mode: every eigencomponent descends
    CALL ef_prfo_step(n, eval, evec, grad, imode, .FALSE., 1.0E6_dq, step, slen, ierr)
    CALL check(DOT_PRODUCT(step, grad) < 0.0_dq, 'prfo: min step descends')
  END SUBROUTINE test_prfo

  SUBROUTINE qr_orthonormalize(n, a)
    ! crude Gram-Schmidt: enough for a test fixture
    INTEGER,INTENT(IN) :: n
    REAL(dq),INTENT(INOUT) :: a(n,n)
    INTEGER :: i, j
    REAL(dq) :: v(n), p
    DO j = 1, n
      DO i = 1, j-1
        p = DOT_PRODUCT(a(:,i), a(:,j))
        a(:,j) = a(:,j) - p*a(:,i)
      END DO
      v = a(:,j)
      a(:,j) = v/SQRT(DOT_PRODUCT(v,v))
    END DO
  END SUBROUTINE qr_orthonormalize

  SUBROUTINE test_updates(nfail)
    INTEGER,INTENT(INOUT) :: nfail
    INTEGER,PARAMETER :: n = 3
    REAL(dq) :: H(n,n), H0(n,n), s(n), y(n), grad(n), x(n), x2(n)
    LOGICAL  :: applied
    INTEGER  :: i, j

    ! Bofill on a quadratic: g(x) = H0 (x - x0); pick x0=0, random H0
    DO j = 1, n
      DO i = 1, j
        CALL RANDOM_NUMBER(H0(i,j)); H0(i,j) = 2*H0(i,j)-1
        H0(j,i) = H0(i,j)
      END DO
    END DO
    CALL RANDOM_NUMBER(x);  grad  = MATMUL(H0, x)
    CALL RANDOM_NUMBER(x2); y     = MATMUL(H0, x2) - grad
    s = x2 - x
    H = eye(n)
    CALL ef_bofill_update(n, H, s, y, applied)
    CALL check(applied, 'bofill: applied')
    CALL check(MAXVAL(ABS(H - TRANSPOSE(H))) < 1.0D-14, 'bofill: symmetric')
    CALL check(ALL(ABS(H) < HUGE(1.0_dq)/2), 'bofill: finite')

    ! n=1 closed form: SR1 weight phi=1 -> H' = h + (y-h)^2/y
    BLOCK
      REAL(dq) :: H1(1,1), s1(1), y1(1), expect
      H1 = RESHAPE((/ 1.0_dq /), (/1,1/))
      s1 = 1.0_dq
      y1 = -2.0_dq
      expect = 1.0_dq + (-3.0_dq)**2/(-2.0_dq)      ! = -3.5
      CALL ef_bofill_update(1, H1, s1, y1, applied)
      CALL check(applied .AND. ABS(H1(1,1) - expect) < 1.0D-12, &
                 'bofill: n=1 closed form')
      ! Powell, negative s.y kept, positive flipped to same result
      H1 = 1.0_dq
      CALL ef_powell_update(1, H1, s1, y1, applied)
      CALL check(applied .AND. ABS(H1(1,1) - (-2.0_dq)) < 1.0D-12, &
                 'powell: n=1 closed form (s.y<0)')
      y1 = 3.0_dq
      H1 = 1.0_dq
      CALL ef_powell_update(1, H1, s1, y1, applied)
      ! y flips to -3, then PSB with eta=-4 gives H' = 1+2*(-4)-(-4) = -3
      CALL check(applied .AND. ABS(H1(1,1) - (-3.0_dq)) < 1.0D-12, &
                 'powell: flips positive curvature along s')
    END BLOCK

    ! degenerate inputs are rejected, not fatal
    CALL ef_bofill_update(n, H, 0.0_dq*s, y, applied)
    CALL check(.NOT. applied, 'bofill: rejects zero step')
  END SUBROUTINE test_updates

  SUBROUTINE test_model_hessian(nfail)
    INTEGER,INTENT(INOUT) :: nfail
    REAL(dq) :: k(3,4)
    CALL ef_model_hessian_diag(4, (/ 1, 6, 7, 8 /), 2.0_dq, k)
    CALL check(ALL(ABS(k(2,:) - k(1,:)) < 1.0D-14 .AND. &
                   ABS(k(3,:) - k(1,:)) < 1.0D-14), 'model: isotropic per atom')
    CALL check(ABS(k(1,1) - 80.0_dq) < 1.0D-14, 'model: H scaled by 2')
    CALL check(k(1,3) > k(1,2), 'model: N stiffer than C')
  END SUBROUTINE test_model_hessian

  SUBROUTINE test_model_bonds(nfail)
    INTEGER,INTENT(INOUT) :: nfail
    INTEGER,PARAMETER :: nat = 2
    INTEGER  :: z(nat)
    REAL(dq) :: cart(3,nat), a(3,3), b(3,3), H(3*nat,3*nat)
    REAL(dq) :: k, dz

    ! cubic 6 A cell, b = a^-1 = diag(1/6)
    a = 0.0_dq
    a(1,1) = 6.0_dq; a(2,2) = 6.0_dq; a(3,3) = 6.0_dq
    b = 0.0_dq
    b(1,1) = 1.0_dq/6.0_dq; b(2,2) = 1.0_dq/6.0_dq; b(3,3) = 1.0_dq/6.0_dq
    z = 1                                  ! two H

    ! bonded pair: d = 0.74 A along z
    cart = 0.0_dq
    cart(3,2) = 0.74_dq
    CALL ef_model_hessian_full(2, nat, z, cart, a, b, H)
    ! SCHLEGEL: k = 18 Z_i Z_j / d^5, H-H d=0.74 -> 18/0.74^5 = 80.0
    CALL check(ABS(H(3,3) - (40.0_dq + 18.0_dq/0.74_dq**5)) < 1.0_dq, &
               'bonds: diagonal gains k along the bond')
    CALL check(ABS(H(1,4)) < 1.0D-12 .AND. ABS(H(1,1) - 40.0_dq) < 1.0D-10, &
               'bonds: no x-coupling')
    CALL check(ABS(H(1,4)) < 1.0D-12, 'bonds: cross xy zero')
    CALL check(ABS(H(1,4) - 0.0_dq) < 1.0D-12, 'bonds: (redundant)')

    ! far pair: d = 3.0 A > rcov sum + tol -> no bond term at all
    cart = 0.0_dq
    cart(3,2) = 3.0_dq
    CALL ef_model_hessian_full(2, nat, z, cart, a, b, H)
    CALL check(ABS(H(3,3) - 40.0_dq) < 1.0D-10, 'bonds: far pair diagonal plain')
    CALL check(ABS(H(1,4)) < 1.0D-12, 'bonds: far pair cross plain')

    ! PBC minimum image: atoms across the cell boundary, 0.2 A apart
    cart = 0.0_dq
    cart(3,1) = 0.1_dq
    cart(3,2) = 5.9_dq
    CALL ef_model_hessian_full(2, nat, z, cart, a, b, H)
    CALL check(ABS(H(3,3) - (40.0_dq + 18.0_dq/0.2_dq**5)) < 1.0_dq, &
               'bonds: PBC minimum-image detects the bond')
  END SUBROUTINE test_model_bonds

  SUBROUTINE test_extpen(nfail)
    INTEGER,INTENT(INOUT) :: nfail
    INTEGER,PARAMETER :: nat = 3
    INTEGER  :: z(nat), ia
    REAL(dq) :: cart(3,nat), m(nat), H(3*nat,3*nat), t(3*nat)
    REAL(dq) :: curv

    z = 1; m = 1.0_dq
    cart = 0.0_dq
    cart(1,2) = 0.8_dq; cart(1,3) = 1.6_dq
    H = 0.0_dq
    DO ia = 1, 3*nat
      H(ia,ia) = 5.0_dq
    END DO
    ! global +x translation (unit normalized)
    t = 0.0_dq
    DO ia = 1, nat
      t(3*(ia-1)+1) = 1.0_dq/SQRT(3.0_dq)
    END DO
    CALL ef_external_penalty(nat, cart, m, .FALSE., 1.0E4_dq, H)
    curv = DOT_PRODUCT(t, MATMUL(H, t))
    CALL check(ABS(curv - 1.0E4_dq - 5.0_dq) < 1.0_dq, &
               'extpen: translation stiffened to ~K')
    ! internal antisymmetric direction must be untouched
    t = 0.0_dq
    t(1) = 1.0_dq/SQRT(2.0_dq); t(4) = -1.0_dq/SQRT(2.0_dq)
    curv = DOT_PRODUCT(t, MATMUL(H, t))
    CALL check(ABS(curv - 5.0_dq) < 0.001_dq, &
               'extpen: internal direction untouched')
  END SUBROUTINE test_extpen

  SUBROUTINE test_freqs(nfail)
    INTEGER,INTENT(INOUT) :: nfail
    INTEGER,PARAMETER :: nat = 2
    REAL(dq) :: m(nat), H(3*nat,3*nat), fr(3*nat)
    INTEGER  :: nimag, info

    ! two uncoupled unit-mass atoms with +k and -k springs
    m = 1.0_dq
    H = 0.0_dq
    H(1,1) = 4.0_dq;  H(2,2) = 4.0_dq;  H(3,3) = 4.0_dq
    H(4,4) = -1.0_dq; H(5,5) = -1.0_dq; H(6,6) = -1.0_dq
    CALL ef_mass_weight_freqs(nat, m, H, fr, nimag, 10.0_dq, info)
    CALL check(info == 0, 'freqs: info==0')
    CALL check(ABS(fr(6) - 2.0_dq*EF_WAVENUMBER_CONST) < 1.0D-8, &
               'freqs: +4 -> +2*const cm^-1')
    CALL check(ABS(fr(1) + 1.0_dq*EF_WAVENUMBER_CONST) < 1.0D-8, &
               'freqs: -1 -> -1*const cm^-1')
    CALL check(nimag == 3, 'freqs: counts 3 imaginary')
  END SUBROUTINE test_freqs

  SUBROUTINE test_incar(nfail)
    INTEGER,INTENT(INOUT) :: nfail
    CHARACTER(*),PARAMETER :: f = 'test_incar.tmp'
    INTEGER :: unit, iv
    REAL(dq) :: rv
    LOGICAL :: lv, found
    CHARACTER(64) :: sv

    OPEN(NEWUNIT=unit, FILE=f, STATUS='REPLACE', ACTION='WRITE')
    WRITE(unit,'(A)') '# comment line'
    WRITE(unit,'(A)') 'EF_TS = .True.   ! enabled'
    WRITE(unit,'(A)') '  ef_maxstep = 0.25  ; trailing comment'
    WRITE(unit,'(A)') 'EF_MODE = track'
    WRITE(unit,'(A)') 'NSW = 200'
    WRITE(unit,'(A)') 'EF_NOCHECK = .FALSE.'
    CLOSE(unit)

    CALL ef_incar_logical(f, 'EF_TS', lv, found)
    CALL check(found .AND. lv, 'incar: logical true')
    CALL ef_incar_logical(f, 'EF_NOCHECK', lv, found)
    CALL check(found .AND. .NOT. lv, 'incar: logical false')
    CALL ef_incar_logical(f, 'EF_ABSENT', lv, found)
    CALL check(.NOT. found, 'incar: missing logical')
    CALL ef_incar_real(f, 'EF_MAXSTEP', rv, found)
    CALL check(found .AND. ABS(rv-0.25_dq) < 1.0D-14, 'incar: real')
    CALL ef_incar_int(f, 'NSW', iv, found)
    CALL check(found .AND. iv == 200, 'incar: int')
    CALL ef_incar_str(f, 'EF_MODE', sv, found)
    CALL check(found .AND. TRIM(sv) == 'TRACK', 'incar: string uppercased')
    ! bad value -> treated as absent
    OPEN(NEWUNIT=unit, FILE=f, STATUS='REPLACE', ACTION='WRITE')
    WRITE(unit,'(A)') 'EF_MAXSTEP = banana'
    CLOSE(unit)
    CALL ef_incar_real(f, 'EF_MAXSTEP', rv, found)
    CALL check(.NOT. found, 'incar: unparsable value rejected')
  END SUBROUTINE test_incar

END PROGRAM test_ef_core
