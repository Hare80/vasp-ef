#include "symbol.inc"
!=====================================================================
! ef_ts.F -- Hessian-based eigenvector-following transition-state
! optimizer for VASP (P-RFO step, Bofill/Powell updates, mode tracking).
!
! This is the whole optimizer.  It is hooked into VASP with two lines
! in src/chain.F (see patches/ and install.sh in the VASP-EF
! distribution):
!
!     USE ef_ts            ! in the header of MODULE chain
!     ...
!     IF (ef_ts_step(nions,posion,toten,force,a,b,iu6)) RETURN
!                                          ! first line of chain_force
!
! Activation protocol (same convention as the VTST optimizers):
!     INCAR:  IBRION = 3, POTIM = 0, EF_TS = .TRUE.
!
! INCAR tags (all optional unless noted):
!     EF_TS          logical, .TRUE. to activate            (required)
!     EF_MAXSTEP     trust radius / max step length, A      (default 0.2)
!     EF_TRUST       FIXED (default) | ADAPTIVE
!     EF_MODE        LOWEST (default) | TRACK | GUESS
!     EF_UPDATE      BOFILL (default) | POWELL
!     EF_NOCHECK     declare convergence on forces alone (Gaussian
!                    NoEigenTest style); default .FALSE.
!     EF_HSCALE      scale of the model Hessian             (default 1.0)
!     EF_DIAG        FULL only (P-RFO needs the full spectrum); other
!                    values are ignored with a note
!     EF_READ_HESS   .TRUE.: read initial Hessian from EFHESSIAN
!     EF_HESS_NWRITE write EFHESSIAN every N steps (default 0: only
!                    at convergence)
!     EF_PRINT       output verbosity 0|1|2                 (default 1)
!
! Files in the run directory:
!     MODECAR       read at startup when EF_MODE=GUESS (vtst-style
!                   format, Cartesian components, 3 numbers per atom)
!     NEWMODECAR    followed mode after every step (vtst convention)
!     EFDAT         one line per step with diagnostics
!     EFHESSIAN     Cartesian Hessian snapshot; read back with
!                   EF_READ_HESS
!
! Limitations (v1): fixed unit cell (ISIF<=2); Selective Dynamics is
! not treated explicitly (same as the VTST optimizers); all MPI ranks
! redundantly perform the (cheap) optimizer arithmetic on the full
! force array, exactly like VTST does.
!
! SPDX-License-Identifier: MIT
!=====================================================================
  MODULE ef_ts
    USE prec
    USE ef_core
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: ef_ts_step

    ! ----- optimizer state ------------------------------------------
    LOGICAL,SAVE :: initialized=.FALSE., active=.FALSE.
    INTEGER,SAVE :: nion=0, nstep=0
    INTEGER,SAVE :: iu6_saved=6                 ! stdout unit from chain_force
    REAL(q),ALLOCATABLE,SAVE :: H(:,:)          ! (3n,3n) Cartesian Hessian
    REAL(q),ALLOCATABLE,SAVE :: g_old(:,:)      ! (3,nions) previous gradient
    REAL(q),ALLOCATABLE,SAVE :: s_try(:)        ! (3n) step of the pending proposal
    REAL(q),ALLOCATABLE,SAVE :: pos_good(:,:), pos_try(:,:)
    REAL(q),ALLOCATABLE,SAVE :: mode_prev(:)    ! (3n) previously followed mode
    REAL(q),ALLOCATABLE,SAVE :: amass(:)        ! (nions) atomic masses, amu
    INTEGER,ALLOCATABLE,SAVE :: znum(:)         ! (nions) atomic numbers
    LOGICAL,SAVE :: have_prev=.FALSE., have_gold=.FALSE.
    LOGICAL,SAVE :: have_eold=.FALSE.
    LOGICAL,SAVE :: have_good=.FALSE., have_try=.FALSE.
    REAL(q),SAVE :: e_old=0.0_q, e_good=0.0_q, emax=0.2_q

    ! ----- INCAR parameters ------------------------------------------
    REAL(q),SAVE :: maxstep=0.2_q, trust=0.2_q, trust0=0.2_q
    REAL(q),SAVE :: hscale=1.0_q
    REAL(q),SAVE :: ediffg=0.0_q
    LOGICAL,SAVE :: have_ediffg=.FALSE.
    INTEGER,SAVE :: hess_nwrite=0, iprint=1
    LOGICAL,SAVE :: ladaptive=.FALSE., lnocheck=.FALSE.
    LOGICAL,SAVE :: lreadhess=.FALSE.
    LOGICAL,SAVE :: climbing=.FALSE.            ! one-way: min phase -> climb
    REAL(q),SAVE :: fswitch=0.5_q               ! eV/A: start climbing below

    ! ----- FD curvature refinement along the TS mode -------------------
    LOGICAL,SAVE :: lfdrefine=.FALSE.
    LOGICAL,SAVE :: lrotfix=.FALSE.             ! project molecular rotations
    INTEGER,SAVE :: fdref_k=10, fdpending=0, steps_since_fd=0
    REAL(q),PARAMETER :: EXTPEN_K = 1.0E4_q     ! external-DOF penalty
    REAL(q),SAVE :: fddelta=0.04_q, fdproj_plus=0.0_q
    REAL(q),ALLOCATABLE,SAVE :: fdv(:)          ! (3n) refinement direction

    ! ----- bond-aware model Hessian --------------------------------------
    LOGICAL,SAVE :: lbonds=.TRUE.
    REAL(q),SAVE :: bondk=80.0_q, btol=0.45_q

    ! ----- CALCFC: finite-difference Hessian at startup ------------------
    ! EF_FCHESS_N > 0 additionally REBUILDS the exact Hessian every N
    ! steps at the current position (Gaussian RecalcFC analog) -- resets
    ! Bofill drift on long or wandering climbs.
    LOGICAL,SAVE :: lcalcfc=.FALSE., fcdone=.FALSE.
    INTEGER,SAVE :: fcstage=0, fck=0, nact=0, fchess_n=0, steps_since_fchess=0
    INTEGER,ALLOCATABLE,SAVE :: actlist(:)      ! active atoms for CALCFC
    REAL(q),SAVE :: fcdelta=0.02_q
    REAL(q),ALLOCATABLE,SAVE :: fcFplus_vec(:), fcFminus_vec(:)

    ! ----- mode-curvature least-squares history (Sella-style) ----------
    INTEGER,PARAMETER :: HISTM = 4
    REAL(q),ALLOCATABLE,SAVE :: sh_hist(:,:)     ! (HISTM,3n) accepted steps
    REAL(q),ALLOCATABLE,SAVE :: yh_hist(:,:)     ! (HISTM,3n) gradient diffs
    INTEGER,SAVE :: hist_n=0, hist_p=0
      INTEGER,SAVE :: iupd=1                      ! 1=Bofill, 2=Powell
      CHARACTER(1),SAVE :: cmode='L'              ! L=LOWEST, T=TRACK, G=GUESS

      CHARACTER(512),SAVE :: line_tmp             ! scratch for file parsing
      INTEGER,SAVE :: eq_index                    ! scratch for file parsing

    REAL(q),PARAMETER :: ZERO_THRESH_CM = 50.0_q ! |nu| below = zero mode
    REAL(q),SAVE :: guess_k = 1.0_q              ! seeded curvature (EF_GUESS_K)

  CONTAINS

!=====================================================================
! main entry: called from chain_force on every ionic step.
! Returns .TRUE. if EF_TS is active and the positions were updated.
!=====================================================================
    FUNCTION ef_ts_step(ni,posion,toten,force,a,b,iu6) RESULT(handled)
      INTEGER,INTENT(IN)      :: ni, iu6
      REAL(q),INTENT(INOUT)   :: posion(3,ni)
      REAL(q),INTENT(IN)      :: force(3,ni), toten, a(3,3), b(3,3)
      LOGICAL                 :: handled

      INTEGER  :: n, i, j, info, imode, ierr, nimag
      LOGICAL  :: lapplied
      REAL(q)  :: sv, yv
      REAL(q),ALLOCATABLE :: vg(:), cart(:,:), Hwork(:,:)
      REAL(q)  :: g(3,ni), s3(3), sflat(3*ni), y(3*ni)
      REAL(q),ALLOCATABLE :: eval(:), evec(:,:)
      REAL(q),ALLOCATABLE :: freqs(:)
      REAL(q)  :: fmax, de, kappa, nu1, slen

      handled = .FALSE.
      IF (.NOT. initialized) THEN
        initialized = .TRUE.
        CALL ef_init(ni, iu6)
      END IF
      IF (.NOT. active) RETURN
      handled = .TRUE.
      n = 3*ni

      ! Cartesian gradient (g = -F); forces are Cartesian eV/A
      g = -force
      fmax = MAXVAL(ABS(force))

      ! --- FD curvature probe state machine ------------------------------
      ! These calls evaluate forces at DISPLACED points, so they must be
      ! handled before the accept/reject logic.
      IF (fdpending == 1) THEN
        ! forces at +delta: record projection, move to -delta
        fdproj_plus = DOT_PRODUCT(RESHAPE(force,(/n/)), fdv)
        DO i = 1, ni
          s3 = fdv(3*(i-1)+1:3*i)
          posion(:,i) = pos_good(:,i) - MATMUL(s3, b)*fddelta
        END DO
        fdpending = 2
        RETURN
      END IF
      IF (fdpending == 2) THEN
        ! forces at -delta: central-difference curvature, lock into H,
        ! restore the accepted position, request one refresh evaluation
        kappa = (DOT_PRODUCT(RESHAPE(force,(/n/)), fdv) - fdproj_plus) &
                / (2.0_q*fddelta)
        CALL ef_seed_mode(n, H, fdv, kappa)
        posion = pos_good
        fdpending = 3
        steps_since_fd = 0
        IF (iu6>=0 .AND. iprint>=1) WRITE(iu6,'(A,F10.3,A)') &
          ' EF: FD probe: curvature along TS mode =',kappa,' eV/A^2'
        RETURN
      END IF
      IF (fdpending == 3) THEN
        ! forces are back at the accepted position: continue normally
        fdpending = 0
      END IF

      ! --- CALCFC: build the exact Hessian by central differences -------
      ! Consumes 6*nact NSW iterations at the start; EF_ACTIVE selects
      ! the active atoms (default: all).
      IF (lcalcfc .AND. .NOT. fcdone) THEN
        IF (fcstage == 1) THEN
          ! forces at +delta of column fck
          fcFplus_vec = RESHAPE(force, (/n/))
          posion = pos_good
          DO i = 1, 3
            posion(i, actlist((fck-1)/3+1)) = pos_good(i, actlist((fck-1)/3+1)) &
              - fcdelta*b(MOD(fck-1,3)+1, i)
          END DO
          fcstage = 2
          RETURN
        ELSE IF (fcstage == 2) THEN
          ! forces at -delta: assemble the column
          fcFminus_vec = RESHAPE(force, (/n/))
          j = 3*(actlist((fck-1)/3+1)-1) + MOD(fck-1,3) + 1
          DO i = 1, n
            H(i,j) = (fcFminus_vec(i) - fcFplus_vec(i)) / (2.0_q*fcdelta)
            H(j,i) = H(i,j)
          END DO
          fck = fck + 1
          IF (fck > 3*nact) THEN
            fcdone = .TRUE.
            fcstage = 0
            posion = pos_good
            IF (iu6>=0 .AND. iprint>=1) WRITE(iu6,'(A)') &
              ' EF: CALCFC Hessian complete'
            ! fall through: normal processing starts below
          ELSE
            posion = pos_good
            DO i = 1, 3
              posion(i, actlist((fck-1)/3+1)) = pos_good(i, actlist((fck-1)/3+1)) &
                + fcdelta*b(MOD(fck-1,3)+1, i)
            END DO
            fcstage = 1
            RETURN
          END IF
        ELSE
          ! fcstage == 0: start column 1 (+delta displacement).
          ! First use of pos_good: anchor it at the incoming geometry.
          pos_good = posion
          e_good  = toten
          fck = 1
          fcstage = 1
          posion = pos_good
          DO i = 1, 3
            posion(i, actlist((fck-1)/3+1)) = pos_good(i, actlist((fck-1)/3+1)) &
              + fcdelta*b(MOD(fck-1,3)+1, i)
          END DO
          RETURN
        END IF
      END IF

      ! --- accept/reject the previous proposal (trust region) ----------
      IF (have_try) THEN
        have_try = .FALSE.
        IF (toten > e_good + emax) THEN
          ! rejected: energy rose too much -- restore and retry smaller.
          ! Hessian/g_old are untouched (updates are deferred to
          ! acceptance), so the retry starts from a clean state.
          posion = pos_good
          IF (trust <= 0.020001_q*trust0) THEN
            ! already at the trust floor and still rising: the model
            ! Hessian is useless here -- rebuild it and re-enter the
            ! minimization phase instead of dead-locking.
            CALL build_model_hessian(ni, pos_good, a, b)
            trust = trust0
            have_prev = .FALSE.
            climbing = .FALSE.
            IF (iu6>=0 .AND. iprint>=1) WRITE(iu6,'(A)') &
              ' EF: trust floor reached -- Hessian reset, back to MIN phase'
          ELSE
            trust = MAX(0.1_q*trust0, 0.5_q*trust)
          END IF
          IF (iu6>=0 .AND. iprint>=1) WRITE(iu6,'(A,F8.3,A,F6.3)') &
            ' EF: step rejected (dE =',toten-e_good,' eV), trust -> ',trust
          RETURN
        END IF
        ! accepted: apply the deferred Hessian update with the (s,y) pair
        ! that produced this position
        y = RESHAPE(g,(/n/)) - RESHAPE(g_old,(/n/))
        IF (SQRT(DOT_PRODUCT(y,y)) <= 1.0E2_q) THEN
          IF (iupd == 1) THEN
            CALL ef_bofill_update(n, H, s_try, y, lapplied)
          ELSE
            CALL ef_powell_update(n, H, s_try, y, lapplied)
          END IF
          ! (mode-curvature calibration is handled by the periodic FD
          ! probes, EF_FDREFINE, which use dedicated displaced points)
        ELSE IF (iu6>=0 .AND. iprint>=1) THEN
          WRITE(iu6,'(A)') ' EF: NOTE: ||y|| too large, Hessian update skipped'
        END IF
        pos_good = posion
        e_good = toten
      END IF

      nstep = nstep + 1

      ALLOCATE(eval(n), evec(n,n), freqs(n))

      ! --- initial Hessian on the first step ---------------------------
      IF (nstep == 1 .AND. .NOT. fcdone) THEN
        IF (lreadhess) THEN
          CALL read_efhessian(ni, posion, a, b)
        ELSE
          CALL build_model_hessian(ni, posion, a, b)
        END IF
        IF (iu6>=0 .AND. iprint>=1) THEN
          IF (lreadhess) THEN
            WRITE(iu6,'(A)') ' EF: initial Hessian read from EFHESSIAN'
          ELSE
            WRITE(iu6,'(A)') ' EF: initial Hessian = element model (EF_HSCALE)'
          END IF
        END IF
        ! EF_MODE=GUESS: seed the Hessian along the MODECAR direction so
        ! that it becomes an actual eigenmode (soft, k = 1 eV/A^2).  The
        ! diagonal model Hessian alone has atom-local modes with no
        ! physical meaning; the guess gives the climb a real direction.
        IF (cmode == 'G') THEN
          ALLOCATE(vg(n)); vg = 0.0_q
          IF (read_modecar('MODECAR', vg)) THEN
            CALL ef_seed_mode(n, H, vg, guess_k)
            IF (iu6>=0 .AND. iprint>=1) WRITE(iu6,'(A,F7.3)') &
              ' EF: Hessian seeded along MODECAR direction, k = ',guess_k
          ELSE
            IF (iu6>=0) WRITE(iu6,'(A)') &
              ' EF: WARNING: EF_MODE=GUESS but no MODECAR; plain model used'
          END IF
          DEALLOCATE(vg)
        END IF
      END IF

      ! --- eigendecomposition -------------------------------------------
      ! The spectrum is computed on a copy with the external (rigid)
      ! directions penalized: their noisy curvatures must never be
      ! amplified by the RFO shift or selected as the TS mode.
      ALLOCATE(Hwork(n,n), cart(3,ni))
      DO i = 1, ni
        cart(:,i) = MATMUL(posion(:,i), a)
      END DO
      Hwork = H
      CALL ef_external_penalty(ni, cart, amass, lrotfix, EXTPEN_K, Hwork)
      CALL ef_sym_eigh(n, Hwork, eval, evec, info)
      IF (info /= 0) THEN
        IF (iu6>=0) WRITE(iu6,'(A,I0)') ' EF: FATAL: eigen-solver failed, info=', info
        STOP
      END IF

      ! --- mode selection -------------------------------------------------
      imode = select_mode(n, eval, evec)
      kappa = eval(imode)
      IF (ALLOCATED(mode_prev)) DEALLOCATE(mode_prev)
      ALLOCATE(mode_prev(n)); mode_prev = evec(:,imode)
      have_prev = .TRUE.

      ! --- two-phase strategy ---------------------------------------------
      ! While forces are large the geometry is not in the quasi-Newton
      ! regime and an uphill step would fight stiff directions: first
      ! minimize (pure RFO step).  Once the residual force drops below
      ! EF_FSWITCH -- or the followed mode has turned negative -- latch
      ! into the climb phase permanently.
      IF (.NOT. climbing .AND. (kappa < 0.0_q .OR. fmax < fswitch)) THEN
        climbing = .TRUE.
        IF (iu6>=0 .AND. iprint>=1) THEN
          IF (kappa < 0.0_q) THEN
            WRITE(iu6,'(A)') ' EF: negative curvature detected -- start climbing'
          ELSE
            WRITE(iu6,'(A,F6.3)') ' EF: forces small enough -- start climbing, Fmax = ',fmax
          END IF
        END IF
      END IF

      ! --- periodic exact-Hessian rebuild (RecalcFC) -----------------------
      IF (lcalcfc .AND. fchess_n > 0 .AND. .NOT. fcdone) THEN
        IF (steps_since_fchess >= fchess_n .AND. fdpending == 0) THEN
          fcdone = .FALSE.
          fcstage = 0
          fck = 1
          steps_since_fchess = 0
          have_try = .FALSE.
          pos_good = posion
          posion = pos_good
          DO i = 1, 3
            posion(i, actlist((fck-1)/3+1)) = pos_good(i, actlist((fck-1)/3+1)) &
              + fcdelta*b(MOD(fck-1,3)+1, i)
          END DO
          fcstage = 1
          IF (iu6>=0 .AND. iprint>=1) WRITE(iu6,'(A)') &
            ' EF: RecalcFC: rebuilding exact Hessian at current position'
          RETURN
        END IF
      END IF

      ! --- FD curvature refinement trigger ---------------------------------
      ! Replaces this iteration's step proposal: displace +delta along
      ! the followed mode; the two following iterations evaluate -delta
      ! and the refresh.  Costs 3 NSW iterations, no main.F patch.
      steps_since_fd = steps_since_fd + 1
      IF (lfdrefine .AND. climbing .AND. have_prev .AND. &
          fdpending == 0 .AND. steps_since_fd > fdref_k) THEN
        fdv = mode_prev
        DO i = 1, ni
          s3 = fdv(3*(i-1)+1:3*i)
          posion(:,i) = pos_good(:,i) + MATMUL(s3, b)*fddelta
        END DO
        fdpending = 1
        steps_since_fd = 0
        IF (iu6>=0 .AND. iprint>=1) WRITE(iu6,'(A)') &
          ' EF: FD probe: displacing +delta along TS mode'
        RETURN
      END IF

      ! --- P-RFO step -------------------------------------------------------
      CALL ef_prfo_step(n, eval, evec, RESHAPE(g,(/n/)), imode, climbing, &
                        trust, sflat, slen, ierr)
      IF (ierr /= 0) THEN
        IF (iu6>=0) WRITE(iu6,'(A)') ' EF: zero force, nothing to do'
        RETURN
      END IF

      ! apply the Cartesian step in direct coordinates and record the
      ! proposal.  The Hessian update is DEFERRED until the proposal is
      ! accepted (next entry), so rejected steps never poison H.
      pos_good = posion                        ! reference = accepted position
      e_good  = toten
      s_try   = sflat
      DO i = 1, ni
        s3 = sflat(3*(i-1)+1:3*i)             ! Cartesian step of atom i
        posion(:,i) = posion(:,i) + MATMUL(s3, b)
      END DO
      pos_try = posion
      have_try = .TRUE.
      g_old = g
      have_gold = .TRUE.
      steps_since_fchess = steps_since_fchess + 1

      ! --- sanity guard on the accumulated Hessian -------------------------
      IF (.NOT. hessian_is_sane(n)) THEN
        IF (iu6>=0) WRITE(iu6,'(A)') &
          ' EF: WARNING: Hessian corrupted (NaN/overflow); reset to model'
        CALL build_model_hessian(ni, pos_good, a, b)
        trust = trust0
        have_prev = .FALSE.
        climbing = .FALSE.
      END IF

      ! --- adaptive trust radius ---------------------------------------------
      de = 0.0_q
      IF (have_eold) de = toten - e_old
      IF (ladaptive .AND. have_eold .AND. .NOT. climbing) THEN
        ! minimization phase: rising energy shrinks the radius
        IF (de < 0.0_q .AND. slen < 0.99_q*trust) THEN
          trust = MIN(1.5_q*trust, 4.0_q*trust0)
        ELSE IF (de >= 0.01_q) THEN
          trust = MAX(0.5_q*trust, 0.02_q*trust0)
        END IF
      ELSE IF (climbing) THEN
        ! climb phase: energy RISES toward the barrier.  Endgame ratchet:
        ! as Fmax approaches the convergence criterion, cap the radius in
        ! two stages so the trajectory settles on the (possibly flat)
        ! saddle top instead of bouncing across it.
        IF (have_ediffg) THEN
          IF (fmax < 3.0_q*ABS(ediffg)) THEN
            trust = MIN(trust, 0.05_q*trust0)
          ELSE IF (fmax < 10.0_q*ABS(ediffg)) THEN
            trust = MIN(trust, 0.10_q*trust0)
          ELSE
            IF (slen < 0.99_q*trust) trust = MIN(1.5_q*trust, 4.0_q*trust0)
          END IF
        ELSE
          IF (slen < 0.99_q*trust) trust = MIN(1.5_q*trust, 4.0_q*trust0)
        END IF
      END IF
      e_old = toten
      have_eold = .TRUE.

      ! --- diagnostics ----------------------------------------------------------
      ! spectrum on the PENALIZED copy: external (rigid) modes sit at
      ! ~+K and are excluded from the imaginary-mode count
      CALL ef_mass_weight_freqs(ni, amass, Hwork, freqs, nimag, ZERO_THRESH_CM, info)
      nu1 = freqs(1)
      IF (iu6>=0 .AND. iprint>=1) THEN
        WRITE(iu6,'(A,I5,A,F14.6,A,F10.5,A,F10.4,A,I0,A,F10.2,A,I0)') &
          ' EF: step ',nstep,' E ',toten,' fmax ',fmax, &
          ' kappa ',kappa,' mode ',imode,' nu1 ',nu1,' nimag ',nimag
        IF (iprint>=2) WRITE(iu6,'(A,F8.4,A,ES10.3,A,A)') &
          ' EF: |step| ',slen,' dE ',de,' phase ', &
          MERGE('CLIMB','MIN  ',climbing)
      END IF
      CALL write_efdat(toten, de, fmax, kappa, imode, nu1, nimag, slen)

      ! --- snapshots and convergence report --------------------------------
      IF (hess_nwrite > 0 .AND. MOD(nstep, hess_nwrite)==0) CALL write_efhessian(n)
      IF (have_ediffg .AND. fmax < ABS(ediffg)) THEN
        CALL write_efhessian(n)
        IF (nimag == 1 .OR. lnocheck) THEN
          IF (iu6>=0) WRITE(iu6,'(A)') &
            ' EF: CONVERGED (force criterion + single imaginary frequency)'
        ELSE IF (iu6>=0) THEN
          WRITE(iu6,'(A,I0,A)') ' EF: WARNING: force converged but ', &
               nimag,' imaginary modes found.'
          WRITE(iu6,'(A)') &
            '     Consider EF_MODE=TRACK/GUESS, a MODECAR from a dimer run, or EF_NOCHECK=.TRUE.'
        END IF
      END IF

      DEALLOCATE(eval, evec, freqs)
    END FUNCTION ef_ts_step

!=====================================================================
! initialization: read INCAR, set parameters, parse species info
!=====================================================================
    SUBROUTINE ef_init(ni, iu6)
      INTEGER,INTENT(IN) :: ni, iu6
      CHARACTER(64)      :: sval
      LOGICAL            :: found
      REAL(q)            :: rval
      INTEGER            :: ival, i

      nion = ni
      iu6_saved = iu6

      CALL ef_incar_logical('INCAR','EF_TS',active,found)
      IF (.NOT. active) RETURN

      ! validate the vtst-style protocol
      CALL ef_incar_int('INCAR','IBRION',ival,found)
      IF (.NOT. found .OR. ival /= 3) THEN
        IF (iu6>=0) WRITE(iu6,'(A)') ' EF: FATAL: EF_TS requires IBRION=3'
        STOP
      END IF
      rval = 1.0_q
      CALL ef_incar_real('INCAR','POTIM',rval,found)
      IF (.NOT. found .OR. ABS(rval) > 1.0E-12_q) THEN
        IF (iu6>=0) WRITE(iu6,'(A)') ' EF: FATAL: EF_TS requires POTIM=0'
        STOP
      END IF

      CALL ef_incar_real('INCAR','EF_MAXSTEP',maxstep,found)
      IF (.NOT. found) maxstep = 0.2_q
      trust = maxstep; trust0 = maxstep

      CALL ef_incar_str('INCAR','EF_TRUST',sval,found)
      IF (found .AND. TRIM(sval)=='ADAPTIVE') ladaptive = .TRUE.

      CALL ef_incar_str('INCAR','EF_MODE',sval,found)
      IF (found) cmode = sval(1:1)

      CALL ef_incar_str('INCAR','EF_UPDATE',sval,found)
      IF (found .AND. TRIM(sval)=='POWELL') iupd = 2

      CALL ef_incar_logical('INCAR','EF_NOCHECK',lnocheck,found)
      IF (.NOT. found) lnocheck = .FALSE.

      CALL ef_incar_real('INCAR','EF_HSCALE',hscale,found)
      IF (.NOT. found) hscale = 1.0_q

      CALL ef_incar_real('INCAR','EF_FSWITCH',fswitch,found)
      IF (.NOT. found) fswitch = 0.5_q

      CALL ef_incar_real('INCAR','EF_GUESS_K',guess_k,found)
      IF (.NOT. found) guess_k = 1.0_q

      CALL ef_incar_real('INCAR','EF_EMAX',emax,found)
      IF (.NOT. found) emax = 0.2_q

      CALL ef_incar_logical('INCAR','EF_ROTFIX',lrotfix,found)
      IF (.NOT. found) lrotfix = .FALSE.       ! enable for gas-phase molecules

      CALL ef_incar_logical('INCAR','EF_FDREFINE',lfdrefine,found)
      IF (.NOT. found) lfdrefine = .FALSE.
      CALL ef_incar_int('INCAR','EF_FDREF_K',fdref_k,found)
      IF (.NOT. found .OR. fdref_k < 1) fdref_k = 10
      CALL ef_incar_real('INCAR','EF_FDDELTA',fddelta,found)
      IF (.NOT. found .OR. fddelta <= 0.0_q) fddelta = 0.04_q

      CALL ef_incar_logical('INCAR','EF_BONDS',lbonds,found)
      IF (.NOT. found) lbonds = .TRUE.
      CALL ef_incar_real('INCAR','EF_BOND_K',bondk,found)
      IF (.NOT. found) bondk = 80.0_q
      CALL ef_incar_real('INCAR','EF_BOND_TOL',btol,found)
      IF (.NOT. found) btol = 0.45_q

      CALL ef_incar_logical('INCAR','EF_CALCFC',lcalcfc,found)
      IF (.NOT. found) lcalcfc = .FALSE.
      CALL ef_incar_int('INCAR','EF_FCHESS_N',fchess_n,found)
      IF (.NOT. found .OR. fchess_n < 0) fchess_n = 0
      CALL ef_incar_real('INCAR','EF_FCDELTA',fcdelta,found)
      IF (.NOT. found .OR. fcdelta <= 0.0_q) fcdelta = 0.02_q
      ALLOCATE(actlist(ni))
      CALL ef_incar_int_list('INCAR','EF_ACTIVE',ni,actlist,nact,found)
      IF (.NOT. found .OR. nact <= 0) THEN
        DO i = 1, ni
          actlist(i) = i
        END DO
        nact = ni
      END IF
      ALLOCATE(fcFplus_vec(3*ni), fcFminus_vec(3*ni))
      ALLOCATE(sh_hist(HISTM,3*ni), yh_hist(HISTM,3*ni))
      sh_hist = 0.0_q; yh_hist = 0.0_q
      IF (lcalcfc .AND. iu6>=0 .AND. iprint>=1) THEN
        WRITE(iu6,'(A,I0,A,F6.3,A,I0,A)') &
          ' EF: CALCFC active-region Hessian: ',nact,' atoms, delta = ', &
          fcdelta,' A (consumes ',6*nact,' NSW iterations)'
      END IF

      ALLOCATE(fdv(3*ni)); fdv = 0.0_q
      IF (lfdrefine .AND. iu6>=0 .AND. iprint>=1) &
        WRITE(iu6,'(A,I0,A,F5.2,A)') &
        ' EF: NOTE: FD refinement every ',fdref_k, &
        ' steps, delta = ',fddelta, &
        ' A; each probe consumes 3 extra NSW iterations'

      CALL ef_incar_logical('INCAR','EF_READ_HESS',lreadhess,found)
      IF (.NOT. found) lreadhess = .FALSE.

      CALL ef_incar_int('INCAR','EF_HESS_NWRITE',hess_nwrite,found)
      IF (.NOT. found) hess_nwrite = 0

      CALL ef_incar_int('INCAR','EF_PRINT',iprint,found)
      IF (.NOT. found) iprint = 1

      CALL ef_incar_real('INCAR','EDIFFG',ediffg,have_ediffg)

      ! EF_DIAG accepted for forward compatibility; P-RFO needs the
      ! full spectrum, so only FULL is implemented.
      CALL ef_incar_str('INCAR','EF_DIAG',sval,found)
      IF (found .AND. TRIM(sval)/='FULL' .AND. iu6>=0) &
        WRITE(iu6,'(A)') ' EF: NOTE: EF_DIAG other than FULL ignored (P-RFO needs all modes)'

      ! species info for the model Hessian and mass weighting
      ALLOCATE(amass(ni), znum(ni))
      CALL parse_species(ni, znum, amass)

      ALLOCATE(H(3*ni,3*ni), g_old(3,ni))
      ALLOCATE(pos_good(3,ni), pos_try(3,ni))
      ALLOCATE(s_try(3*ni))
      H = 0.0_q; g_old = 0.0_q

      IF (iu6>=0 .AND. iprint>=1) THEN
        WRITE(iu6,'(A)') ' '
        WRITE(iu6,'(A)') ' EF: VASP-EF eigenvector-following TS optimizer active'
        WRITE(iu6,'(A,F7.3,A,A,A,A,A,A)') ' EF: maxstep ',maxstep, &
             ' trust ',TRIM(trust_name()),' mode ',TRIM(mode_name()), &
             ' update ',TRIM(upd_name())
      END IF
    END SUBROUTINE ef_init

    FUNCTION trust_name() RESULT(s)
      CHARACTER(16) :: s
      IF (ladaptive) THEN; s='ADAPTIVE'; ELSE; s='FIXED'; END IF
    END FUNCTION trust_name

    ! NaN / overflow check on the approximate Hessian
    FUNCTION hessian_is_sane(n) RESULT(ok)
      INTEGER,INTENT(IN) :: n
      LOGICAL :: ok
      INTEGER :: i
      ok = .FALSE.
      DO i = 1, n
        IF (ANY(H(i,1:n) /= H(i,1:n))) RETURN            ! NaN
        IF (MAXVAL(ABS(H(i,1:n))) > 1.0E8_q) RETURN      ! overflow
      END DO
      ok = .TRUE.
    END FUNCTION hessian_is_sane

    FUNCTION mode_name() RESULT(s)
      CHARACTER(16) :: s
      IF (cmode=='T') THEN
        s='TRACK'
      ELSE IF (cmode=='G') THEN
        s='GUESS'
      ELSE
        s='LOWEST'
      END IF
    END FUNCTION mode_name

    FUNCTION upd_name() RESULT(s)
      CHARACTER(16) :: s
      IF (iupd==1) THEN; s='BOFILL'; ELSE; s='POWELL'; END IF
    END FUNCTION upd_name

!=====================================================================
! species parsing: POTCAR TITEL lines give the element of each type
! (in POSCAR order); POSCAR line 6/7 give the counts.
!=====================================================================
    SUBROUTINE parse_species(ni, z, mass)
      INTEGER,INTENT(IN)  :: ni
      INTEGER,INTENT(OUT) :: z(ni)
      REAL(q),INTENT(OUT) :: mass(ni)

      CHARACTER(2)   :: sym(120)
      INTEGER        :: ntyp, counts(120), ia, it, ios, unit, nt, n2
      LOGICAL        :: v5format, isal

      z = 0; mass = 40.0_q
      CALL read_potcar_titel(sym, ntyp)
      IF (ntyp == 0) THEN
        IF (iu6_loc()>=0) WRITE(iu6_loc(),'(A)') &
          ' EF: WARNING: could not parse POTCAR; default masses/stiffness used'
        RETURN
      END IF

      OPEN(NEWUNIT=unit, FILE='POSCAR', STATUS='OLD', ACTION='READ', IOSTAT=ios)
      IF (ios /= 0) THEN
        IF (iu6_loc()>=0) WRITE(iu6_loc(),'(A)') &
          ' EF: WARNING: no POSCAR; default masses/stiffness used'
        RETURN
      END IF
      DO it = 1, 5
        READ(unit,'(A)',IOSTAT=ios) line_tmp   ! comment, scale, lattice
        IF (ios /= 0) EXIT
      END DO
      READ(unit,'(A)',IOSTAT=ios) line_tmp     ! line 6: symbols (v5) or counts (v4)
      CALL parse_counts(line_tmp, counts, nt, v5format)
      IF (v5format) THEN
        READ(unit,'(A)',IOSTAT=ios) line_tmp   ! line 7: counts
        CALL parse_counts(line_tmp, counts, nt, isal)
      END IF
      CLOSE(unit)
      IF (nt /= ntyp .AND. iu6_loc()>=0) &
        WRITE(iu6_loc(),'(A)') ' EF: WARNING: POTCAR/POSCAR species count mismatch'

      ia = 0
      DO it = 1, MIN(nt, ntyp, 120)
        DO n2 = 1, counts(it)
          ia = ia + 1
          IF (ia > ni) EXIT
          z(ia)    = symbol_z(sym(it))
          mass(ia) = symbol_mass(sym(it))
        END DO
      END DO
      DO ia = 1, ni                            ! fill any remainder defensively
        IF (z(ia) == 0) THEN
          z(ia) = 6; mass(ia) = 12.0_q
        END IF
      END DO
    END SUBROUTINE parse_species

    FUNCTION iu6_loc() RESULT(iu)
      INTEGER :: iu
      iu = iu6_saved
    END FUNCTION iu6_loc

    SUBROUTINE read_potcar_titel(sym, ntyp)
      CHARACTER(2),INTENT(OUT) :: sym(120)
      INTEGER,INTENT(OUT)      :: ntyp
      CHARACTER(512)           :: line, rhs
      INTEGER                  :: unit, ios, n
      LOGICAL                  :: exists
      ntyp = 0; sym = '  '
      INQUIRE(FILE='POTCAR', EXIST=exists)
      IF (.NOT. exists) RETURN
      OPEN(NEWUNIT=unit, FILE='POTCAR', STATUS='OLD', ACTION='READ', IOSTAT=ios)
      IF (ios /= 0) RETURN
      DO
        READ(unit,'(A)',IOSTAT=ios) line
        IF (ios /= 0) EXIT
        IF (INDEX(line, 'TITEL') == 0) CYCLE
        eq_index = INDEX(line, '=')
        IF (eq_index == 0) CYCLE
        rhs = ADJUSTL(line(eq_index+1:))
        ! rhs = "<library> <symbol> <date>" -> drop token 1, keep token 2
        n = SCAN(rhs, ' '//CHAR(9))
        IF (n == 0) CYCLE
        rhs = ADJUSTL(rhs(n:))
        n = SCAN(rhs, ' '//CHAR(9))
        IF (n > 1) rhs = rhs(1:n-1)
        ntyp = ntyp + 1
        IF (ntyp <= 120) sym(ntyp) = upcase2(TRIM(ADJUSTL(rhs)))
      END DO
      CLOSE(unit)
    END SUBROUTINE read_potcar_titel

    SUBROUTINE parse_counts(line, counts, nt, is_alpha)
      CHARACTER(*),INTENT(IN)  :: line
      INTEGER,INTENT(OUT)      :: counts(:)
      INTEGER,INTENT(OUT)      :: nt
      LOGICAL,INTENT(OUT)      :: is_alpha
      INTEGER :: ios, pos, v
      CHARACTER(64) :: tok
      counts = 0; nt = 0; is_alpha = .FALSE.
      pos = 1
      DO
        CALL next_token(line, pos, tok)
        IF (LEN_TRIM(tok) == 0) EXIT
        READ(tok, '(I12)', IOSTAT=ios) v
        IF (ios /= 0) THEN
          is_alpha = .TRUE.; RETURN
        END IF
        nt = nt + 1
        IF (nt <= SIZE(counts)) counts(nt) = v
      END DO
    END SUBROUTINE parse_counts

    SUBROUTINE next_token(line, pos, tok)
      CHARACTER(*),INTENT(IN)    :: line
      INTEGER,INTENT(INOUT)      :: pos
      CHARACTER(*),INTENT(OUT)   :: tok
      INTEGER :: p0, ln
      tok = ''; ln = LEN(line)
      DO WHILE (pos <= ln .AND. line(pos:pos)==' ')
        pos = pos + 1
      END DO
      p0 = pos
      DO WHILE (pos <= ln .AND. line(pos:pos)/=' ')
        pos = pos + 1
      END DO
      IF (pos-1 >= p0) tok = line(p0:pos-1)
    END SUBROUTINE next_token

!=====================================================================
! mode selection: LOWEST | TRACK | GUESS (MODECAR)
!=====================================================================
    FUNCTION select_mode(n, eval, evec) RESULT(imode)
      INTEGER,INTENT(IN) :: n
      REAL(q),INTENT(IN) :: eval(n), evec(n,n)
      INTEGER            :: imode

      INTEGER      :: i, ibest, nim
      REAL(q)      :: ov, ovp
      REAL(q),ALLOCATABLE :: vg(:)

      ! penalized external modes (ef_external_penalty, eval ~ +K) are
      ! never candidates for the TS mode; LOWEST remains mode 1 (the
      ! eigenvalues are ascending, internal modes come first)
      nim = n
      DO i = 1, n
        IF (eval(i) > 1.0E3_q) THEN
          nim = i - 1
          EXIT
        END IF
      END DO
      IF (nim < 1) nim = 1
      imode = 1                                   ! LOWEST (eval ascending)
      IF (cmode == 'G' .AND. nstep <= 1) THEN
        ALLOCATE(vg(n)); vg = 0.0_dq
        IF (read_modecar('MODECAR', vg)) THEN
          ibest = 1; ovp = -1.0_q
          DO i = 1, nim
            ov = ABS(DOT_PRODUCT(vg, evec(:,i)))
            IF (ov > ovp) THEN; ovp = ov; ibest = i; END IF
          END DO
          imode = ibest
          IF (iu6_loc()>=0) WRITE(iu6_loc(),'(/A,I0,A,F6.3)') &
            ' EF: mode guess from MODECAR -> mode ',imode,' overlap ',ovp
        ELSE
          IF (iu6_loc()>=0) WRITE(iu6_loc(),'(/A)') ' EF: WARNING: no/bad MODECAR, using LOWEST'
        END IF
        DEALLOCATE(vg)
      ELSE IF (cmode == 'T' .AND. have_prev) THEN
        ibest = 1; ovp = -1.0_q
        DO i = 1, nim
          ov = ABS(DOT_PRODUCT(mode_prev, evec(:,i)))
          IF (eval(i) < 0.0_q) ov = ov + 0.5_q    ! prefer still-negative modes
          IF (ov > ovp) THEN; ovp = ov; ibest = i; END IF
        END DO
        imode = ibest
      END IF
    END FUNCTION select_mode

    FUNCTION read_modecar(fname, vg) RESULT(ok)
      CHARACTER(*),INTENT(IN) :: fname
      REAL(q),INTENT(OUT)     :: vg(3*nion)
      LOGICAL                 :: ok
      INTEGER :: unit, ios, i
      REAL(q) :: m(3,nion), nrm
      LOGICAL :: exists
      ok = .FALSE.; vg = 0.0_q
      INQUIRE(FILE=fname, EXIST=exists)
      IF (.NOT. exists) RETURN
      OPEN(NEWUNIT=unit, FILE=fname, STATUS='OLD', ACTION='READ', IOSTAT=ios)
      IF (ios /= 0) RETURN
      READ(unit,*,IOSTAT=ios) (m(1:3,i), i=1,nion)
      CLOSE(unit)
      IF (ios /= 0) RETURN
      vg = RESHAPE(m, (/3*nion/))
      nrm = SQRT(DOT_PRODUCT(vg,vg))
      IF (nrm < 1.0E-12_q) RETURN
      vg = vg/nrm
      ok = .TRUE.
    END FUNCTION read_modecar

    SUBROUTINE write_mode(fname, v)
      CHARACTER(*),INTENT(IN) :: fname
      REAL(q),INTENT(IN)      :: v(3*nion)
      REAL(q) :: m(3,nion), nrm
      INTEGER :: unit, i, ios
      IF (iu6_saved < 0) RETURN               ! master rank only
      m = RESHAPE(v, (/3,nion/))
      nrm = SQRT(SUM(m**2))
      IF (nrm > 1.0E-12_q) m = m/nrm
      OPEN(NEWUNIT=unit, FILE=fname, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
      IF (ios /= 0) RETURN
      WRITE(unit,'(3ES20.10)') (m(1:3,i), i=1,nion)
      CLOSE(unit)
    END SUBROUTINE write_mode

!=====================================================================
! EFDAT log (one line per step)
!=====================================================================
    SUBROUTINE write_efdat(toten, de, fmax, kappa, imode, nu1, nimag, slen)
      REAL(q),INTENT(IN) :: toten, de, fmax, kappa, nu1, slen
      INTEGER,INTENT(IN) :: imode, nimag
      INTEGER :: unit, ios
      LOGICAL :: exists
      IF (iu6_saved < 0) RETURN               ! master rank only
      INQUIRE(FILE='EFDAT', EXIST=exists)
      OPEN(NEWUNIT=unit, FILE='EFDAT', POSITION='APPEND', ACTION='WRITE', IOSTAT=ios)
      IF (ios /= 0) RETURN
      IF (.NOT. exists) WRITE(unit,'(A)') &
        '#  step         E(eV)        dE(eV)    Fmax(eV/A) kappa(eV/A^2)'// &
        '  mode   nu1(cm^-1)  nimag  |step|(A)'
      WRITE(unit,'(I6,F14.6,ES12.4,F11.5,ES11.3,I6,ES12.4,I5,F8.4)') &
        nstep, toten, de, fmax, kappa, imode, nu1, nimag, slen
      CLOSE(unit)
    END SUBROUTINE write_efdat

!=====================================================================
! EFHESSIAN snapshot (header + full symmetric matrix, row-major text)
!=====================================================================
    SUBROUTINE write_efhessian(n)
      INTEGER,INTENT(IN) :: n
      INTEGER :: unit, ios, i
      IF (iu6_saved < 0) RETURN               ! master rank only
      OPEN(NEWUNIT=unit, FILE='EFHESSIAN', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
      IF (ios /= 0) RETURN
      WRITE(unit,'(A)') '# EFHESSIAN v1 (VASP-EF)'
      WRITE(unit,'(A,I0,A,I0)') '# NIONS= ',nion,'  DOF= ',n
      WRITE(unit,'(A)') '# UNITS= eV/A^2 Cartesian, full symmetric, row-major'
      DO i = 1, n
        WRITE(unit,'(6ES24.16)') H(i,1:n)
      END DO
      CLOSE(unit)
    END SUBROUTINE write_efhessian

    SUBROUTINE read_efhessian(ni, posion, a, b)
      INTEGER,INTENT(IN)      :: ni
      REAL(q),INTENT(IN)      :: posion(3,ni), a(3,3), b(3,3)
      CHARACTER(512) :: line
      INTEGER :: n, unit, ios, i, k, cnt, nread
      LOGICAL :: exists
      n = 3*nion
      INQUIRE(FILE='EFHESSIAN', EXIST=exists)
      IF (.NOT. exists) THEN
        IF (iu6_loc()>=0) WRITE(iu6_loc(),'(/A)') &
          ' EF: WARNING: EF_READ_HESS set but EFHESSIAN missing; using model Hessian'
        CALL build_model_hessian(ni, posion, a, b)
        RETURN
      END IF
      OPEN(NEWUNIT=unit, FILE='EFHESSIAN', STATUS='OLD', ACTION='READ', IOSTAT=ios)
      IF (ios /= 0) THEN
        CALL build_model_hessian(ni, posion, a, b); RETURN
      END IF
      DO                                            ! skip header lines
        READ(unit,'(A)',IOSTAT=ios) line
        IF (ios /= 0) GO TO 900
        IF (line(1:1) /= '#') THEN
          BACKSPACE(unit)
          EXIT
        END IF
      END DO
      DO i = 1, n
        cnt = 0
        DO WHILE (cnt < n)
          nread = MIN(6, n-cnt)
          READ(unit,*,IOSTAT=ios) (H(i,cnt+k), k=1,nread)
          IF (ios /= 0) GO TO 900
          cnt = cnt + nread
        END DO
      END DO
      CLOSE(unit)
      RETURN
900   IF (iu6_loc()>=0) WRITE(iu6_loc(),'(/A)') &
        ' EF: FATAL: could not parse EFHESSIAN'
      STOP
    END SUBROUTINE read_efhessian

    SUBROUTINE build_model_hessian(ni, posion, a, b)
      INTEGER,INTENT(IN)      :: ni
      REAL(q),INTENT(IN)      :: posion(3,ni), a(3,3), b(3,3)
      INTEGER :: n, ia, i
      REAL(q) :: kd(3,nion), cart(3,nion)
      n = 3*nion
      IF (lbonds) THEN
        DO ia = 1, nion
          cart(:,ia) = MATMUL(posion(:,ia), a)
        END DO
        CALL ef_model_hessian_full(nion, znum, cart, a, b, bondk, btol, &
                                   hscale, H)
      ELSE
        CALL ef_model_hessian_diag(nion, znum, hscale, kd)
        H = 0.0_q
        DO ia = 1, nion
          DO i = 1, 3
            H(3*(ia-1)+i, 3*(ia-1)+i) = kd(i,ia)
          END DO
        END DO
      END IF
    END SUBROUTINE build_model_hessian

!=====================================================================
! element tables: symbol -> atomic number / standard atomic mass
!=====================================================================
    FUNCTION symbol_z(sym) RESULT(zz)
      CHARACTER(2),INTENT(IN) :: sym
      INTEGER :: zz
      zz = symbol_lookup(sym, 1)
    END FUNCTION symbol_z

    FUNCTION symbol_mass(sym) RESULT(mm)
      CHARACTER(2),INTENT(IN) :: sym
      REAL(q) :: mm
      mm = symbol_lookup(sym, 2) * 0.001_q     ! table stores mass x 1000
    END FUNCTION symbol_mass

    ! what=1: atomic number; what=2: atomic mass in amu x 1000
    FUNCTION symbol_lookup(sym, what) RESULT(val)
      CHARACTER(2),INTENT(IN) :: sym
      INTEGER,INTENT(IN)      :: what
      INTEGER                 :: val
      CHARACTER(2),PARAMETER  :: syl(1:60) = (/ &
        'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
        'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca', &
        'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
        'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr', &
        'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
        'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd' /)
      INTEGER,PARAMETER :: zv(1:60) = (/ &
         1,  2,  3,  4,  5,  6,  7,  8,  9, 10, &
        11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
        21, 22, 23, 24, 25, 26, 27, 28, 29, 30, &
        31, 32, 33, 34, 35, 36, 37, 38, 39, 40, &
        41, 42, 43, 44, 45, 46, 47, 48, 49, 50, &
        51, 52, 53, 54, 55, 56, 57, 58, 59, 60 /)
      INTEGER,PARAMETER :: mx(1:60) = (/ &      ! (amu*1000) as integer
         1008, 4003, 6940, 9012,10810,12011,14007,15999,18998,20180, &
        22990,24305,26982,28085,30974,32060,35450,39948,39098,40078, &
        44956,47867,50942,51996,54938,55845,58933,58693,63546,65380, &
        69723,72630,74922,78971,79904,83798,85468,87620,88906,91224, &
        92906,95950,98000,101070,102910,106420,107870,112410,114820, &
        118710,121760,127600,126900,131290,132910,137330,138910,140120,140910, &
        144240 /)
      INTEGER :: i
      val = 0
      DO i = 1, 60
        IF (TRIM(sym) == TRIM(syl(i))) THEN
          IF (what == 1) THEN
            val = zv(i)
          ELSE
            val = mx(i)
          END IF
          RETURN
        END IF
      END DO
      IF (what == 1) THEN
        val = 6
      ELSE
        val = 12011
      END IF
    END FUNCTION symbol_lookup

    FUNCTION upcase2(s) RESULT(u)
      CHARACTER(*),INTENT(IN) :: s
      CHARACTER(LEN(s))       :: u
      INTEGER :: i, ic
      u = s
      DO i = 1, LEN(s)
        ic = IACHAR(s(i:i))
        IF (ic >= IACHAR('a') .AND. ic <= IACHAR('z')) u(i:i) = ACHAR(ic-32)
      END DO
    END FUNCTION upcase2

  END MODULE ef_ts
