# VASP-EF Developer Guide

## Architecture

```
src/ef_ts.F
  MODULE ef_core    -- pure math: eigensolver, P-RFO, Bofill/Powell,
                       model Hessian builders, INCAR readers,
                       mass-weighted frequencies
  MODULE ef_ts      -- VASP glue: state machine, species parsing,
                       file I/O, optimizer entry point

patches/            chain.F hook (stock and VTST variants)
scripts/            Python workflow tools
tests/              standalone unit tests
```

One source file (`ef_ts.F`) contains both the optimizer and its
numerical core. It is hooked into VASP via a single call at the top
of `chain_force` in `chain.F`:

```fortran
IF (ef_ts_step(nions,posion,toten,force,a,b,iu6)) RETURN
```

When `EF_TS=.TRUE.` in INCAR, this returns `.TRUE.` and VASP's own
position update is skipped (via `IBRION=3, POTIM=0`).

## State machines

The optimizer uses three interlocking state machines driven by the
`optflag` protocol (VASP re-evaluates forces when positions change):

### CALCFC (exact FD Hessian)

Consumes `6*Nact` NSW iterations. For each active atom and each
Cartesian direction, two displaced single points give one column of
the exact Hessian by central differences. The assembled H is exact
at the current geometry (Gaussian CalcFC analog).

### FD curvature probe

Two displaced single points along the followed TS mode give its
curvature by central force difference. Locked into H via
`ef_seed_mode` (projection deflation). Costs 3 NSW iterations.
Triggered every `EF_FDREF_K` steps during the climb phase.

### Accept/reject

Each proposal is checked on the next force evaluation. If the energy
rose by more than `EF_EMAX`, the step is rejected and retried with a
halved trust radius. Hessian updates are deferred until acceptance,
so rejected steps never poison H.

## Phase logic

1. **MIN phase**: damped quasi-Newton (RFO shift below b_min). Reduces
   Fmax until < EF_FSWITCH (0.5 eV/Å) or negative curvature appears.
2. **CLIMB phase**: P-RFO with the followed mode uphill. Trust radius
   ratchets down near convergence (3×/10× EDIFFG thresholds).

Transitions are one-way: MIN → CLIMB. A Hessian reset returns to MIN.

## Mode selection

`select_mode` returns the eigenmode index to follow:

* LOWEST (default): first internal mode (most negative eigenvalue)
* TRACK: mode with maximum overlap to the previous step's mode,
  with a +0.5 bias toward still-negative modes
* GUESS: step 1 only, aligned to the MODECAR direction via maximum
  overlap; seeds the Hessian along the guess with curvature GUESS_K

Penalized external modes (eval > 10³) are excluded from selection.

## Key algorithms

### Single-shift eigenvector following

One shift λ ∈ (b_m, b_next) solves:

```
Σᵢ gᵢ²/(bᵢ−λ)² = R²
```

by ternary search for the U-shaped minimum, then bisection on the
upper branch. λ near bₘ → climb-dominant; λ near b_next → relax-
dominant. The balance is automatic.

### Projection deflation (mode seeding)

```fortran
H' = k·v·vᵀ + (I − v·vᵀ)·H·(I − v·vᵀ)
```

Makes v an exact eigenvector with eigenvalue k. Used to inject a
user-supplied MODECAR direction into the model Hessian.

### Bofill update

```fortran
φ = (sᵀy)² / ((sᵀs)(yᵀy))
H ← H + (ηsᵀ+sηᵀ)/(sᵀs) − (sᵀη)ssᵀ/(sᵀs)²
      + φ[ηηᵀ/(sᵀy) − (sᵀη)ssᵀ/(sᵀs)²]
```

### Trust-region rejection

Proposals that raise E by more than `EF_EMAX` are rejected. Hessian
updates are deferred until acceptance (rejected steps never poison H).
Trust radius: grow ×1.5 on accepted sub-cap steps, halve on rejections.
Endgame ratchet caps the radius when Fmax approaches the criterion.

## File formats

### EFHESSIAN

```
# EFHESSIAN v1 (VASP-EF)
# NIONS= <n>  DOF= <3n>
# UNITS= eV/A^2 Cartesian, full symmetric, row-major
<3n rows × 6 values per row>
```

### EFDAT

One line per step, space-separated:

```
step  E  dE  Fmax  kappa  mode  nu1  nimag  |step|
```

### MODECAR / NEWMODECAR

3N Cartesian components, three per line, format `3ES20.10`,
unit-normalized. Same convention as VTST dimer.

## Adding a new feature

1. **New INCAR tag**: add a read call in `ef_init`, declare the SAVE
   variable at module scope, document in the file header and
   `docs/MANUAL.md`.
2. **New mode selection strategy**: add a branch in `select_mode`,
   guarded by a new `cmode` character.
3. **New Hessian update**: add a subroutine in `ef_core.F` following
   the Bofill/Powell pattern (accept/reject with `applied` flag),
   add a case in `ef_ts_step` and an `EF_UPDATE` value.
4. **New output file**: add a writer subroutine with the
   `IF (iu6<0) RETURN` master-rank guard and call it from the
   diagnostics section.

## Testing

```bash
# Fortran core (standalone, no VASP needed)
bash tests/run_tests.sh

# Python workflow scripts
python tests/test_scripts.py
```

Tests use a stub `prec` module so `ef_ts.F` compiles outside VASP.
The math core (`ef_core.F`) compiles standalone; only `ef_ts.F`
needs VASP's `prec` for the `q` kind parameter.

## Code style

* All-English comments and docs
* Free-form Fortran, `.F` extension (cpp preprocessing)
* `REAL(q)` from VASP `prec` in `ef_ts.F`; `REAL(dq)` in `ef_core.F`
  (both resolve to double precision)
* Prefix `ef_` on all public procedures and variables
* Master-rank guard `IF (iu6>=0)` on all file writes and stdout
* One subroutine = one clear purpose; prefer many small routines
* Tags use `EF_` prefix; single-word or compound words in CAPS

## Dependencies

* VASP ≥ 6.x source tree
* BLAS/LAPACK (already linked by VASP)
* Python 3 for workflow scripts (stdlib only)
