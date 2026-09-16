# VASP-EF Developer Guide

## Architecture

```
src/ef_ts.F          MODULE ef_core   -- pure math (no VASP deps)
                       MODULE ef_ts    -- VASP glue + optimizer
patches/             chain.F hook variants (stock / vtst)
scripts/             Python workflow tools
tests/               standalone unit tests
```

One source file contains two modules. `ef_core` is compiled standalone
for unit tests; `ef_ts` is compiled inside VASP.

## Hook

`chain.F` `chain_force` calls `ef_ts_step(nions,posion,toten,force,a,b,iu6)`.
Returns `.TRUE.` → VASP skips its own update (IBRION=3, POTIM=0).

## State machines

### CALCFC (exact FD Hessian)

`fcstage` 0→1→2→(0) per column. Two displaced single points per column.
36 NSW iterations for 6 atoms. Assembles exact H by central differences.

### FD curvature probe

`fdpending` 0→1→2→3→(0). Two displaced single points along the TS mode
give its curvature by central force difference. Locked via ef_seed_mode.
3 NSW iterations per probe.

### Accept/reject

Each proposal checked on next force call. Rejection restores pos_good,
halves trust, does NOT apply Hessian update (deferred to acceptance).

## Phase logic

MIN phase: damped quasi-Newton (RFO shift below b_min).
CLIMB phase: P-RFO with followed mode uphill.
One-way latch: MIN → CLIMB when Fmax < EF_FSWITCH or κ < 0.

## Mode selection

LOWEST: first internal mode (most negative eigenvalue).
TRACK: max overlap to previous mode, +0.5 bias for negative modes.
GUESS: max overlap with MODECAR at step 1; seeds Hessian via
ef_seed_mode, then switches to TRACK behavior.

## Model Hessian

`ef_model_hessian_full(inhess, nat, zat, cart, a, b, H)`:
diagonal element terms (Fischer-Almlof table) plus springs along
interatomic directions. Parameterized by `inhess`:
1=LINDH springs, 2=SCHLEGEL springs, 3=ALMLOF diagonal only.

## Key algorithms

### Single-shift EF

One shift λ ∈ (b_m, b_next) solves Σg²/(b−λ)² = min(1,R)².
Ternary search for interior minimum + bisection for the upper root.
Falls back to minimum-length shift when unreachable.

### Projection deflation (mode seeding)

H' = k·v·vᵀ + (I − v·vᵀ)·H·(I − v·vᵀ). Makes v an exact eigenvector.

### Trust-region rejection

Proposals raising E by more than EF_EMAX are rejected. Trust halves.
At the trust floor, H is reset and phase returns to MIN.

## Species parsing

POTCAR TITEL lines give element symbols. POSCAR line 6/7 gives counts.
Both parsed by ef_ts.F directly (no VASP RDATAB).

## File formats

### EFHESSIAN

Header (# NIONS=, # DOF=, # UNITS=) then full symmetric matrix,
row-major, 6 values per line.

### EFDAT

```
step  E  dE  Fmax  kappa  mode  nu1  nimag  |step|
```

### MODECAR / NEWMODECAR

3N Cartesian components, three per line (3ES20.10), unit-normalized.

## MPI model

All ranks redundantly do the optimizer arithmetic (as in VTST).
No additional communication. Bit-reproducible across rank counts.

## Adding a feature

1. New INCAR tag: read in ef_init, declare SAVE variable.
2. New mode selection: add branch in select_mode, guard by cmode.
3. New Hessian update: add subroutine in ef_core following Bofill/Powell.
4. New output file: add writer with iu6>=0 master-rank guard.
5. Add unit test in tests/test_ef_core.F90 or test_scripts.py.

## Testing

```bash
bash tests/run_tests.sh        # Fortran core (gfortran + LAPACK)
python tests/test_scripts.py   # Python scripts (stdlib only)
```

## Code style

* All-English comments
* Free-form Fortran, .F extension (cpp preprocessing)
* REAL(q) from VASP prec in ef_ts.F; REAL(dq) in ef_core.F
* Prefix ef_ on public procedures
* Master-rank guard IF (iu6>=0) on stdout and file writes
* One subroutine = one purpose

## Known limitations

* Fixed cell (ISIF ≤ 2)
* Selective Dynamics not handled
* Cartesian coordinates only (no internal coordinate optimization)
* All ranks hold redundant copy of the Hessian (O((3N)²) memory)
