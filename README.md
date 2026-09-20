# VASP-EF

Eigenvector-following transition-state optimizer for VASP.
Partitioned-RFO step, Bofill update, mode tracking, curvature probing.

References: VTST tools (Henkelman group), Gaussian `Opt=TS`, ORCA `!OptTS`.

Author: Yihong Lian

## Install

Run the install script from the repo root, pointing to your VASP source:

```
bash install.sh --src /path/to/vasp.x.x.x
```

This copies `src/ef_ts.F` into your VASP `src/` directory, patches
`src/chain.F` with a hook, and adds `ef_ts.o` to `src/.objects`.
Rebuild VASP as usual (`make std` or `make gam` / `make ncl`).

Alternatively, apply `patches/chain_F_hook.stock.patch` (or `.vtst.patch`
if your chain.F already has VTST) and copy `src/ef_ts.F` manually.

## Quick start

INCAR:

```
IBRION = 3
POTIM  = 0
EF_TS  = .TRUE.
EDIFFG = -0.05
```

Run VASP as usual. Output: `EFDAT` (per-step log), `NEWMODECAR` (followed
mode), `EFHESSIAN` (Hessian snapshot), stdout `EF:` lines.

Monitor: `python scripts/efstat.py`

Verify frequencies: `python scripts/hessfreq.py`

## All tags

| tag | default | description |
|---|---|---|
| `EF_TS` | .FALSE. | activate |
| `EF_MAXSTEP` | 0.2 | trust radius (Å) |
| `EF_TRUST` | ADAPTIVE | ADAPTIVE or FIXED |
| `EF_MODE` | LOWEST | LOWEST / TRACK / GUESS |
| `EF_UPDATE` | BOFILL | BOFILL or POWELL |
| `EF_NOCHECK` | .FALSE. | skip imaginary-mode check |
| `EF_HSCALE` | 1.0 | scale model Hessian |
| `EF_BONDS` | .TRUE. | bond-aware model Hessian |
| `EF_BOND_K` | 80 | bond spring constant (eV/Å²) |
| `EF_BOND_TOL` | 0.45 | bond detection tolerance (Å) |
| `EF_FSWITCH` | 0.5 | MIN→CLIMB force gate (eV/Å) |
| `EF_EMAX` | 0.2 | proposal rejection threshold (eV) |
| `EF_GUESS_K` | −0.5 | seeded MODECAR curvature (eV/Å²) |
| `EF_CALCFC` | .FALSE. | exact Hessian at start |
| `EF_FCDELTA` | 0.05 | FD Hessian displacement (Å) |
| `EF_FDREFINE` | .FALSE. | TS-mode curvature probe |
| `EF_FDREF_K` | 10 | probe interval (steps) |
| `EF_FDDELTA` | 0.04 | probe displacement (Å) |
| `EF_FCHESS_N` | 0 | rebuild exact Hessian every N steps |
| `EF_READ_HESS` | .FALSE. | read initial Hessian from EFHESSIAN |
| `EF_HESS_NWRITE` | 0 | write EFHESSIAN every N steps |
| `EF_ACTIVE` | all | active atoms for FD Hessian |
| `EF_DIAG` | FULL | P-RFO needs the full spectrum |
| `EF_ROTFIX` | .FALSE. | project rotations (gas-phase molecules) |
| `EDIFFG` | — | VASP force criterion (unchanged) |

See `docs/MANUAL.md` for the full user manual and workflow guidance.

## Examples

| directory | system | notes |
|---|---|---|
| `examples/hcn-isomerization/` | HCN ↔ HNC | 3 atoms, barrier 2.0 eV |
| `examples/h3-exchange/` | H + H₂ → H + H₂ | 3 atoms, 7-step convergence |
| `examples/nh3-inversion/` | NH₃ umbrella | 4 atoms, barrier 0.21 eV |
| `examples/sn2-cl-exchange/` | SN2 Cl exchange | 6 atoms, shallow double-well |
| `examples/bench-set/` | batch builder | all five setups |
| `examples/h-si001/` | H/Si(001) demo | structural demo only, not validated |

## Scripts

| script | purpose |
|---|---|
| `efstat.py` | print convergence status from EFDAT |
| `modemake.py` | build MODECAR from two structures |
| `hessfreq.py` | mass-weighted frequencies from EFHESSIAN |
| `efrestart.py` | restart from CONTCAR + NEWMODECAR + EFHESSIAN |

## Testing

```
bash tests/run_tests.sh
python tests/test_scripts.py
```

## Citation

Y.H. Lian, VASP-EF: eigenvector-following transition-state optimizer for
VASP (2026). https://github.com/Hare80/vasp-ef

## License

MIT
