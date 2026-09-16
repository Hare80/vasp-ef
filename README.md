# VASP-EF

**Eigenvector-following transition-state optimizer for VASP.**
Partitioned-RFO step, Bofill update, mode tracking, curvature probing.

Reference: VTST tools (Henkelman group), Gaussian `Opt=TS`, ORCA `!OptTS`.

Author: Yihong Lian — [github.com/Hare80](https://github.com/Hare80)

## Install

Copy into your VASP `src/` directory (VASP 6.6.0 tested):

```
cp src/ef_ts.F  <vasp>/src/
cp patches/chain.F <vasp>/src/chain.F
```

Add one line before `chain.o` in `<vasp>/src/.objects`:

```
  ef_ts.o \
```

Rebuild: `make std` (or `make gam` / `make ncl`).

For CMake builds the `.objects` change is picked up automatically.
See `INSTALL.md` for detailed instructions and troubleshooting.

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

## Tags

| tag | default | description |
|---|---|---|
| `EF_TS` | .FALSE. | activate |
| `EF_INHESS` | 1 | 1=LINDH 2=SCHLEGEL 3=ALMLOF 4=READ |
| `EF_CALCFC` | .FALSE. | exact Hessian at start |
| `EF_FCHESS` | 0 | rebuild exact Hessian every N steps |
| `EF_ACTIVE` | all | active atoms for FD Hessian |
| `EF_MODE` | 1 | 1=LOWEST 2=TRACK 3=GUESS |
| `EF_UPDATE` | 1 | 1=Bofill 2=Powell |
| `EF_MAXSTEP` | 0.2 | trust radius (Å) |
| `EF_FDREF_K` | 10 | curvature probe interval (0=off) |
| `EF_NOCHECK` | .FALSE. | skip imaginary-mode check |
| `EF_ROTFIX` | .FALSE. | project rotations |
| `EDIFFG` | — | VASP force criterion (unchanged) |

See `docs/MANUAL.md` for the full user manual with workflow guidance.

## Examples

| directory | system | notes |
|---|---|---|
| `examples/hcn-isomerization/` | HCN ↔ HNC | 3 atoms, barrier 2.0 eV |
| `examples/h3-exchange/` | H + H₂ → H + H₂ | 3 atoms, converged in 7 steps |
| `examples/nh3-inversion/` | NH₃ umbrella | 4 atoms, barrier 0.21 eV |
| `examples/sn2-cl-exchange/` | SN2 Cl exchange | 6 atoms, shallow double-well |
| `examples/h-si001/` | H/Si(001) diffusion | 17 atoms, surface system |
| `examples/bench-set/` | batch builder | constructs all inputs |

## Testing

```
bash tests/run_tests.sh          # Fortran math core
python tests/test_scripts.py     # workflow scripts
```

## Citation

Y. Lian, VASP-EF: eigenvector-following transition-state optimizer for
VASP, https://github.com/Hare80/vasp-ef (2026).

Please also cite VTST if you use NEB/dimer workflows alongside EF.

## License

MIT
