# VASP-EF User Manual

## Protocol

INCAR must contain (same convention as VTST):

```
IBRION = 3
POTIM  = 0
EF_TS  = .TRUE.
```

VASP's own optimizers are disabled; VASP-EF updates positions.

## Workflow

### 1. Prepare a starting structure

Start near the saddle point. Good sources: relaxed scan endpoint, NEB
image near the barrier, displaced minimum along the soft mode, or a
previous TS structure (CONTCAR).

### 2. Choose the mode to follow

| EF_MODE | when to use |
|---|---|
| LOWEST | follow the most negative eigenvalue; simple barriers |
| TRACK | follow the mode most similar to the previous step; robust against eigenvalue crossings |
| GUESS | lock to the MODECAR direction at step 1; best when you know the reaction coordinate |

For `EF_MODE=GUESS`, place a `MODECAR` file in the run directory
(3N Cartesian components, one line per atom, same format as VTST).
Generate with `scripts/modemake.py POS1 POS2 MODECAR`.

### 3. Run

```
mpirun -n N vasp_std
```

### 4. Monitor

```
python scripts/efstat.py
```

### 5. Verify

```
python scripts/hessfreq.py
```

## Tag reference

### Activation

| tag | default | description |
|---|---|---|
| `EF_TS` | .FALSE. | activate optimizer |
| `EF_MAXSTEP` | 0.2 | trust radius (Å) |
| `EF_TRUST` | ADAPTIVE | or FIXED |
| `EF_PRINT` | 1 | output verbosity (1=standard) |

### Mode following

| tag | default | description |
|---|---|---|
| `EF_MODE` | LOWEST | LOWEST, TRACK, or GUESS |
| `EF_GUESS_K` | −0.5 | curvature of the seeded MODECAR direction |

### Hessian

| tag | default | description |
|---|---|---|
| `EF_INHESS` | 1 | 1=LINDH, 2=SCHLEGEL, 3=ALMLOF, 4=READ |
| `EF_HSCALE` | 1.0 | scale model Hessian |
| `EF_BONDS` | .TRUE. | bond-aware springs |
| `EF_BOND_K` | 80 | bond spring constant (eV/Å²) |
| `EF_BOND_TOL` | 0.45 | bond detection tolerance (Å) |
| `EF_CALCFC` | .FALSE. | exact FD Hessian at start |
| `EF_FCDELTA` | 0.05 | FD displacement (Å) |
| `EF_FCHESS_N` | 0 | rebuild every N steps |
| `EF_ACTIVE` | all | active atoms for FD Hessian |
| `EF_READ_HESS` | .FALSE. | read from EFHESSIAN |
| `EF_HESS_NWRITE` | 0 | write every N steps |

### Update

| tag | default | description |
|---|---|---|
| `EF_UPDATE` | BOFILL | or POWELL |

### Step control

| tag | default | description |
|---|---|---|
| `EF_FSWITCH` | 0.5 | MIN→CLIMB force gate (eV/Å) |
| `EF_EMAX` | 0.2 | proposal rejection threshold (eV) |
| `EF_NOCHECK` | .FALSE. | suppress imaginary-mode warning |
| `EF_ROTFIX` | .FALSE. | project rotations |

### Misc

| tag | default | description |
|---|---|---|
| `EF_FDREF_K` | 10 | curvature probe interval |
| `EF_PRINT` | 1 | output verbosity |

## Output files

| file | content |
|---|---|
| `EFDAT` | step, E, dE, Fmax, κ, mode, ν₁, nimag, \|step\| |
| `NEWMODECAR` | followed TS mode (unit vector, Cartesian) |
| `EFHESSIAN` | Hessian snapshot |
| `MODECAR` | read at startup with EF_MODE=GUESS |

## Interpreting EFDAT

* **kappa**: curvature along the followed mode. Negative near a TS.
* **nu1**: lowest mass-weighted frequency. Negative = imaginary.
* **nimag**: number of imaginary modes. Should converge to 1.
* **|step|**: actual displacement. Decreases near convergence.

## Troubleshooting

| problem | fix |
|---|---|
| converges to a minimum (nimag=0) | displace along the reaction coordinate first, or use EF_MODE=GUESS with a MODECAR |
| bounces without settling | reduce EF_MAXSTEP; ensure EDIFF is tight enough (1e-7) |
| imaginary modes appear then vanish | mode crossing; try EF_MODE=TRACK |
| forces don't converge | check SCF convergence (EDIFF), increase NSW |
| stuck bouncing on a wall | reduce EF_MAXSTEP to 0.05; check for atom overlaps |

## Limitations

* Fixed unit cell (ISIF ≤ 2)
* Selective Dynamics not explicitly handled
* All MPI ranks do redundant optimizer arithmetic (same as VTST)
* Cartesian coordinates only
