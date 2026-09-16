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

Start near the saddle point. Good sources:

- relaxed scan endpoint
- NEB image near the barrier
- displaced minimum along the soft mode
- previous TS structure (CONTCAR)

If starting from a minimum, first displace along the expected reaction
coordinate by ~0.1–0.3 Å. A minimum itself has no negative mode to
climb, and the optimizer may wander.

### 2. Choose the mode to follow

| EF_MODE | when to use |
|---|---|
| 1 LOWEST | follow the most negative eigenvalue; simple barriers |
| 2 TRACK | follow the mode most similar to the previous step; robust against eigenvalue crossings |
| 3 GUESS | lock to the MODECAR direction at step 1; best when you know the reaction coordinate |

For `EF_MODE=3`, place a `MODECAR` file in the run directory
(3N Cartesian components, one line per atom, same format as VTST).
Generate with `scripts/modemake.py POS1 POS2 MODECAR`.

### 3. Run

```
mpirun -n N vasp_std
```

Monitor convergence:

```
python scripts/efstat.py
```

The optimization stops when Fmax < |EDIFFG| (VASP criterion). VASP-EF
then checks the mass-weighted spectrum: one imaginary mode = TS found.

### 4. Verify

```
python scripts/hessfreq.py
```

Confirms the converged structure has exactly one imaginary frequency
and prints all mass-weighted frequencies from EFHESSIAN.

### 5. Restart

```
python scripts/efrestart.py old_run new_run
```

Copies CONTCAR → POSCAR, NEWMODECAR → MODECAR, EFHESSIAN, and sets
EF_READ_HESS = .TRUE. in INCAR.

## Tag reference

### Activation

| tag | default | description |
|---|---|---|
| `EF_TS` | .FALSE. | activate optimizer |
| `EF_MAXSTEP` | 0.2 | trust radius (Å); reduce to 0.05–0.10 for hard cases |
| `EF_PRINT` | — | removed; output is always printed |

### Mode following

| tag | default | description |
|---|---|---|
| `EF_MODE` | 1 | 1=lowest eigenvalue, 2=track previous mode, 3=align to MODECAR |
| `EF_GUESS_K` | −0.5 | curvature of the seeded MODECAR direction; negative = climb immediately |

### Hessian

| tag | default | description |
|---|---|---|
| `EF_INHESS` | 1 | initial Hessian: 1=LINDH, 2=SCHLEGEL, 3=ALMLOF, 4=READ |
| `EF_CALCFC` | .FALSE. | compute exact FD Hessian at start (Gaussian CalcFC) |
| `EF_FCHESS` | 0 | rebuild exact Hessian every N steps (Gaussian RecalcFC) |
| `EF_ACTIVE` | all | active atoms for FD Hessian (e.g. `1 2 5` for adsorbate + nearby substrate) |

The exact Hessian (CALCFC or FCHESS) always supersedes the model
Hessian. For large systems use EF_ACTIVE to restrict the FD to the
chemically relevant atoms.

### Hessian update

| tag | default | description |
|---|---|---|
| `EF_UPDATE` | 1 (Bofill) | SR1/PSB interpolation; keeps negative curvature |
| — | — | 2 (Powell) flips positive curvature along the step |

### Step control

| tag | default | description |
|---|---|---|
| `EF_MAXSTEP` | 0.2 | maximum step length per iteration (Å) |
| `EF_TRUST` | ADAPTIVE | radius grows on accepted steps, shrinks on rejections |
| `EF_FDREF_K` | 10 | curvature probe interval; 0 = off |

### Convergence reporting

| tag | default | description |
|---|---|---|
| `EF_NOCHECK` | .FALSE. | if .TRUE., suppress imaginary-mode warning |
| `EDIFFG` | (VASP) | force criterion; VASP stops on Fmax < \|EDIFFG\| |
| — | — | imaginary-mode count is advisory (doesn't stop VASP) |

### Misc

| tag | default | description |
|---|---|---|
| `EF_ROTFIX` | .FALSE. | project global rotations (gas-phase molecules only) |
| `EF_FSWITCH` | 0.5 | internal: force threshold for phase transition |

## Output files

| file | content |
|---|---|
| `EFDAT` | one line per step: step, E, dE, Fmax, κ, mode, ν₁, nimag, \|step\| |
| `NEWMODECAR` | followed TS mode (unit vector, Cartesian, VTST format) |
| `EFHESSIAN` | Cartesian Hessian; written at convergence and on request |
| `MODECAR` | read at startup with EF_MODE=3 |

## Interpreting EFDAT

```
#  step      E(eV)     dE(eV)   Fmax   kappa    mode  nu1(cm-1) nimag |step|
   15   -17.7053   -0.0024   0.3006  -2.389     1   -410.4      1  0.0196
```

* **kappa**: curvature along the followed mode. Negative near a TS.
* **nu1**: lowest mass-weighted frequency. Becomes imaginary (negative)
  as the structure approaches the saddle.
* **nimag**: number of imaginary modes. Should converge to 1 for a TS.
* **|step|**: actual displacement. Decreases as the trust radius shrinks
  near convergence.

## Troubleshooting

| problem | fix |
|---|---|
| converges to a minimum (nimag=0) | displace along the reaction coordinate first, or use EF_MODE=3 with a MODECAR |
| bounces without settling | reduce EF_MAXSTEP; ensure EDIFF is tight enough (1e-7) |
| imaginary modes appear then vanish | mode crossing; try EF_MODE=TRACK |
| forces don't converge | check SCF convergence (EDIFF), increase NSW |
| stuck bouncing on a wall | reduce EF_MAXSTEP to 0.05; check for atom overlaps in CONTCAR |

## Limitations

* Fixed unit cell (ISIF ≤ 2)
* Selective Dynamics not explicitly handled
* All MPI ranks do redundant optimizer arithmetic (same as VTST)
* Cartesian coordinates only (no internal coordinate optimization)
