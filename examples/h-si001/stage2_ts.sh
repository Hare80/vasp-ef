#!/usr/bin/env bash
# stage2: from the relaxed minimum, place H at the dimer bridge and run
# the VASP-EF transition-state search.  Run after INCAR.min finishes.
set -euo pipefail
cd "$(dirname "$0")"
python - <<'PY'
from ase.io import read, write
import numpy as np
at = read('CONTCAR', format='vasp')
tops = [a.index for a in at if a.symbol=='Si' and a.position[2] > 12.0]
h = at.positions[0]
a_idx = min(tops, key=lambda i: np.linalg.norm(at.positions[i]-h))
b_idx = min((i for i in tops if i != a_idx),
            key=lambda i: np.linalg.norm(at.positions[i]-at.positions[a_idx]))
pair = [a_idx, b_idx]
print('dimer Si:', pair)
mid = (at.positions[pair[0]] + at.positions[pair[1]])/2
at.positions[0] = mid + np.array([0.0, 0.0, 1.10])
write('POSCAR', at, format='vasp', direct=True, vasp5=True, sort=True)
print('H -> bridge', at.positions[0])
PY
cp MODECAR MODECAR.keep 2>/dev/null || true
cat > MODECAR <<'MODE'
  0.0000000000E+00  0.0000000000E+00  1.0000000000E+00
MODE
cat > INCAR <<'INCEOF'
SYSTEM = H diffusion on Si(001), bridging TS (VASP-EF)
ENCUT  = 250
EDIFF  = 1E-6
ISMEAR = 0
SIGMA  = 0.05
ISPIN  = 1
LWAVE  = .FALSE.
LCHARG = .FALSE.
ISYM   = 0
IBRION = 3
POTIM  = 0
EF_TS  = .TRUE.
EDIFFG = -0.05
NSW    = 400
EF_MAXSTEP  = 0.15
EF_TRUST    = ADAPTIVE
EF_MODE     = GUESS
EF_GUESS_K  = -0.5
EF_UPDATE   = BOFILL
EF_BONDS    = .TRUE.
EF_CALCFC   = .TRUE.
EF_ACTIVE   = 1 3 4
EF_FCDELTA  = 0.05
EF_FDREFINE = .TRUE.
EF_FDREF_K  = 10
EF_FDDELTA  = 0.05
EF_PRINT    = 2
INCEOF
rm -f OSZICAR OUTCAR CONTCAR XDATCAR EFDAT NEWMODECAR EFHESSIAN REPORT vasprun.xml vaspout.h5 CHG WAVECAR PCDAT DOSCAR EIGENVAL IBZKPT
"/d/vasp-win-build/repo/vasp-windows-msys2-build/build_work/stock/20260915-050351/vasp-6.6.0-msys2-portable/run.bat"
