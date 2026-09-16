# VASP-EF transition-state search: HCN -> HNC hydrogen shift
# Workflow: relax HCN/HNC minima (IBRION=2), displace H to the bent
# TS neighbourhood (POSCAR.bent), then polish with VASP-EF.
SYSTEM = HCN isomerization TS search (VASP-EF)
ENCUT  = 400
EDIFF  = 1E-7
ISMEAR = 0
SIGMA  = 0.05
ISPIN  = 1
LWAVE  = .FALSE.
LCHARG = .FALSE.
ISYM   = 0        # recommended near saddle points

# --- VASP-EF activation protocol (same convention as VTST) ---
IBRION = 3
POTIM  = 0
EF_TS  = .TRUE.

# --- optimizer settings: polish from the bent seed ---
EDIFFG      = -0.05
NSW         = 300
EF_MAXSTEP  = 0.08
EF_TRUST    = ADAPTIVE
EF_MODE     = GUESS
EF_GUESS_K  = -0.5
EF_UPDATE   = BOFILL
EF_BONDS    = .TRUE.
EF_CALCFC   = .TRUE.     # exact Hessian at start (small system: all atoms)
EF_FCDELTA  = 0.05
EF_FDREFINE = .TRUE.     # keep the TS-mode curvature calibrated
EF_FDREF_K  = 5
EF_FDDELTA  = 0.04
EF_PRINT    = 2
EF_ROTFIX   = .TRUE.
