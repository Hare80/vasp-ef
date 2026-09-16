# H/Si(001) diffusion — structural demonstration only

This example shows the file layout for a surface TS search, but the
slab structure has NOT been validated against literature. Use as a
template for setting up your own system, not as a quantitative
benchmark.

## What is here

| file | content |
|---|---|
| `POSCAR` | Si(001) slab (16 Si) + H at a bridge site |
| `MODECAR` | unit mode: H moving vertically |
| `INCAR` | VASP-EF tags for TS search |
| `make_slab.py` | regenerates the slab (requires ASE) |
| `stage2_ts.sh` | two-step workflow: relax → TS search |

## To make this a proper benchmark

1. Obtain a properly relaxed Si(001)-(2×1) reconstructed surface
2. Add H at a known adsorption site (monohydride)
3. Relax the full slab + H (IBRION=2, no EF)
4. Then run the TS search for H hopping to an adjacent site
5. Compare barrier with literature (~1.5–1.9 eV)

Without these steps the forces and energies are not meaningful.
