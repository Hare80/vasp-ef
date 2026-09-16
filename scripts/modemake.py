#!/usr/bin/env python3
"""modemake.py -- build a MODECAR initial-mode guess from two structures.

This is the VASP-EF counterpart of vtstscripts' modemake.pl.  The MODECAR
it writes is read by VASP-EF when INCAR has EF_MODE=GUESS: the optimizer
follows the Hessian eigenmode with the largest overlap with this vector.

Usage:
    python modemake.py POSCAR_REACTANT POSCAR_PRODUCT [MODECAR]

The mode is the minimum-image displacement (product - reactant) in
Cartesian coordinates, normalized to unit length -- the same convention
as NEWMODECAR written by VASP-EF.  Both structures must share the same
lattice (the lattice of the first file is used).

Example:
    python modemake.py 00/CONTCAR 04/CONTCAR MODECAR
"""
from __future__ import annotations

import sys

from vefcommon import (direct_to_cart, min_image_delta, read_poscar,
                       write_modecar)


def main(argv):
    if len(argv) < 3 or len(argv) > 4:
        print(__doc__)
        return 2
    pos1, pos2 = argv[1], argv[2]
    out = argv[3] if len(argv) > 3 else "MODECAR"

    a = read_poscar(pos1)
    b = read_poscar(pos2)
    if a.natoms != b.natoms:
        print(f"error: atom count mismatch ({pos1}: {a.natoms}, "
              f"{pos2}: {b.natoms})", file=sys.stderr)
        return 1

    mode = []
    for i in range(a.natoms):
        d = min_image_delta(a.frac[i], b.frac[i])
        cart = direct_to_cart(d, a.lattice)
        mode.append(cart)

    norm = sum(c**2 for atom in mode for c in atom) ** 0.5
    if norm < 1e-12:
        print("error: structures are identical; no displacement to follow",
              file=sys.stderr)
        return 1
    mode = [[c / norm for c in atom] for atom in mode]

    write_modecar(out, mode)
    print(f"wrote {out}: unit mode from |{pos2} - {pos1}|, "
          f"{a.natoms} atoms (norm {norm:.4f} A)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
