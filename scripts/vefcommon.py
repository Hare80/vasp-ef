#!/usr/bin/env python3
"""vefcommon.py -- shared helpers for the VASP-EF workflow scripts.

Pure Python 3 standard library.  Provides:

  read_poscar(path)      -> Poscar(lattice, frac, symbols, counts)
  title_elements(path)   -> list of element symbols from POTCAR TITEL lines
  element_mass(sym)      -> standard atomic mass in amu
  element_z(sym)         -> atomic number
  direct_to_cart / cart_to_direct / min_image

Element tables cover Z = 1..60; unknown symbols fall back to carbon-like
values with a warning, mirroring the Fortran side of VASP-EF.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field
from typing import List, Optional

# symbol, atomic number, standard atomic mass (amu)
_TABLE = [
    ("H", 1, 1.008), ("He", 2, 4.003), ("Li", 3, 6.940), ("Be", 4, 9.012),
    ("B", 5, 10.810), ("C", 6, 12.011), ("N", 7, 14.007), ("O", 8, 15.999),
    ("F", 9, 18.998), ("Ne", 10, 20.180), ("Na", 11, 22.990),
    ("Mg", 12, 24.305), ("Al", 13, 26.982), ("Si", 14, 28.085),
    ("P", 15, 30.974), ("S", 16, 32.060), ("Cl", 17, 35.450),
    ("Ar", 18, 39.948), ("K", 19, 39.098), ("Ca", 20, 40.078),
    ("Sc", 21, 44.956), ("Ti", 22, 47.867), ("V", 23, 50.942),
    ("Cr", 24, 51.996), ("Mn", 25, 54.938), ("Fe", 26, 55.845),
    ("Co", 27, 58.933), ("Ni", 28, 58.693), ("Cu", 29, 63.546),
    ("Zn", 30, 65.380), ("Ga", 31, 69.723), ("Ge", 32, 72.630),
    ("As", 33, 74.922), ("Se", 34, 78.971), ("Br", 35, 79.904),
    ("Kr", 36, 83.798), ("Rb", 37, 85.468), ("Sr", 38, 87.620),
    ("Y", 39, 88.906), ("Zr", 40, 91.224), ("Nb", 41, 92.906),
    ("Mo", 42, 95.950), ("Tc", 43, 98.000), ("Ru", 44, 101.070),
    ("Rh", 45, 102.910), ("Pd", 46, 106.420), ("Ag", 47, 107.870),
    ("Cd", 48, 112.410), ("In", 49, 114.820), ("Sn", 50, 118.710),
    ("Sb", 51, 121.760), ("Te", 52, 127.600), ("I", 53, 126.900),
    ("Xe", 54, 131.290), ("Cs", 55, 132.910), ("Ba", 56, 137.330),
    ("La", 57, 138.910), ("Ce", 58, 140.120), ("Pr", 59, 140.910),
    ("Nd", 60, 144.240),
]
_MASS = {s: m for s, _, m in _TABLE}
_Z = {s: z for s, z, _ in _TABLE}

_WARNED = set()


def element_mass(sym: str) -> float:
    s = sym.strip().capitalize()
    if s not in _MASS:
        if s not in _WARNED:
            print(f"warning: unknown element '{sym}', using C-like defaults",
                  file=sys.stderr)
            _WARNED.add(s)
        return 12.011
    return _MASS[s]


def element_z(sym: str) -> int:
    s = sym.strip().capitalize()
    if s not in _Z:
        return 6
    return _Z[s]


@dataclass
class Poscar:
    comment: str
    lattice: List[List[float]]              # 3x3, rows are lattice vectors (A)
    frac: List[List[float]]                 # n x 3 direct coordinates
    symbols: List[str] = field(default_factory=list)  # per species
    counts: List[int] = field(default_factory=list)   # per species

    @property
    def natoms(self) -> int:
        return len(self.frac)

    def symbols_per_atom(self) -> List[str]:
        out = []
        for s, c in zip(self.symbols, self.counts):
            out += [s] * c
        if len(out) != self.natoms:            # v4 POSCAR without symbols
            out = ["X"] * self.natoms
        return out

    def masses(self) -> List[float]:
        return [element_mass(s) for s in self.symbols_per_atom()]


def _floats(tokens):
    return [float(t) for t in tokens]


def read_poscar(path: str = "POSCAR") -> Poscar:
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    comment = lines[0]
    scale = float(lines[1].split()[0])
    lat = [_floats(lines[2 + i].split()[:3]) for i in range(3)]
    if scale < 0.0:
        # negative scale = target volume: scale = (V/|scale|)^(1/3)
        import math
        vol = (lat[0][0]*lat[1][1]*lat[2][2]
               + lat[0][1]*lat[1][2]*lat[2][0]
               + lat[0][2]*lat[1][0]*lat[2][1]
               - lat[0][2]*lat[1][1]*lat[2][0]
               - lat[0][1]*lat[1][0]*lat[2][2]
               - lat[0][0]*lat[1][2]*lat[2][1])
        scale = abs(scale / vol) ** (1.0 / 3.0)
    lattice = [[c * scale for c in row] for row in lat]

    counts: List[int] = []
    symbols: List[str] = []
    idx = 6
    tok6 = lines[5].split()
    if tok6 and all(t.isalpha() for t in tok6):
        symbols = tok6
        counts = [int(t) for t in lines[6].split()[:len(symbols)]]
        idx = 7
    else:
        counts = [int(t) for t in tok6]

    ln = lines[idx].strip().lower()
    if ln.startswith("s"):
        idx += 1                                # Selective dynamics line
    idx += 1                                    # Direct/Cartesian line
    coord_type = lines[idx - 1].strip().lower()[:1]

    frac: List[List[float]] = []
    for i in range(sum(counts)):
        toks = lines[idx + i].split()[:3]
        if coord_type == "c" or coord_type == "k":
            cart = _floats(toks)
            frac.append(cart_to_direct(cart, lattice))
        else:
            frac.append(_floats(toks))
    return Poscar(comment, lattice, frac, symbols, counts)


def title_elements(path: str = "POTCAR") -> List[str]:
    """Element symbol of each TITEL block in a (concatenated) POTCAR."""
    out: List[str] = []
    try:
        with open(path, errors="replace") as fh:
            for line in fh:
                if "TITEL" in line and "=" in line:
                    rhs = line.split("=", 1)[1].split()
                    if len(rhs) >= 2:
                        out.append(rhs[1])
    except OSError:
        pass
    return out


def direct_to_cart(direct, lattice) -> List[float]:
    """x_cart(alpha) = sum_d x_d * lattice[d][alpha]  (VASP convention)."""
    return [sum(direct[d] * lattice[d][a] for d in range(3))
            for a in range(3)]


def cart_to_direct(cart, lattice) -> List[float]:
    """Solve x_cart = x_direct . lattice for x_direct (3x3 Cramer)."""
    import copy
    a = copy.deepcopy(lattice)
    rhs = list(cart)
    # Gaussian elimination with partial pivoting: a^T x = rhs
    for i in range(3):
        a[i] = [a[j][i] for j in range(3)]
    m = [a[i][:] + [rhs[i]] for i in range(3)]
    for col in range(3):
        piv = max(range(col, 3), key=lambda r: abs(m[r][col]))
        if abs(m[piv][col]) < 1e-300:
            raise ValueError("singular lattice")
        m[col], m[piv] = m[piv], m[col]
        for r in range(3):
            if r != col and abs(m[col][col]) > 0:
                f = m[r][col] / m[col][col]
                m[r] = [m[r][k] - f * m[col][k] for k in range(4)]
    return [m[i][3] / m[i][i] for i in range(3)]


def min_image_delta(direct_i, direct_j) -> List[float]:
    """Minimum-image difference of two direct-coordinate vectors."""
    return [(direct_j[d] - direct_i[d]) - round(direct_j[d] - direct_i[d])
            for d in range(3)]


def read_modecar(path: str = "MODECAR"):
    """Return the (n,3) mode array as a flat list of per-atom [x,y,z]."""
    with open(path) as fh:
        tokens = fh.readline().split() + \
            [t for ln in fh for t in ln.split()]
    vals = [float(t) for t in tokens]
    if len(vals) % 3 != 0:
        raise ValueError(f"{path}: component count not a multiple of 3")
    return [vals[i:i + 3] for i in range(0, len(vals), 3)]


def write_modecar(path: str, mode_flat) -> None:
    """mode_flat: list of per-atom [x,y,z] (written vtst-style, 3 per line)."""
    with open(path, "w") as fh:
        for atom in mode_flat:
            fh.write("  {:20.10E}  {:20.10E}  {:20.10E}\n".format(*atom))


if __name__ == "__main__":                     # tiny self-test
    p = read_poscar(sys.argv[1] if len(sys.argv) > 1 else "POSCAR")
    print(p.comment, p.natoms, "atoms", p.symbols, p.counts)
