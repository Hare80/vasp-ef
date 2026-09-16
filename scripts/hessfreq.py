#!/usr/bin/env python3
"""hessfreq.py -- mass-weighted frequency analysis of an EFHESSIAN file.

Reads the Cartesian Hessian snapshot written by VASP-EF, mass-weights it
with the atomic masses (from POSCAR+POTCAR, or --masses), and prints the
normal-mode frequencies in cm-1.  Negative (imaginary) frequencies are
printed with a minus sign, VTST/ASE style.  This lets you verify the
"exactly one imaginary mode" criterion offline, without rerunning
IBRION=6/8.

Usage:
    python hessfreq.py [EFHESSIAN] [options]

    --poscar P    structure file (default POSCAR); supplies species/masses
    --masses LIST comma-separated masses in atom order (overrides POTCAR)
    --modes FILE  also write the mass-weighted eigenvectors (default: none)
    --thresh T    |nu| below T cm-1 counts as zero mode (default 10)

Notes:
    * uses numpy when available, otherwise a built-in Jacobi solver
      (fine for a few hundred degrees of freedom)
    * mass-weighting: H_mw = M^-1/2 H M^-1/2; nu = sign*sqrt(|lam|)*521.47
"""
from __future__ import annotations

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vefcommon import read_poscar, title_elements, element_mass  # noqa: E402

WAVENUMBER_CONST = 521.4700        # sqrt(eV/(A^2 amu)) -> cm^-1


def read_efhessian(path):
    """Return (nions, dof, H[n][n])."""
    with open(path) as fh:
        lines = fh.readlines()
    body = []
    nions = dof = None
    for ln in lines:
        s = ln.strip()
        if not s:
            continue
        if s.startswith("#"):
            if "NIONS=" in s:
                nions = int(s.split("NIONS=")[1].split()[0])
            if "DOF=" in s:
                dof = int(s.split("DOF=")[1].split()[0])
            continue
        body.extend(float(t) for t in s.split())
    if dof is None:
        n = math.isqrt(len(body))
        if n * n != len(body):
            raise ValueError(f"{path}: {len(body)} values is not a square matrix")
        dof = n
    if len(body) != dof * dof:
        raise ValueError(f"{path}: expected {dof}x{dof} values, got {len(body)}")
    H = [body[i * dof:(i + 1) * dof] for i in range(dof)]
    return nions, dof, H


def jacobi_eigh(a, sweeps=100, tol=1e-11):
    """Eigen-decomposition of a symmetric matrix by cyclic Jacobi rotations.

    Returns (eigenvalues ascending, eigenvectors as columns v[i][j] =
    component j of mode i).  Pure python; O(n^2) memory, adequate for
    the moderate sizes of TS calculations.
    """
    n = len(a)
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(sweeps):
        off = sum(a[i][j] ** 2 for i in range(n) for j in range(n) if i != j)
        if off < tol:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if abs(a[p][q]) < 1e-300:
                    continue
                theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q])
                t = (1.0 if theta >= 0 else -1.0) / \
                    (abs(theta) + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = v[k][p], v[k][q]
                    v[k][p] = c * vkp - s * vkq
                    v[k][q] = s * vkp + c * vkq
    evals = [a[i][i] for i in range(n)]
    order = sorted(range(n), key=lambda i: evals[i])
    return ([evals[i] for i in order],
            [[v[k][i] for k in range(n)] for i in order])


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("hessian", nargs="?", default="EFHESSIAN")
    ap.add_argument("--poscar", default="POSCAR")
    ap.add_argument("--masses", default=None)
    ap.add_argument("--modes", default=None)
    ap.add_argument("--thresh", type=float, default=10.0)
    args = ap.parse_args(argv)

    if not os.path.isfile(args.hessian):
        print(f"error: {args.hessian} not found", file=sys.stderr)
        return 1
    nions, dof, H = read_efhessian(args.hessian)

    if args.masses:
        masses = [float(m) for m in args.masses.split(",")]
    else:
        syms_atoms = None
        if os.path.isfile(args.poscar):
            p = read_poscar(args.poscar)
            syms_atoms = p.symbols_per_atom()
            if p.natoms != nions:
                print(f"warning: POSCAR has {p.natoms} atoms, EFHESSIAN "
                      f"says {nions}; using POSCAR", file=sys.stderr)
        if syms_atoms is None or any(s == "X" for s in syms_atoms):
            titel = title_elements("POTCAR") if os.path.isfile("POTCAR") else []
            if titel and os.path.isfile(args.poscar):
                counts = read_poscar(args.poscar).counts
                syms_atoms = [e for e, c in zip(titel, counts) for _ in range(c)]
        if syms_atoms is None:
            print("error: cannot determine masses (no POSCAR/POTCAR); "
                  "use --masses", file=sys.stderr)
            return 1
        masses = [element_mass(s) for s in syms_atoms]
    if len(masses) != nions:
        print(f"error: got {len(masses)} masses for {nions} atoms",
              file=sys.stderr)
        return 1

    dminv = [1.0 / math.sqrt(max(m, 1e-8)) for m in masses for _ in range(3)]
    hmw = [[dminv[i] * H[i][j] * dminv[j] for j in range(dof)]
           for i in range(dof)]

    try:
        import numpy as np
        evals, evecs = np.linalg.eigh(np.array(hmw))
        evals = evals.tolist()
        evecs = evecs.T.tolist()
        used = "numpy"
    except ImportError:
        evals, evecs = jacobi_eigh(hmw)
        used = "built-in Jacobi"

    freqs = [(WAVENUMBER_CONST * math.sqrt(abs(l)) * (1 if l >= 0 else -1))
             for l in evals]
    nimag = sum(1 for f in freqs if f < -args.thresh)

    print(f"EFHESSIAN: {nions} atoms, {dof} dof   "
          f"solver: {used}   zero-mode threshold: {args.thresh} cm-1")
    print(f"{'#':>4}  {'nu [cm-1]':>12}")
    for i, f in enumerate(freqs):
        tag = ""
        if f < -args.thresh:
            tag = "   <-- imaginary"
        elif abs(f) <= args.thresh:
            tag = "   (zero mode)"
        print(f"{i:>4}  {f:>12.2f}{tag}")
    print(f"\nimaginary modes (|nu| > {args.thresh} cm-1): {nimag}")
    if nimag == 1:
        print("VERDICT: single imaginary frequency -- consistent with a "
              "transition state")
    elif nimag == 0:
        print("VERDICT: no imaginary frequency -- this is a minimum, not a "
              "saddle point")
    else:
        print(f"VERDICT: {nimag} imaginary frequencies -- higher-order saddle")

    if args.modes:
        with open(args.modes, "w") as fh:
            for i, vec in enumerate(evecs):
                fh.write(f"# mode {i}  nu = {freqs[i]:.2f} cm-1\n")
                for ia in range(nions):
                    comp = [vec[3 * ia + k] for k in range(3)]
                    fh.write("  {:20.10E}  {:20.10E}  {:20.10E}\n".format(*comp))
        print(f"wrote eigenvectors to {args.modes}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
