#!/usr/bin/env python3
"""test_scripts.py -- end-to-end tests for the VASP-EF workflow scripts.

Run from the repository root or anywhere:
    python tests/test_scripts.py

Exit status 0 = all tests passed.  Uses only the standard library.
"""
import math
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "..", "scripts")
sys.path.insert(0, SCRIPTS)

import vefcommon  # noqa: E402
import modemake   # noqa: E402
import efstat     # noqa: E402

FAILURES = []


def check(cond, name):
    print(("  ok   " if cond else "  FAIL ") + name)
    if not cond:
        FAILURES.append(name)


def make_poscar(path, lattice, fracs, symbols, comment="test"):
    with open(path, "w") as fh:
        fh.write(comment + "\n1.0\n")
        for row in lattice:
            fh.write("  %.10f  %.10f  %.10f\n" % tuple(row))
        fh.write(" ".join(symbols) + "\n")
        fh.write(" ".join(str(c) for c in
                          [symbols.count(s) for s in dict.fromkeys(symbols)]) + "\n")
        fh.write("Direct\n")
        for f in fracs:
            fh.write("  %.10f  %.10f  %.10f\n" % tuple(f))


LAT = [[6.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 6.0]]


def test_modemake(tmp):
    p1 = os.path.join(tmp, "POSCAR1")
    p2 = os.path.join(tmp, "POSCAR2")
    mc = os.path.join(tmp, "MODECAR")
    make_poscar(p1, LAT, [[0.0, 0.0, 0.0], [0.5, 0.5, 0.4]], ["H", "O"])
    make_poscar(p2, LAT, [[0.0, 0.0, 0.0], [0.5, 0.5, 0.6]], ["H", "O"])
    rc = subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, "modemake.py"), p1, p2, mc],
        capture_output=True, text=True)
    check(rc.returncode == 0, "modemake: exit 0")
    mode = vefcommon.read_modecar(mc)
    # displacement of atom 2: +0.2*6.0 = 1.2 A along z; atom 1: 0
    check(abs(mode[1][2] - 1.0) < 1e-8, "modemake: unit mode along z")
    check(all(abs(c) < 1e-8 for c in mode[0]), "modemake: atom 1 fixed")
    # identical structures must be rejected
    rc = subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, "modemake.py"), p1, p1, mc],
        capture_output=True, text=True)
    check(rc.returncode != 0, "modemake: rejects identical structures")


def test_efstat(tmp):
    efdat = os.path.join(tmp, "EFDAT")
    with open(efdat, "w") as fh:
        fh.write("#  step         E(eV)        dE(eV)    Fmax(eV/A) "
                 "kappa(eV/A^2)  mode   nu1(cm^-1)  nimag  |step|(A)\n")
        fh.write("   1     -10.000000    1.0000E-01    0.90000    "
                 "-0.5000     1   -123.45     1    0.2000\n")
        fh.write("   2     -10.300000   -3.0000E-01    0.04000    "
                 "-0.6000     1   -152.30     1    0.1500\n")
    with open(os.path.join(tmp, "INCAR"), "w") as fh:
        fh.write("EDIFFG = -0.05\nEF_TS = .TRUE.\n")
    rc = subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, "efstat.py"), tmp],
        capture_output=True, text=True)
    check(rc.returncode == 0, "efstat: exit 0")
    check("CONVERGED" in rc.stdout, "efstat: verdict CONVERGED")
    check("-152.30" in rc.stdout, "efstat: reads nu1")
    # not converged: forces too high
    with open(os.path.join(tmp, "INCAR"), "w") as fh:
        fh.write("EDIFFG = -0.005\n")
    rc = subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, "efstat.py"), tmp],
        capture_output=True, text=True)
    check("not converged yet" in rc.stdout, "efstat: detects unconverged")


def test_hessfreq(tmp):
    # 2 atoms, masses 1 (H) and 16 (O): diagonal blocks k_H=+4, k_O=-1
    hess = os.path.join(tmp, "EFHESSIAN")
    dof = 6
    H = [[0.0] * dof for _ in range(dof)]
    for i in range(3):
        H[i][i] = 4.0
        H[3 + i][3 + i] = -1.0
    with open(hess, "w") as fh:
        fh.write("# EFHESSIAN v1 (VASP-EF)\n# NIONS= 2  DOF= 6\n")
        fh.write("# UNITS= eV/A^2 Cartesian, full symmetric, row-major\n")
        for row in H:
            fh.write("".join(f"{v:24.16E}" for v in row) + "\n")
    poscar = os.path.join(tmp, "POSCAR")
    make_poscar(poscar, LAT, [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], ["H", "O"])
    rc = subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, "hessfreq.py"), hess,
         "--poscar", poscar], capture_output=True, text=True)
    check(rc.returncode == 0, "hessfreq: exit 0")
    check(rc.stderr.strip() == "", "hessfreq: no warnings")
    expected = 521.4700

    def freq_values(stdout):
        vals = []
        for line in stdout.splitlines():
            tok = line.split()
            if len(tok) >= 2:
                try:
                    vals.append(float(tok[1]))
                except ValueError:
                    pass
        return vals
    # mass-weighted: lambda = k/m -> nu = const*sqrt(|k|/m)
    exp_h = expected * math.sqrt(4.0 / 1.008)     # H atom, k = +4
    exp_o = -expected * math.sqrt(1.0 / 15.999)   # O atom, k = -1
    check(any(abs(v - exp_h) < 0.05 for v in freq_values(rc.stdout)),
          "hessfreq: +4 on H -> const*sqrt(4/1.008) cm-1")
    check(any(abs(v - exp_o) < 0.05 for v in freq_values(rc.stdout)),
          "hessfreq: -1 on O -> -const*sqrt(1/15.999) cm-1")
    check(sum(1 for line in rc.stdout.splitlines()
              if "<-- imaginary" in line) == 3, "hessfreq: 3 imaginary modes")
    check("single imaginary" not in rc.stdout, "hessfreq: verdict not TS")


def test_efrestart(tmp):
    src = os.path.join(tmp, "run1")
    dst = os.path.join(tmp, "run2")
    os.makedirs(src)
    make_poscar(os.path.join(src, "CONTCAR"), LAT,
                [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], ["H", "O"])
    with open(os.path.join(src, "NEWMODECAR"), "w") as fh:
        fh.write("  0.0000000000E+00  0.0000000000E+00  1.0000000000E+00\n"
                 "  0.0000000000E+00  0.0000000000E+00 -1.0000000000E+00\n")
    with open(os.path.join(src, "EFHESSIAN"), "w") as fh:
        fh.write("# EFHESSIAN v1\n# NIONS= 2  DOF= 6\n")
    with open(os.path.join(src, "INCAR"), "w") as fh:
        fh.write("IBRION = 3\nPOTIM = 0\nEF_TS = .TRUE.\nEDIFFG = -0.05\n")
    rc = subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, "efrestart.py"), src, dst],
        capture_output=True, text=True)
    check(rc.returncode == 0, "efrestart: exit 0")
    check(os.path.isfile(os.path.join(dst, "POSCAR")), "efrestart: POSCAR placed")
    check(os.path.isfile(os.path.join(dst, "MODECAR")), "efrestart: MODECAR placed")
    check(os.path.isfile(os.path.join(dst, "EFHESSIAN")), "efrestart: EFHESSIAN copied")
    with open(os.path.join(dst, "INCAR")) as fh:
        content = fh.read()
    check("EF_READ_HESS = .TRUE." in content, "efrestart: INCAR gains EF_READ_HESS")
    check("EF_MODE = GUESS" in content, "efrestart: INCAR gains EF_MODE=GUESS")
    # rerunning must not duplicate tags
    subprocess.run([sys.executable, os.path.join(SCRIPTS, "efrestart.py"),
                    src, dst], capture_output=True, text=True)
    with open(os.path.join(dst, "INCAR")) as fh:
        content = fh.read()
    check(content.count("EF_READ_HESS") == 1, "efrestart: idempotent INCAR")


def test_vefcommon(tmp):
    p = os.path.join(tmp, "P4")
    make_poscar(p, LAT, [[0.1, 0.2, 0.3]], ["Si"])
    pos = vefcommon.read_poscar(p)
    cart = vefcommon.direct_to_cart(pos.frac[0], pos.lattice)
    back = vefcommon.cart_to_direct(cart, pos.lattice)
    check(all(abs(a - b) < 1e-10 for a, b in zip(pos.frac[0], back)),
          "vefcommon: direct<->cartesian roundtrip")
    els = ["H", "O"]
    check(abs(vefcommon.element_mass("h") - 1.008) < 1e-9, "vefcommon: mass H")
    check(vefcommon.element_z("si") == 14, "vefcommon: Z Si (case-insensitive)")


def main():
    with tempfile.TemporaryDirectory() as tmp:
        test_vefcommon(tmp)
        test_modemake(tmp)
        test_efstat(tmp)
        test_hessfreq(tmp)
        test_efrestart(tmp)
    if FAILURES:
        print(f"FAILED: {len(FAILURES)} check(s)")
        return 1
    print("ALL SCRIPT TESTS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
