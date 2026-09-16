#!/usr/bin/env python3
"""efrestart.py -- prepare a restarted VASP-EF run from a stopped one.

Copies the continuation state of a previous (interrupted or force-
converged-but-not-saddle) run into the current directory, the vtst way:

    CONTCAR     -> POSCAR      (previous POSCAR kept as POSCAR.pre_ef)
    NEWMODECAR  -> MODECAR     (last followed mode becomes the guess)
    EFHESSIAN   -> EFHESSIAN   (kept; INCAR gains EF_READ_HESS=.TRUE.)

Usage:
    python efrestart.py [src-dir] [dst-dir]     (defaults: . .)

If dst differs from src, files are copied (src is left untouched);
if they are the same, in-place backups are made.  INCAR handling:
EF_READ_HESS is appended to INCAR in dst if not already present.
"""
from __future__ import annotations

import os
import shutil
import sys


def has_tag(incar, tag):
    if not os.path.isfile(incar):
        return True                        # nothing to do, don't touch
    with open(incar, errors="replace") as fh:
        for line in fh:
            line = line.split("!")[0].split("#")[0].split(";")[0]
            if "=" in line and line.split("=")[0].strip().upper() == tag:
                return True
    return False


def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    src = args[0] if len(args) > 0 else "."
    dst = args[1] if len(args) > 1 else src

    need = ["CONTCAR", "NEWMODECAR"]
    for f in need:
        if not os.path.isfile(os.path.join(src, f)):
            print(f"error: {os.path.join(src, f)} not found -- nothing to "
                  "restart from", file=sys.stderr)
            return 1

    os.makedirs(dst, exist_ok=True)
    same = os.path.abspath(src) == os.path.abspath(dst)

    def transfer(name, dst_name=None, backup=True):
        dst_name = dst_name or name
        s = os.path.join(src, name)
        d = os.path.join(dst, dst_name)
        if same and backup and os.path.isfile(d):
            shutil.copy2(d, d + ".pre_ef")
        shutil.copy2(s, d)
        print(f"  {name} -> {d}")

    transfer("CONTCAR", "POSCAR")
    transfer("NEWMODECAR", "MODECAR")
    if os.path.isfile(os.path.join(src, "EFHESSIAN")):
        transfer("EFHESSIAN")
    else:
        print("  (no EFHESSIAN in src; the restart will start from the "
              "model Hessian)")

    incar_src = os.path.join(src, "INCAR")
    incar_dst = os.path.join(dst, "INCAR")
    if os.path.isfile(incar_src):
        if same:
            shutil.copy2(incar_dst, incar_dst + ".pre_ef")
        shutil.copy2(incar_src, incar_dst)
    if os.path.isfile(incar_dst) and not has_tag(incar_dst, "EF_READ_HESS"):
        with open(incar_dst, "a") as fh:
            fh.write("\nEF_READ_HESS = .TRUE.\n")
        print("  INCAR += EF_READ_HESS=.TRUE.")
    if os.path.isfile(incar_dst) and not has_tag(incar_dst, "EF_MODE"):
        with open(incar_dst, "a") as fh:
            fh.write("EF_MODE = GUESS\n")
        print("  INCAR += EF_MODE=GUESS (follow the MODECAR you just installed)")

    print("restart directory ready.  Run VASP again from there.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
