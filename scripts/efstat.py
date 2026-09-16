#!/usr/bin/env python3
"""efstat.py -- print the status of a VASP-EF transition-state run.

The counterpart of vtstscripts' nebstat.pl / dimstat.pl.  Reads EFDAT
(written by VASP-EF every ionic step) and INCAR (for EDIFFG), then
prints the step history and a verdict.

Usage:
    python efstat.py [run-dir] [-n LAST]

    run-dir     directory holding EFDAT/INCAR  (default: .)
    -n LAST     show only the last LAST steps
"""
from __future__ import annotations

import argparse
import os
import sys

HEADER = ("step", "E (eV)", "dE (eV)", "Fmax", "kappa", "mode",
          "nu1 (cm-1)", "nimag", "|step|")


def read_efdat(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tok = line.split()
            if len(tok) < 9:
                continue
            rows.append(dict(step=int(tok[0]), e=float(tok[1]),
                             de=float(tok[2]), fmax=float(tok[3]),
                             kappa=float(tok[4]), mode=int(tok[5]),
                             nu1=float(tok[6]), nimag=int(tok[7]),
                             step_len=float(tok[8])))
    return rows


def read_ediffg(path):
    if not os.path.isfile(path):
        return None
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.split("!")[0].split("#")[0].split(";")[0]
            if "=" not in line:
                continue
            key, _, val = line.partition("=")
            if key.strip().upper() == "EDIFFG":
                try:
                    return float(val.split()[0])
                except (ValueError, IndexError):
                    return None
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundir", nargs="?", default=".")
    ap.add_argument("-n", type=int, default=None, help="show last N steps")
    args = ap.parse_args(argv)

    efdat = os.path.join(args.rundir, "EFDAT")
    if not os.path.isfile(efdat):
        print(f"error: {efdat} not found (is this a VASP-EF run directory?)",
              file=sys.stderr)
        return 1
    rows = read_efdat(efdat)
    if not rows:
        print("EFDAT is empty")
        return 0

    show = rows if args.n is None else rows[-args.n:]
    print(f"{'':>6}  {'E (eV)':>14}  {'dE (eV)':>11}  {'Fmax':>9}  "
          f"{'kappa':>9}  {'mode':>4}  {'nu1 (cm-1)':>11}  {'nimag':>5}  "
          f"{'|step|':>7}")
    for r in show:
        print(f"{r['step']:>6}  {r['e']:>14.6f}  {r['de']:>11.3E}  "
              f"{r['fmax']:>9.5f}  {r['kappa']:>9.4f}  {r['mode']:>4}  "
              f"{r['nu1']:>11.2f}  {r['nimag']:>5}  {r['step_len']:>7.4f}")

    last = rows[-1]
    ediffg = read_ediffg(os.path.join(args.rundir, "INCAR"))
    print("")
    print(f"steps: {last['step']}   energy: {last['e']:.6f} eV   "
          f"Fmax: {last['fmax']:.5f} eV/A   imaginary modes: {last['nimag']}")
    force_ok = (ediffg is not None and last["fmax"] < abs(ediffg))
    if force_ok:
        print(f"force criterion: MET (Fmax < |EDIFFG| = {abs(ediffg):.5f})")
    else:
        print("force criterion: not met" +
              ("" if ediffg is not None else " (EDIFFG not found in INCAR)"))
    if force_ok and (last["nimag"] == 1):
        print("VERDICT: CONVERGED -- stationary point with a single "
              "imaginary frequency (a transition state)")
    elif force_ok:
        print("VERDICT: force-converged, but not a 1-mode saddle; consider "
              "EF_MODE=TRACK/GUESS or restarting from a dimer run")
    else:
        print("VERDICT: not converged yet")
    return 0


if __name__ == "__main__":
    sys.exit(main())
