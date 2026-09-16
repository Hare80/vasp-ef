#!/usr/bin/env bash
# =============================================================================
# install.sh -- install (or remove) VASP-EF in a VASP source tree.
#
# VASP-EF adds an eigenvector-following transition-state optimizer
# (INCAR: EF_TS=.TRUE. with IBRION=3, POTIM=0).  The installer
#
#   * copies src/ef_core.F and src/ef_ts.F into $SRC/src/
#   * patches src/chain.F with a 3-line hook (auto-detects whether the
#     tree already carries the VTST code and uses the matching patch)
#   * registers ef_core.o / ef_ts.o in src/.objects before chain.o
#
# It never touches main.F and contains no VASP source.
#
# Usage:
#   bash install.sh --src /path/to/vasp.x.y.z            # install
#   bash install.sh --src /path/to/vasp.x.y.z --uninstall
#
# Works in bash on Linux and in the MSYS2 UCRT64 shell on Windows.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC=""
UNINSTALL=0

log()  { printf '\033[1;32m[vasp-ef]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[vasp-ef]\033[0m %s\n' "$*"; }
die()  { printf '\033[1;31m[vasp-ef]\033[0m %s\n' "$*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --src)       SRC="${2:-}"; shift 2 ;;
    --uninstall) UNINSTALL=1; shift ;;
    -h|--help)   sed -n '2,20p' "$0"; exit 0 ;;
    *) die "unknown argument: $1 (try --help)" ;;
  esac
done

[[ -n "$SRC" ]] || die "--src <vasp-source-dir> is required"
[[ -f "$SRC/src/main.F" && -f "$SRC/src/.objects" ]] || \
  die "not a VASP source tree: $SRC (missing src/main.F or src/.objects)"
[[ -f "$SRC/src/chain.F" ]] || die "missing $SRC/src/chain.F"

CHAIN_F="$SRC/src/chain.F"
OBJECTS_F="$SRC/src/.objects"
MARKER="# --- VASP-EF (install.sh) ---"

# ---------------------------------------------------------------------------
if [[ "$UNINSTALL" == 1 ]]; then
  if [[ -f "$SRC/src/chain.F.pre_ef" ]]; then
    cp -f "$SRC/src/chain.F.pre_ef" "$CHAIN_F"
    log "restored chain.F from chain.F.pre_ef"
  else
    warn "no chain.F.pre_ef backup found; chain.F left untouched"
  fi
  if grep -qF "$MARKER" "$OBJECTS_F"; then
    tmp="$(mktemp)"
    awk -v m="$MARKER" 'BEGIN{skip=0} $0==m{skip=1;next} skip>0{skip--;next} {print}' \
      "$OBJECTS_F" > "$tmp" && mv -f "$tmp" "$OBJECTS_F"
    log "removed VASP-EF entries from .objects"
  fi
  rm -f "$SRC/src/ef_core.F" "$SRC/src/ef_ts.F"
  log "uninstall complete"
  exit 0
fi

# ---------------------------------------------------------------------------
# idempotency check
if grep -q "ef_ts_step" "$CHAIN_F"; then
  die "chain.F already contains the VASP-EF hook (run --uninstall first?)"
fi
if grep -q "ef_core.o" "$OBJECTS_F"; then
  warn ".objects already lists ef_core.o (continuing)"
fi

# flavor detection: VTST trees call opt_step from chain.F
FLAVOR="stock"
if grep -q "opt_step" "$CHAIN_F"; then
  FLAVOR="vtst"
fi
log "detected chain.F flavor: $FLAVOR"

# 1) copy our sources
cp -f "$HERE/src/ef_core.F" "$HERE/src/ef_ts.F" "$SRC/src/"
log "copied ef_core.F, ef_ts.F -> $SRC/src/"

# 2) patch chain.F (keep a one-time backup)
if [[ ! -f "$SRC/src/chain.F.pre_ef" ]]; then
  cp -f "$CHAIN_F" "$SRC/src/chain.F.pre_ef"
fi
if command -v patch >/dev/null 2>&1; then
  (cd "$SRC" && patch -p1 -N --quiet < "$HERE/patches/chain_F_hook.$FLAVOR.patch") || \
    die "patch failed on src/chain.F (see messages above)"
else
  die "'patch' not found -- install patch (Linux: coreutils distro package, MSYS2: pacman -S patch)"
fi
grep -q "ef_ts_step" "$CHAIN_F" || die "patch verification failed (no ef_ts_step in chain.F)"
log "patched src/chain.F (hook installed)"

# 3) register objects before chain.o (idempotent)
if ! grep -q "ef_core.o" "$OBJECTS_F"; then
  tmp="$(mktemp)"
  awk -v m="$MARKER" '
    !done && /(^|[[:space:]\\])chain\.o([[:space:]\\]|$)/ {
      print m
      print "  ef_core.o ef_ts.o \\"
      done = 1
    }
    { print }
    END { if (!done) exit 2 }
  ' "$OBJECTS_F" > "$tmp" || { rm -f "$tmp"; die "could not find chain.o in .objects"; }
  mv -f "$tmp" "$OBJECTS_F"
fi
awk '/ef_core\.o/{ef=NR} /(^|[[:space:]\\])chain\.o([[:space:]\\]|$)/{c=NR} END{exit !(ef && c && ef<c)}' \
  "$OBJECTS_F" || die ".objects verification failed (ef_core.o must precede chain.o)"
log "registered ef_core.o ef_ts.o in src/.objects"

log "done.  Next steps:"
log "  1. rebuild VASP (make / cmake) in $SRC"
log "  2. run with INCAR:  IBRION=3, POTIM=0, EF_TS=.TRUE."
log "  3. monitor with:    python scripts/efstat.py   (see README.md)"
