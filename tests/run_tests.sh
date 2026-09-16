#!/usr/bin/env bash
# run_tests.sh -- build and run the standalone ef_core unit tests.
#
# Requires: gfortran, a BLAS/LAPACK implementation.
#   Linux (Debian/Ubuntu):  sudo apt install gfortran liblapack-dev
#   MSYS2 UCRT64:           pacman -S mingw-w64-ucrt-x86_64-gcc-fortran \
#                                 mingw-w64-ucrt-x86_64-openblas
#
# Usage:  bash tests/run_tests.sh
set -euo pipefail
cd "$(dirname "$0")/build"

FC="${FC:-gfortran}"
# OpenBLAS bundles LAPACK; fall back to a separate liblapack if needed.
LIBS="${LIBS:--lopenblas}"

rm -f *.o *.mod test_ef_core test_ef_core.exe
"$FC" -O2 -ffree-form -ffree-line-length-none -Wall -c ../../src/ef_core.F
"$FC" -O2 -Wall -c ../test_ef_core.F90
"$FC" -o test_ef_core ef_core.o test_ef_core.o ${LIBS}
./test_ef_core
