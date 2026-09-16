# Installing VASP-EF

VASP-EF adds two files and a small patch to your VASP `src/` directory.
No configure step, no dependency changes, no build system modification.

## Files

| source file | destination |
|---|---|
| `src/ef_ts.F` | `<vasp>/src/ef_ts.F` |
| `patches/chain.F` | replaces `<vasp>/src/chain.F` (backup first) |
| `src/` (`.objects` edit) | add `ef_ts.o \` before `chain.o` |

## Steps

1. Copy `src/ef_ts.F` into your VASP source:

   ```
   cp vasp-ef/src/ef_ts.F <vasp>/src/
   ```

2. Replace `chain.F` with the hooked version from `patches/`:

   ```
   cp vasp-ef/patches/chain.F <vasp>/src/chain.F
   ```

   Two variants are provided:
   - `chain_F_hook.stock.patch` for stock VASP
   - `chain_F_hook.vtst.patch` for VTST-patched VASP

   Or apply manually: add `USE ef_ts` to the top of `MODULE chain`
   and `IF (ef_ts_step(nions,posion,toten,force,a,b,iu6)) RETURN`
   as the first executable line of `SUBROUTINE chain_force`.

3. Add `ef_ts.o \` immediately before `chain.o` in `src/.objects`.

4. Rebuild:

   ```
   make std          # or make gam / make ncl
   ```

   For CMake builds the `.objects` change is picked up automatically.

## INCAR

Add to your INCAR:

```
IBRION = 3
POTIM  = 0
EF_TS  = .TRUE.
```

See `docs/MANUAL.md` for all available tags and workflow guidance.

## Verify

After rebuild, a short test run should print to stdout:

```
EF: VASP-EF eigenvector-following TS optimizer active
```

## Uninstall

Restore the backup:

```
cp <vasp>/src/chain.F.pre_ef <vasp>/src/chain.F
```

Remove `ef_ts.o` from `.objects` and delete `src/ef_ts.F`.

## Notes

* Windows/MSYS2: same steps; the vasp-windows-msys2-build pipeline
  handles the build system automatically after the files are copied.
* The two C files compile with any Fortran compiler that builds VASP.
  No additional libraries or flags are required.
