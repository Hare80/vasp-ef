#!/usr/bin/env python3
"""make_slab.py -- build the H/Si(001) H-diffusion TS-search inputs.

Creates a dimer-reconstructed Si(001)-(2x1) slab (4 layers, H-terminated
bottom), one H atom on one dimer Si, and writes:

  POSCAR.ts      H at the bridging position between the two dimer Si
                 (starting guess for the saddle)
  MODECAR.ts     unit mode: H moving vertically (the transfer coordinate)

Requires ASE.  The second dimer in the cell is a spectator.
"""
import numpy as np
from ase.build import diamond100

slab = diamond100('Si', size=(2, 2, 4), a=5.4309, vacuum=9)
cell = np.array(slab.cell)
top_z = slab.positions[:, 2].max()
top = [a.index for a in slab if abs(a.position[2] - top_z) < 1e-6]

# pair the 4 top Si into 2 dimers along x (same y, |dx| = half the cell)
pairs = []
used = set()
for i in top:
    if i in used:
        continue
    for j in top:
        if j == i or j in used:
            continue
        d = slab.positions[j] - slab.positions[i]
        if abs(d[1]) < 0.1 and abs(abs(d[0]) - slab.cell[0, 0] / 2) < 0.6:
            pairs.append((i, j))
            used.update((i, j))
            break
assert len(pairs) == 2, f"expected 2 dimer pairs, got {pairs}"

# reconstruct: pull each pair together (Si-Si dimer ~2.4 A) and down
for i, j in pairs:
    pi, pj = slab.positions[i].copy(), slab.positions[j].copy()
    mid = (pi + pj) / 2
    for idx, p, other in ((i, pi, pj), (j, pj, pi)):
        direction = mid - p
        direction[2] = 0.0
        slab.positions[idx] = p + direction * 0.28      # ~0.55 A each
        slab.positions[idx][2] -= 0.12

# dimer 0 carries the migrating H; dimer 1 gets H termination (passivate)
i, j = pairs[0]
si_a, si_b = slab.positions[i].copy(), slab.positions[j].copy()
d = si_b - si_a
d /= np.linalg.norm(d)

# monohydride H on the OUTER Si of dimer 0 (along its dangling bond: up+outward)
out = np.array([d[0], 0.0, 1.0])
out /= np.linalg.norm(out)
h_start = si_a + out * 1.49

# bottom passivation: not done (fixed-cell demo; bottom Si keep dangling)
# --- write POSCAR.ts: H at the bridge above the dimer centre ---
bridge = (si_a + si_b) / 2 + np.array([0.0, 0.0, 1.55])
at = slab.copy()
at.append('H')
at.positions[-1] = bridge
at.write('POSCAR.ts', format='vasp', direct=True, vasp5=True, sort=True)

# --- MODECAR.ts: unit mode = H vertical motion (x aligned to cell) ---
mode = np.zeros((len(at), 3))
mode[-1] = [0.0, 0.0, 1.0]
with open('MODECAR.ts', 'w') as fh:
    for row in mode:
        fh.write('  {:20.10E}  {:20.10E}  {:20.10E}\n'.format(*row))

print('top pairs:', pairs)
print('H start (dimer Si):', h_start)
print('H bridge (TS guess):', bridge)
print('wrote POSCAR.ts ({} atoms) and MODECAR.ts'.format(len(at)))
