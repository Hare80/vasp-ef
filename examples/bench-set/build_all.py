import numpy as np
from ase import Atoms
from ase.io import write

def rot(axis, ang):
    a = np.array(axis, float); a /= np.linalg.norm(a)
    K = np.array([[0,-a[2],a[1]],[a[2],0,-a[0]],[-a[1],a[0],0]])
    return np.eye(3) + np.sin(ang)*K + (1-np.cos(ang))*(K@K)

def bondcheck(at, pairs, name):
    p = at.get_positions()
    out = []
    for (i,j,d) in pairs:
        out.append((name, i, j, round(float(np.linalg.norm(p[i]-p[j])),3), 'want', d))
    print(out)

# ============ 1) CH3NC -> CH3CN (7 atoms, 正交长胞 8x11x16) ============
# 异腈: CH3-N≡C 沿 +z
ch3 = []
cn, nc, chh = 1.17, 1.42, 1.09
mC  = [0,0,0]; mN = [0,0,nc]; mCt = [0,0,nc+cn]
for a in np.radians([0,120,240]):
    ch3.append(mC + chh*np.array([np.sin(1.91)*np.cos(a), np.sin(1.91)*np.sin(a), -np.cos(1.91)]))
iso = Atoms('CH3NC', positions=[mC]+ch3+[mN, mCt], cell=[8,11,16], pbc=True)
iso.center()
iso.write('ch3nc-run/POSCAR', format='vasp', direct=True, vasp5=True, sort=False)
# TS 初猜: N(0,0,0), Ct(0,0,1.18), 甲基摆到侧面 (0,1.913,0.59)
N0 = np.array([0,0,0]); Ct = np.array([0,0,1.18]); M = np.array([0,1.913,0.59])
Hs = [M + chh*np.array([0.30, np.cos(a)*0.95, np.sin(a)*0.95]) for a in np.radians([90,210,330])]
ts = Atoms('CH3NC', positions=[N0]+[M]+[Ct]+Hs, cell=[8,11,16], pbc=True)
ts.center()
ts.write('ch3nc-run/POSCAR.ts', format='vasp', direct=True, vasp5=True, sort=False)
# MODECAR: 甲基从 N 端摆向 TS 位置
d = (np.array(M) - np.array(mC)); d /= np.linalg.norm(d)
md = np.zeros((6,3)); md[0:4] = d
md[4] = [0,0,0]; md[5] = [0,0,0]
with open('ch3nc-run/MODECAR','w') as f:
    for r in md: f.write('  {:20.10E}  {:20.10E}  {:20.10E}\n'.format(*r))
bondcheck(iso, [(0,4,1.42),(4,5,1.17),(0,1,1.09)], 'CH3NC')

# ============ 2) 丙二醛质子转移 (9 原子, 斜三斜胞) ============
O1 = [-1.28, 1.28, 0]; O3 = [1.28, 1.28, 0]          # H 键供体/受体(顶部)
u = np.array([0.82, -0.57, 0])                        # O3->C3 方向
C3 = np.array(O3) + 1.32*u
u2 = np.array([-0.82, -0.57, 0])
C1 = np.array(O1) + 1.32*u2                           # O1->C1
C2 = np.array([0, -1.15, 0])
HO = [0.30, 1.23, 0]                                  # 转移质子(靠 O3, 井)
HC2 = [1.12, -2.24, 0]
HC1 = [-0.41, -0.38, 0]
HC3 = [3.44, 0.81, 0]
mal = Atoms(symbols=['C','C','C','O','O','H','H','H','H'],
            positions=[C1,C2,C3,O1,O3,HO,HC2,HC1,HC3])
mal.set_cell([[10,0,0],[-2.5,9.5,0],[-1.2,-0.8,15]]); mal.center()
mal.write('mal-run/POSCAR', format='vasp', direct=True, vasp5=True, sort=False)
# MODECAR: 质子从 O3 侧移向中央(TS)
md = np.zeros((9,3)); md[5] = [0.7, -0.3, 0]
md /= np.linalg.norm(md)
with open('mal-run/MODECAR','w') as f:
    for r in md: f.write('  {:20.10E}  {:20.10E}  {:20.10E}\n'.format(*r))

# ============ 3) 环丁烯开环 (10 原子, 另一形状斜胞) ============
C1 = [-0.67,0,0.12]; C2 = [0.67,0,0.12]; C3 = [0.75,1.40,-0.28]; C4 = [-0.75,1.40,-0.28]
Hv1 = [-1.42,-0.45,0.62]; Hv2 = [1.42,-0.45,0.62]
H3a = [0.85,2.15,0.44]; H3b = [0.70,1.65,-1.30]
H4a = [-0.85,2.15,0.44]; H4b = [-0.70,1.65,-1.30]
cyc = Atoms(symbols=['C','C','C','C','H','H','H','H','H','H'],
            positions=[C1,C2,C3,C4,Hv1,Hv2,H3a,H3b,H4a,H4b])
cyc.set_cell([[11,0,0],[-3,10.5,0],[0.8,-1.2,14]]); cyc.center()
cyc.write('cyc-run/POSCAR', format='vasp', direct=True, vasp5=True, sort=False)
# MODECAR: C3-C4 键断开 + 双 H 外旋
md = np.zeros((10,3))
md[2] = [0.35, 0.45, -0.30]; md[3] = [-0.35, 0.45, -0.30]
md[6] = [0.35, 0.45, -0.30]; md[7] = [0.35, 0.45, -0.30]
md[8] = [-0.35, 0.45, -0.30]; md[9] = [-0.35, 0.45, -0.30]
md /= np.linalg.norm(md)
with open('cyc-run/MODECAR','w') as f:
    for r in md: f.write('  {:20.10E}  {:20.10E}  {:20.10E}\n'.format(*r))
bondcheck(cyc, [(0,1,1.34),(1,2,1.51),(2,3,1.50),(3,0,1.51)], 'cyclobutene')

# ============ 4) 鲁棒性冒烟: NH3 与 H3 在斜三斜胞 ============
rho, h = 0.9421, 0.3829
pos = [[0,0,0],[rho,0,-h],[rho*np.cos(2*np.pi/3), rho*np.sin(2*np.pi/3), -h],
       [rho*np.cos(4*np.pi/3), rho*np.sin(4*np.pi/3), -h]]
nh3 = Atoms('NH3', positions=pos, cell=[[9,0,0],[1.5,8.5,0],[-1.0,0.6,11]], pbc=True)
nh3.center()
nh3.write('nh3-tri/POSCAR', format='vasp', direct=True, vasp5=True, sort=False)
md = np.zeros((4,3)); md[0] = [0,0,1.0]; md[1:4] = [0,0,-0.12]
md /= np.linalg.norm(md)
with open('nh3-tri/MODECAR','w') as f:
    for r in md: f.write('  {:20.10E}  {:20.10E}  {:20.10E}\n'.format(*r))
h3 = Atoms('H3', positions=[[0,0,-0.85],[0,0,0.0],[0,0,1.0]],
           cell=[[9,0,0],[-1.5,8,0],[0.5,0.8,10]], pbc=True)
h3.center()
h3.write('h3-tri/POSCAR', format='vasp', direct=True, vasp5=True, sort=False)
md = np.zeros((3,3)); md[0] = [0,0,-0.37]; md[2] = [0,0,0.93]
md /= np.linalg.norm(md)
with open('h3-tri/MODECAR','w') as f:
    for r in md: f.write('  {:20.10E}  {:20.10E}  {:20.10E}\n'.format(*r))
print('ALL BUILT')
