# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""K12 validation of three coupling laws."""
import numpy as np
from collections import defaultdict
from fractions import Fraction
import time

def gen_matchings(nv):
    verts = list(range(nv))
    matchings = []
    def _build(rem, cur):
        if not rem:
            matchings.append(frozenset(cur))
            return
        first = rem[0]
        for i in range(1, len(rem)):
            _build([v for j,v in enumerate(rem) if j!=0 and j!=i],
                   cur + [(first, rem[i])])
    _build(verts, [])
    return matchings

def edge_dir(a, b, p):
    a, b = min(a,b), max(a,b)
    if a >= p or b >= p: return -1
    d = (b - a) % p
    return min(d, p - d) - 1

def torus_ms(matching, p):
    dirs = []
    for e in sorted(matching):
        d = edge_dir(e[0], e[1], p)
        if d >= 0: dirs.append(d)
    return tuple(sorted(dirs))

def build_overlap(ms_a, ms_b, nv, same=False):
    na, nb = len(ms_a), len(ms_b)
    O = np.zeros((na, nb))
    for i in range(na):
        for j in range(nb):
            O[i,j] = nv if (same and i==j) else 2*len(ms_a[i] & ms_b[j])
    return O

def compute_dir_perm(p, g=2):
    n_dirs = (p-1)//2
    def dc(d):
        d = d%p
        return min(d, p-d) - 1
    perm = {}
    for k in range(n_dirs):
        nd = dc((g*(k+1)) % p)
        if nd < 0: return None
        perm[k] = nd
    return perm

nv = 12
p = 11
print(f"K_12 LAW VALIDATION (nv={nv}, p={p})")
print("="*60)

t0 = time.time()
raw = gen_matchings(nv)
match_sets = [set(m) for m in raw]
N = len(raw)
print(f"  {N} matchings ({time.time()-t0:.1f}s)")

blocks = defaultdict(list)
for i, m in enumerate(raw):
    blocks[torus_ms(m, p)].append(i)
print(f"  {len(blocks)} direction blocks")

perm = compute_dir_perm(p)
visited = set()
orbits = []
for ms in sorted(blocks.keys()):
    if ms in visited: continue
    orbit_ms = []
    cur = ms
    while cur not in visited:
        visited.add(cur)
        orbit_ms.append(cur)
        cur = tuple(sorted(perm[d] for d in cur))
    orbits.append(orbit_ms)

orbit_data = []
fp_idx = None
for oi, oms in enumerate(orbits):
    idx = []
    for ms in oms: idx.extend(blocks[ms])
    orbit_data.append((oms, idx))
    fp = ""
    if len(oms) == 1:
        fp = " (FIXED POINT)"
        fp_idx = oi
    print(f"  Orbit {oi}: |orb|={len(oms)}, ms={len(idx)}, bs={len(blocks[oms[0]])}{fp}")

# Find vacuum orbit by bs=2p pattern
vac_orbit_idx = None
for oi, (oms, idx) in enumerate(orbit_data):
    if len(blocks[oms[0]]) == 2*p:
        vac_orbit_idx = oi
        print(f"\n  Vacuum: Orbit {oi} (bs_block=2p={2*p})")
        break

vac_oms, vac_idx = orbit_data[vac_orbit_idx]
vac_ms_list = [match_sets[i] for i in vac_idx]
bs_total = len(vac_idx)
print(f"  bs_total = {bs_total}")

# LAW 1
print(f"\n{'='*60}")
print("LAW 1: VACUUM HALF-RANK")
print(f"{'='*60}")
t1 = time.time()
O_vac = build_overlap(vac_ms_list, vac_ms_list, nv, same=True)
evals = np.linalg.eigvalsh(O_vac)
rank = np.sum(np.abs(evals) > 1e-6)
null_dim = bs_total - rank
print(f"  rank = {rank}, null = {null_dim}, bs/2 = {bs_total//2}")
law1 = null_dim == bs_total//2
print(f"  LAW 1: {'CONFIRMED' if law1 else 'VIOLATED'} ({time.time()-t1:.1f}s)")

# Spectrum
evals_s = np.sort(evals)[::-1]
groups = []
i = 0
while i < len(evals_s):
    v = evals_s[i]; c = 1
    while i+c < len(evals_s) and abs(evals_s[i+c]-v) < 1e-6: c += 1
    groups.append((v, c)); i += c
print(f"\n  Spectrum:")
for v, c in groups:
    print(f"    {v:12.4f} (x{c})")

# Null projector diagonal
evals_f, evecs_f = np.linalg.eigh(O_vac)
null_v = evecs_f[:, np.abs(evals_f) < 1e-8]
proj = null_v @ null_v.T
diag = np.diag(proj)
print(f"\n  Null projector diagonal: min={diag.min():.6f} max={diag.max():.6f}")
print(f"  All exactly 1/2? {np.allclose(diag, 0.5)}")

# LAW 2
print(f"\n{'='*60}")
print("LAW 2: FP RANK DEFICIT")
print(f"{'='*60}")
if fp_idx is not None:
    fp_oms, fp_list = orbit_data[fp_idx]
    fp_ms = [match_sets[i] for i in fp_list]
    O_cross = build_overlap(vac_ms_list, fp_ms, nv)
    _, sigma, _ = np.linalg.svd(O_cross, full_matrices=False)
    fp_rank = np.sum(sigma > 1e-6)
    deficit = rank - fp_rank
    print(f"  Vacuum rank: {rank}, FP rank: {fp_rank}, deficit: {deficit}")
    law2 = deficit == 2
    print(f"  LAW 2: {'CONFIRMED' if law2 else 'VIOLATED'}")

# LAW 3
print(f"\n{'='*60}")
print("LAW 3: p-RATIONAL COUPLING")
print(f"{'='*60}")
for oi in [0, 1, fp_idx]:
    if oi is None or oi == vac_orbit_idx: continue
    oms, oidx = orbit_data[oi]
    other_ms = [match_sets[i] for i in oidx]
    O_c = build_overlap(vac_ms_list, other_ms, nv)
    _, sig, _ = np.linalg.svd(O_c, full_matrices=False)
    f2 = round(np.sum(sig**2))
    ratio = Fraction(f2, bs_total * len(oidx))
    d = ratio.denominator
    while d % p == 0: d //= p
    fp_str = " (FP)" if len(oms)==1 else ""
    print(f"  Orb {oi}{fp_str}: F^2={f2}, ratio={ratio}, p-rational: {d==1}")

print(f"\nTotal: {time.time()-t0:.1f}s")
