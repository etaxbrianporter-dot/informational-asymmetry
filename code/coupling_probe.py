#!/usr/bin/env python3
"""
Coupling Probe: Inter-orbit coupling (Sense 1) + Spectral gaps (Sense 2)

Sense 1: SVD of off-diagonal overlap blocks O(vacuum, other_orbit)
  - Singular values encode coupling strengths between vacuum and excitations
  - Ratios may reveal Yukawa hierarchy or mixing structure

Sense 2: Full eigenvalue spectrum within vacuum block
  - Gap ratios, degeneracy patterns, spectral fingerprint

Starting at K₈ (105 matchings, 4 orbits, vacuum = Orbit 2)
"""

import numpy as np
from itertools import combinations
import time

def generate_matchings(n_vertices):
    """Generate all perfect matchings of K_{n_vertices}."""
    verts = list(range(n_vertices))
    matchings = []
    def _build(remaining, current):
        if not remaining:
            matchings.append(frozenset(current))
            return
        first = remaining[0]
        for i in range(1, len(remaining)):
            partner = remaining[i]
            new_remaining = [v for j, v in enumerate(remaining) if j != 0 and j != i]
            _build(new_remaining, current + [(first, partner)])
    _build(verts, [])
    return matchings

def direction_label(matching, p):
    """Assign direction multiset to a matching (matching extractor convention).
    Edge (a,b): d = min((b-a)%p, (a-b)%p), direction = d-1.
    Edges involving vertex >= p are excluded."""
    dirs = []
    for e in sorted(matching):
        a, b = min(e), max(e)
        if a >= p or b >= p:
            continue
        d = (b - a) % p
        d = min(d, p - d)
        dirs.append(d - 1)
    return tuple(sorted(dirs))

def build_overlap(matchings_a, matchings_b, nv, same_block=False):
    """Build overlap matrix between two sets of matchings.
    O[i,j] = 2*|shared edges| for i≠j, nv for i=j (when same_block)."""
    na, nb = len(matchings_a), len(matchings_b)
    O = np.zeros((na, nb), dtype=np.float64)
    for i in range(na):
        for j in range(nb):
            if same_block and i == j:
                O[i, j] = nv
            else:
                shared = len(matchings_a[i] & matchings_b[j])
                O[i, j] = 2 * shared
    return O

# ======================================================================
# K₈ PROBE
# ======================================================================
print("=" * 72)
print("COUPLING PROBE: K₈")
print("=" * 72)

nv = 8
p = 7
n = 4

t0 = time.time()
raw = generate_matchings(nv)
matchings = [frozenset(m) for m in raw]  # already frozensets
# Convert to sets for intersection
match_sets = [set(m) for m in raw]
N = len(matchings)
print(f"  {N} matchings generated ({time.time()-t0:.2f}s)")

# Direction labeling
dir_labels = [direction_label(m, p) for m in matchings]

# Group into blocks
from collections import defaultdict
blocks = defaultdict(list)
for i, dl in enumerate(dir_labels):
    blocks[dl].append(i)

# ×2 mod 7 orbit structure (known from extractor)
# Dir perm: 0→1, 1→2, 2→0
def apply_perm_k8(ms):
    """×2 mod 7 on direction classes: 0→1, 1→2, 2→0"""
    perm = {0: 1, 1: 2, 2: 0}
    return tuple(sorted(perm[d] for d in ms))

# Build orbits
visited = set()
orbits = []
for ms in sorted(blocks.keys()):
    if ms in visited:
        continue
    orbit_ms = []
    cur = ms
    while cur not in visited:
        visited.add(cur)
        orbit_ms.append(cur)
        cur = apply_perm_k8(cur)
    orbits.append(orbit_ms)

print(f"  {len(orbits)} orbits")

# Identify each orbit's matchings
orbit_matchings = []
for oms in orbits:
    idx = []
    for ms in oms:
        idx.extend(blocks[ms])
    orbit_matchings.append(idx)

# Report orbit structure
for oi, (oms, oidx) in enumerate(zip(orbits, orbit_matchings)):
    fp = " (fixed pt)" if len(oms) == 1 else ""
    print(f"  Orbit {oi}: {oms}, |orbit|={len(oms)}, matchings={len(oidx)}{fp}")

# ======================================================================
# SENSE 2: Vacuum block spectral structure
# ======================================================================
print(f"\n{'='*72}")
print("SENSE 2: VACUUM BLOCK SPECTRAL STRUCTURE")
print(f"{'='*72}")

# Vacuum is Orbit 2 at K₈: multisets (0,0,2), (0,1,1), (1,2,2), bs=14
vac_orbit = 2
vac_idx = orbit_matchings[vac_orbit]
vac_matchings = [match_sets[i] for i in vac_idx]
bs = len(vac_idx)

print(f"\n  Vacuum orbit: {orbits[vac_orbit]}")
print(f"  Block size: {bs}")

O_vac = build_overlap(vac_matchings, vac_matchings, nv, same_block=True)
evals_vac = np.sort(np.linalg.eigvalsh(O_vac))[::-1]  # descending

print(f"\n  Full eigenvalue spectrum (descending):")
for i, ev in enumerate(evals_vac):
    print(f"    λ_{i:2d} = {ev:12.6f}")

# Gap analysis
print(f"\n  Spectral gaps:")
for i in range(len(evals_vac)-1):
    gap = evals_vac[i] - evals_vac[i+1]
    ratio = evals_vac[i+1] / evals_vac[i] if evals_vac[i] > 1e-10 else 0
    print(f"    λ_{i}→λ_{i+1}: gap = {gap:8.4f}, ratio λ_{i+1}/λ_{i} = {ratio:.6f}")

# Key ratios
print(f"\n  Key ratios:")
print(f"    λ_0/λ_1 = {evals_vac[0]/evals_vac[1]:.6f}")
if evals_vac[1] > 1e-10:
    print(f"    λ_1/λ_2 = {evals_vac[1]/evals_vac[2]:.6f}")
print(f"    λ_0/Tr = {evals_vac[0]/sum(evals_vac):.6f}")
print(f"    Tr(O_vac) = {sum(evals_vac):.6f} (= bs × nv = {bs*nv})")

# Also do K₆ for comparison
print(f"\n  --- K₆ vacuum for comparison ---")
nv6 = 6; p6 = 5
raw6 = generate_matchings(nv6)
ms6 = [set(m) for m in raw6]
dl6 = [direction_label(frozenset(m), p6) for m in raw6]
blocks6 = defaultdict(list)
for i, d in enumerate(dl6):
    blocks6[d].append(i)
# Vacuum at K₆: orbit containing (0,0) and (1,1), bs=5
vac6_idx = blocks6[(0,0)] + blocks6[(1,1)]
vac6_ms = [ms6[i] for i in vac6_idx]
O6 = build_overlap(vac6_ms, vac6_ms, nv6, same_block=True)
ev6 = np.sort(np.linalg.eigvalsh(O6))[::-1]
print(f"  Eigenvalues: {[f'{e:.6f}' for e in ev6]}")
for i in range(len(ev6)-1):
    r = ev6[i+1]/ev6[i] if ev6[i] > 1e-10 else 0
    print(f"    λ_{i}/λ_{i+1} = {ev6[i]/ev6[i+1]:.6f}  (ratio next/this = {r:.6f})")

# ======================================================================
# SENSE 1: INTER-ORBIT COUPLING
# ======================================================================
print(f"\n{'='*72}")
print("SENSE 1: INTER-ORBIT COUPLING (SVD of off-diagonal blocks)")
print(f"{'='*72}")

for oi in range(len(orbits)):
    if oi == vac_orbit:
        continue
    
    other_idx = orbit_matchings[oi]
    other_matchings = [match_sets[i] for i in other_idx]
    no = len(other_idx)
    
    # Off-diagonal block: O(vacuum, other)
    O_cross = build_overlap(vac_matchings, other_matchings, nv, same_block=False)
    
    # SVD
    U, sigma, Vt = np.linalg.svd(O_cross, full_matrices=False)
    sigma_nz = sigma[sigma > 1e-10]
    
    fp = " (FIXED PT)" if len(orbits[oi]) == 1 else ""
    print(f"\n  Vacuum ↔ Orbit {oi} ({orbits[oi]}, size={no}){fp}")
    print(f"    O_cross shape: {O_cross.shape}")
    print(f"    Rank: {len(sigma_nz)}")
    print(f"    Singular values:")
    for si, sv in enumerate(sigma_nz):
        print(f"      σ_{si:2d} = {sv:12.6f}")
    
    if len(sigma_nz) >= 2:
        print(f"    Ratios (σ_i/σ_0):")
        for si in range(1, len(sigma_nz)):
            print(f"      σ_{si}/σ_0 = {sigma_nz[si]/sigma_nz[0]:.6f}")
    
    if len(sigma_nz) >= 3:
        print(f"    Ratios (consecutive σ_i/σ_{i-1}):")
        for si in range(1, len(sigma_nz)):
            print(f"      σ_{si}/σ_{si-1} = {sigma_nz[si]/sigma_nz[si-1]:.6f}")
    
    # Frobenius norm = sqrt(sum σ²) = coupling strength
    frob = np.sqrt(np.sum(sigma**2))
    print(f"    ||O_cross||_F = {frob:.6f}")
    print(f"    ||O_cross||_F / (bs_vac × bs_other) = {frob/(bs*no):.6f}")

# ======================================================================
# CROSS-ORBIT COUPLING SUMMARY
# ======================================================================
print(f"\n{'='*72}")
print("COUPLING STRENGTH SUMMARY")
print(f"{'='*72}")

strengths = []
for oi in range(len(orbits)):
    if oi == vac_orbit:
        continue
    other_idx = orbit_matchings[oi]
    other_matchings = [match_sets[i] for i in other_idx]
    no = len(other_idx)
    O_cross = build_overlap(vac_matchings, other_matchings, nv, same_block=False)
    U, sigma, Vt = np.linalg.svd(O_cross, full_matrices=False)
    frob = np.sqrt(np.sum(sigma**2))
    sigma_max = sigma[0]
    rank = np.sum(sigma > 1e-10)
    strengths.append((oi, orbits[oi], no, sigma_max, frob, rank, sigma[sigma > 1e-10]))

print(f"\n  {'Orbit':<8} {'|orb|':<6} {'ms':<6} {'σ_max':<12} {'||F||':<12} {'rank':<6} {'σ_max/σ_max(best)'}")
print(f"  {'-'*65}")
best_sigma = max(s[3] for s in strengths)
for oi, oms, no, sm, frob, rank, _ in strengths:
    fp = " *" if len(oms) == 1 else ""
    print(f"  {oi:<8} {len(oms):<6} {no:<6} {sm:<12.4f} {frob:<12.4f} {rank:<6} {sm/best_sigma:.6f}{fp}")

# Hierarchy check: ratios between coupling strengths
print(f"\n  Coupling hierarchy (σ_max ratios):")
sorted_s = sorted(strengths, key=lambda x: x[3], reverse=True)
for i, (oi, oms, no, sm, frob, rank, _) in enumerate(sorted_s):
    if i > 0:
        ratio = sm / sorted_s[0][3]
        print(f"    Orbit {oi} / Orbit {sorted_s[0][0]} = {ratio:.6f} = 1/{1/ratio:.2f}")
    else:
        print(f"    Orbit {oi}: σ_max = {sm:.6f} (strongest)")

print(f"\n  Done in {time.time()-t0:.2f}s")
