#!/usr/bin/env python3
"""
Coupling probe extended to K₁₀.
Tests: does vacuum null space = 21 pattern generalize?
Does universal σ_max persist? Does degeneracy structure scale?
"""
import numpy as np
from collections import defaultdict
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
    O = np.zeros((na, nb), dtype=np.float64)
    for i in range(na):
        for j in range(nb):
            if same and i == j:
                O[i,j] = nv
            else:
                O[i,j] = 2 * len(ms_a[i] & ms_b[j])
    return O

def compute_dir_perm(p, g=2):
    n_dirs = (p - 1) // 2
    def dir_class(d):
        d = d % p
        return min(d, p - d) - 1
    perm = {}
    for k in range(n_dirs):
        diff = k + 1
        new_diff = (g * diff) % p
        new_dir = dir_class(new_diff)
        if new_dir < 0: return None
        perm[k] = new_dir
    return perm

def get_orbits(blocks, perm):
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
    return orbits

# ======================================================================
def analyze_level(nv, vac_orbit_idx):
    p = nv - 1
    n = nv // 2
    
    print(f"\n{'='*72}")
    print(f"COUPLING PROBE: K_{nv} (p={p})")
    print(f"{'='*72}")
    
    t0 = time.time()
    raw = gen_matchings(nv)
    match_sets = [set(m) for m in raw]
    N = len(raw)
    print(f"  {N} matchings ({time.time()-t0:.1f}s)")
    
    # Direction blocks
    blocks = defaultdict(list)
    for i, m in enumerate(raw):
        blocks[torus_ms(m, p)].append(i)
    
    # Orbits
    perm = compute_dir_perm(p)
    orbits = get_orbits(blocks, perm)
    
    # Collect orbit matchings
    orbit_data = []
    for oms in orbits:
        idx = []
        for ms in oms: idx.extend(blocks[ms])
        orbit_data.append((oms, idx))
        fp = " (fp)" if len(oms)==1 else ""
        print(f"  Orbit {len(orbit_data)-1}: |orb|={len(oms)}, "
              f"ms={len(idx)}, bs={len(blocks[oms[0]])}{fp}")
    
    vac_oms, vac_idx = orbit_data[vac_orbit_idx]
    vac_ms = [match_sets[i] for i in vac_idx]
    bs_total = len(vac_idx)
    
    # --- SENSE 2: Vacuum spectral structure ---
    print(f"\n  VACUUM SPECTRUM (orbit {vac_orbit_idx}, {bs_total} matchings):")
    O_vac = build_overlap(vac_ms, vac_ms, nv, same=True)
    evals = np.sort(np.linalg.eigvalsh(O_vac))[::-1]
    
    # Compress: group by value
    tol = 1e-6
    groups = []
    i = 0
    while i < len(evals):
        val = evals[i]
        count = 1
        while i + count < len(evals) and abs(evals[i+count] - val) < tol:
            count += 1
        groups.append((val, count))
        i += count
    
    print(f"    {'Value':<16} {'Mult':<6} {'Cumul'}")
    print(f"    {'-'*35}")
    cumul = 0
    for val, mult in groups:
        cumul += mult
        print(f"    {val:>14.6f}  {mult:>4}   {cumul}")
    
    rank = sum(m for v,m in groups if abs(v) > tol)
    null_dim = bs_total - rank
    print(f"    Rank: {rank}, Null space: {null_dim}")
    print(f"    C(p,2) = C({p},2) = {p*(p-1)//2}")
    print(f"    bs_total/2 = {bs_total/2}")
    
    # --- SENSE 1: Inter-orbit coupling ---
    print(f"\n  INTER-ORBIT COUPLING:")
    for oi, (oms, oidx) in enumerate(orbit_data):
        if oi == vac_orbit_idx: continue
        other_ms = [match_sets[i] for i in oidx]
        O_cross = build_overlap(vac_ms, other_ms, nv)
        U, sigma, Vt = np.linalg.svd(O_cross, full_matrices=False)
        sigma_nz = sigma[sigma > tol]
        
        # Compress singular values
        sg = []
        j = 0
        while j < len(sigma_nz):
            val = sigma_nz[j]
            cnt = 1
            while j+cnt < len(sigma_nz) and abs(sigma_nz[j+cnt]-val) < tol:
                cnt += 1
            sg.append((val, cnt))
            j += cnt
        
        frob2 = np.sum(sigma**2)
        fp = " (FIXED PT)" if len(oms)==1 else ""
        print(f"\n    Vac ↔ Orbit {oi} (|orb|={len(oms)}, ms={len(oidx)}){fp}")
        print(f"    Rank: {len(sigma_nz)}, ||F||² = {frob2:.1f}")
        for sv, cnt in sg:
            sv2 = sv**2
            print(f"      σ={sv:10.4f} (×{cnt})  σ²={sv2:10.2f}")
    
    print(f"\n  Completed in {time.time()-t0:.2f}s")
    return groups

# ======================================================================
# K₆: vacuum = orbit 0
print("="*72)
print("K₆")
analyze_level(6, 0)

# K₈: vacuum = orbit 2
print("\n" + "="*72)
print("K₈")  
analyze_level(8, 2)

# K₁₀: vacuum = orbit 3
print("\n" + "="*72)
print("K₁₀")
analyze_level(10, 3)
