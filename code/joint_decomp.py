#!/usr/bin/env python3
"""Compute exact joint C₂×C₃ decomposition of null and live spaces at K₈."""
import numpy as np
from collections import defaultdict

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

nv, p = 8, 7
raw = gen_matchings(nv)
match_sets = [set(m) for m in raw]

blocks = defaultdict(list)
for i, m in enumerate(raw):
    blocks[torus_ms(m, p)].append(i)

vac_blocks = [(0,0,2), (0,1,1), (1,2,2)]
vac_idx = []
for vb in vac_blocks:
    vac_idx.extend(blocks[vb])
bs = len(vac_idx)

# Overlap matrix
O = np.zeros((bs, bs))
for i in range(bs):
    for j in range(bs):
        O[i,j] = nv if i==j else 2*len(match_sets[vac_idx[i]] & match_sets[vac_idx[j]])

evals, evecs = np.linalg.eigh(O)
order = np.argsort(evals)[::-1]
evals, evecs = evals[order], evecs[:, order]
null_vecs = evecs[:, np.abs(evals) < 1e-8]
live_vecs = evecs[:, np.abs(evals) >= 1e-8]

# Build P (×2 mod 7) and Q (×6 mod 7)
match_to_vi = {raw[vac_idx[vi]]: vi for vi in range(bs)}

def apply_mult(matching, g, p):
    new = set()
    for a,b in matching:
        na = (g*a)%p if a<p else a
        nb = (g*b)%p if b<p else b
        new.add((min(na,nb), max(na,nb)))
    return frozenset(new)

def perm_matrix(g):
    M = np.zeros((bs, bs))
    for vi in range(bs):
        m_new = apply_mult(raw[vac_idx[vi]], g, p)
        M[match_to_vi[m_new], vi] = 1.0
    return M

P = perm_matrix(2)  # C₃ generator
Q = perm_matrix(6)  # C₂ generator (= ×(-1) mod 7)

# Joint eigenspace decomposition
# P has eigenvalues {1, ω, ω²}, Q has eigenvalues {+1, -1}
omega = np.exp(2j*np.pi/3)

# For each (C₂ sign, C₃ eigenvalue), project and find dimension in null/live
print("JOINT C₂ × C₃ IRREP DECOMPOSITION")
print("="*60)
print(f"\n{'Irrep':<15} {'dim(null)':>10} {'dim(live)':>10} {'dim(total)':>10}")
print("-"*50)

total_null = 0
total_live = 0

for c2_label, c2_val in [('+1', 1.0), ('-1', -1.0)]:
    for c3_label, c3_val in [('1', 1.0), ('ω', omega), ('ω²', omega**2)]:
        # Projector onto joint eigenspace
        # P_c3 = (1/3)(I + ω^{-k}P + ω^{-2k}P²)
        # P_c2 = (1/2)(I + s·Q)  where s = ±1
        
        P2 = P @ P
        proj_c3 = (np.eye(bs) + np.conj(c3_val)*P + np.conj(c3_val)**2 * P2) / 3.0
        proj_c2 = (np.eye(bs) + c2_val * Q) / 2.0
        proj = np.real(proj_c3 @ proj_c2)
        
        # Dimension = trace of projector
        dim_total = round(np.real(np.trace(proj)))
        
        # Project onto null space
        proj_null = null_vecs.T @ proj @ null_vecs
        dim_null = round(np.real(np.trace(proj_null)))
        
        proj_live = live_vecs.T @ proj @ live_vecs
        dim_live = round(np.real(np.trace(proj_live)))
        
        label = f"({c2_label}, {c3_label})"
        print(f"{label:<15} {dim_null:>10} {dim_live:>10} {dim_total:>10}")
        total_null += dim_null
        total_live += dim_live

print("-"*50)
print(f"{'Total':<15} {total_null:>10} {total_live:>10} {total_null+total_live:>10}")

# Now compute the eigenvalue content of each irrep in the live space
print(f"\n\nLIVE SPACE: EIGENVALUE × IRREP DECOMPOSITION")
print("="*60)
# Group live eigenvalues
tol = 1e-6
live_evals = evals[np.abs(evals) >= 1e-8]
unique_evals = []
i = 0
while i < len(live_evals):
    val = live_evals[i]
    j = i+1
    while j < len(live_evals) and abs(live_evals[j]-val) < tol:
        j += 1
    unique_evals.append((val, i, j))  # value, start, end in live_vecs
    i = j

print(f"\n{'λ':>12} {'mult':>5} ", end='')
for c2l in ['+1', '-1']:
    for c3l in ['1', 'ω', 'ω²']:
        print(f" ({c2l},{c3l})", end='')
print()
print("-"*80)

for val, si, ei in unique_evals:
    sub = live_vecs[:, si:ei]
    mult = ei - si
    print(f"{val:>12.4f} {mult:>5} ", end='')
    
    for c2_val in [1.0, -1.0]:
        for c3_val in [1.0, omega, omega**2]:
            P2 = P @ P
            proj_c3 = (np.eye(bs) + np.conj(c3_val)*P + np.conj(c3_val)**2 * P2) / 3.0
            proj_c2 = (np.eye(bs) + c2_val * Q) / 2.0
            proj = np.real(proj_c3 @ proj_c2)
            proj_sub = sub.T @ proj @ sub
            dim = round(np.real(np.trace(proj_sub)))
            print(f"    {dim:>4}", end='')
    print()
