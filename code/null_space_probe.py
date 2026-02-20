#!/usr/bin/env python3
"""
Null space identification at K₈.
Q1: Do the 21 null vectors live within direction blocks or mix across them?
Q2: What wreath product irreps do they carry?
Q3: Is there a pattern connecting null vectors to specific matching structures?
"""
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

nv = 8
p = 7

raw = gen_matchings(nv)
match_sets = [set(m) for m in raw]
N = len(raw)

# Group by direction
blocks = defaultdict(list)
for i, m in enumerate(raw):
    blocks[torus_ms(m, p)].append(i)

# Vacuum orbit: blocks (0,0,2), (0,1,1), (1,2,2)
vac_blocks = [(0,0,2), (0,1,1), (1,2,2)]
vac_idx = []
block_ranges = {}  # which indices in vac_idx come from which block
for vb in vac_blocks:
    start = len(vac_idx)
    vac_idx.extend(blocks[vb])
    block_ranges[vb] = (start, len(vac_idx))

bs = len(vac_idx)
print(f"Vacuum orbit: {bs} matchings across 3 blocks")
for vb in vac_blocks:
    s, e = block_ranges[vb]
    print(f"  Block {vb}: indices [{s}:{e}], size {e-s}")

# Build overlap matrix
O = np.zeros((bs, bs))
for i in range(bs):
    for j in range(bs):
        if i == j:
            O[i,j] = nv
        else:
            O[i,j] = 2 * len(match_sets[vac_idx[i]] & match_sets[vac_idx[j]])

# Eigendecomposition
evals, evecs = np.linalg.eigh(O)
# Sort descending
order = np.argsort(evals)[::-1]
evals = evals[order]
evecs = evecs[:, order]

# Identify null space
tol = 1e-8
null_mask = np.abs(evals) < tol
null_vecs = evecs[:, null_mask]
live_vecs = evecs[:, ~null_mask]
n_null = null_vecs.shape[1]
n_live = live_vecs.shape[1]
print(f"\nNull space: dim {n_null}")
print(f"Live space: dim {n_live}")

# ======================================================================
# Q1: Block structure of null space
# ======================================================================
print(f"\n{'='*60}")
print("Q1: BLOCK SUPPORT OF NULL VECTORS")
print(f"{'='*60}")

# For each null vector, compute weight in each block
print("\nWeight distribution across blocks (sum of squared components):")
block_weights_null = np.zeros((n_null, 3))
for k in range(n_null):
    v = null_vecs[:, k]
    for bi, vb in enumerate(vac_blocks):
        s, e = block_ranges[vb]
        block_weights_null[k, bi] = np.sum(v[s:e]**2)

# Summary: total weight per block across all null vectors
total_null_weight = np.sum(block_weights_null, axis=0)
print(f"\n  Total null space weight per block:")
for bi, vb in enumerate(vac_blocks):
    print(f"    Block {vb}: {total_null_weight[bi]:.4f} / {n_null:.0f} = {total_null_weight[bi]/n_null:.4f}")

# Same for live space
block_weights_live = np.zeros((n_live, 3))
for k in range(n_live):
    v = live_vecs[:, k]
    for bi, vb in enumerate(vac_blocks):
        s, e = block_ranges[vb]
        block_weights_live[k, bi] = np.sum(v[s:e]**2)

total_live_weight = np.sum(block_weights_live, axis=0)
print(f"\n  Total live space weight per block:")
for bi, vb in enumerate(vac_blocks):
    print(f"    Block {vb}: {total_live_weight[bi]:.4f} / {n_live:.0f} = {total_live_weight[bi]/n_live:.4f}")

# Check if any null vector is block-pure (supported on single block)
print(f"\n  Block purity of each null vector (max block weight):")
max_weights = np.max(block_weights_null, axis=1)
pure_count = np.sum(max_weights > 0.999)
mixed_count = np.sum(max_weights < 0.501)
print(f"    Pure (>99.9% in one block): {pure_count}")
print(f"    Maximally mixed (<50.1% in any block): {mixed_count}")
print(f"    In between: {n_null - pure_count - mixed_count}")

# ======================================================================
# Q2: C₃ action on null space  
# ======================================================================
print(f"\n{'='*60}")
print("Q2: C₃ (DIRECTION PERMUTATION) ACTION ON NULL SPACE")
print(f"{'='*60}")

# The C₃ action permutes: (0,0,2)→(0,1,1)→(1,2,2)→(0,0,2)
# This means it permutes the 3 blocks. We need the matching-level permutation.
# Direction perm: 0→1, 1→2, 2→0
dir_perm = {0: 1, 1: 2, 2: 0}

# Build the vertex permutation that induces this direction perm
# ×2 mod 7: vertex v → 2v mod 7 (vertex 7 is fixed as "extra")
def apply_vertex_perm(matching, p):
    """Apply ×2 mod p to matching edges."""
    new_edges = set()
    for a, b in matching:
        na = (2 * a) % p if a < p else a
        nb = (2 * b) % p if b < p else b
        new_edges.add((min(na, nb), max(na, nb)))
    return frozenset(new_edges)

# Build permutation matrix on vacuum orbit
# First, create lookup: matching → index in vac_idx
match_to_vacidx = {}
for vi, gi in enumerate(vac_idx):
    match_to_vacidx[raw[gi]] = vi

# Permutation matrix P (42×42)
P = np.zeros((bs, bs))
for vi, gi in enumerate(vac_idx):
    m = raw[gi]
    m_perm = apply_vertex_perm(m, p)
    if m_perm in match_to_vacidx:
        vj = match_to_vacidx[m_perm]
        P[vj, vi] = 1.0
    else:
        print(f"  WARNING: permuted matching not in vacuum orbit!")

# Verify P is a permutation matrix
assert np.allclose(P @ P.T, np.eye(bs)), "P is not orthogonal!"
# Verify P³ = I (C₃ action)
P3 = P @ P @ P
print(f"  P³ = I: {np.allclose(P3, np.eye(bs))}")
# Verify P commutes with O
print(f"  POP⁻¹ = O: {np.allclose(P @ O @ P.T, O)}")

# Project null space onto C₃ eigenspaces
# C₃ eigenvalues: 1, ω, ω² where ω = e^{2πi/3}
omega = np.exp(2j * np.pi / 3)

# Restrict P to null space
P_null = null_vecs.T @ P @ null_vecs  # 21×21 matrix

# Eigenvalues of P restricted to null space
p_evals = np.linalg.eigvals(P_null)
print(f"\n  Eigenvalues of C₃ on null space:")
# Count by C₃ eigenvalue
n_trivial = np.sum(np.abs(p_evals - 1.0) < 0.01)
n_omega = np.sum(np.abs(p_evals - omega) < 0.01)
n_omega2 = np.sum(np.abs(p_evals - omega**2) < 0.01)
print(f"    ω⁰ = 1 (trivial):  {n_trivial}")
print(f"    ω¹ (rotation):      {n_omega}")
print(f"    ω² (rotation):      {n_omega2}")
print(f"    Total: {n_trivial + n_omega + n_omega2}")

# Same analysis for live space
P_live = live_vecs.T @ P @ live_vecs
p_evals_live = np.linalg.eigvals(P_live)
n_triv_l = np.sum(np.abs(p_evals_live - 1.0) < 0.01)
n_om_l = np.sum(np.abs(p_evals_live - omega) < 0.01)
n_om2_l = np.sum(np.abs(p_evals_live - omega**2) < 0.01)
print(f"\n  Eigenvalues of C₃ on live space:")
print(f"    ω⁰ = 1 (trivial):  {n_triv_l}")
print(f"    ω¹ (rotation):      {n_om_l}")
print(f"    ω² (rotation):      {n_om2_l}")

print(f"\n  C₃ character table:")
print(f"    {'Space':<15} {'dim':>4} {'χ(1)':>6} {'χ(ω)':>6} {'χ(ω²)':>6}")
print(f"    {'Null':<15} {n_null:>4} {n_trivial:>6} {n_omega:>6} {n_omega2:>6}")
print(f"    {'Live':<15} {n_live:>4} {n_triv_l:>6} {n_om_l:>6} {n_om2_l:>6}")
print(f"    {'Total':<15} {bs:>4} {n_trivial+n_triv_l:>6} {n_omega+n_om_l:>6} {n_omega2+n_om2_l:>6}")

# ======================================================================
# Q2b: C₂ action (edge swap within direction class)
# ======================================================================
print(f"\n{'='*60}")
print("Q2b: FULL WREATH PRODUCT C₂ ≀ C₃ ANALYSIS")
print(f"{'='*60}")

# C₂ acts as ×(p-1) mod p = ×6 mod 7: vertex v → 6v mod 7 (= -v mod 7)
# This is the "conjugation" that swaps d and p-d (but since we already folded, 
# it acts as identity on direction classes). Let's check.
def apply_conj(matching, p):
    """Apply ×(-1) mod p to matching."""
    new_edges = set()
    for a, b in matching:
        na = ((p-1) * a) % p if a < p else a
        nb = ((p-1) * b) % p if b < p else b
        new_edges.add((min(na, nb), max(na, nb)))
    return frozenset(new_edges)

Q = np.zeros((bs, bs))
for vi, gi in enumerate(vac_idx):
    m = raw[gi]
    m_conj = apply_conj(m, p)
    if m_conj in match_to_vacidx:
        vj = match_to_vacidx[m_conj]
        Q[vj, vi] = 1.0

print(f"  Q² = I: {np.allclose(Q @ Q, np.eye(bs))}")
print(f"  QOQ⁻¹ = O: {np.allclose(Q @ O @ Q.T, O)}")
print(f"  PQ = QP? {np.allclose(P @ Q, Q @ P)}")

# Q on null/live spaces
Q_null = null_vecs.T @ Q @ null_vecs
q_evals_null = np.linalg.eigvals(Q_null)
n_plus = np.sum(np.abs(q_evals_null - 1.0) < 0.01)
n_minus = np.sum(np.abs(q_evals_null + 1.0) < 0.01)
print(f"\n  C₂ (conjugation) on null space:")
print(f"    +1 eigenvalue: {n_plus}")
print(f"    -1 eigenvalue: {n_minus}")

Q_live = live_vecs.T @ Q @ live_vecs
q_evals_live = np.linalg.eigvals(Q_live)
n_plus_l = np.sum(np.abs(q_evals_live - 1.0) < 0.01)
n_minus_l = np.sum(np.abs(q_evals_live + 1.0) < 0.01)
print(f"\n  C₂ (conjugation) on live space:")
print(f"    +1 eigenvalue: {n_plus_l}")
print(f"    -1 eigenvalue: {n_minus_l}")

# ======================================================================
# Q3: Matching structure of null vectors
# ======================================================================
print(f"\n{'='*60}")
print("Q3: WHAT MAKES A MATCHING 'SILENT'?")
print(f"{'='*60}")

# For each matching in the vacuum orbit, compute its total weight in null space
null_projector = null_vecs @ null_vecs.T  # 42×42 projector onto null space
null_weight = np.diag(null_projector)

print(f"\n  Null space weight per matching (diagonal of projector):")
print(f"  (1.0 = fully in null space, 0.0 = fully live)")

for bi, vb in enumerate(vac_blocks):
    s, e = block_ranges[vb]
    print(f"\n  Block {vb}:")
    for vi in range(s, e):
        gi = vac_idx[vi]
        m = sorted(raw[gi])
        w = null_weight[vi]
        marker = " ◄ NULL" if w > 0.9 else (" ◄ LIVE" if w < 0.1 else "")
        print(f"    {m}  w_null = {w:.4f}{marker}")

# Summary statistics
print(f"\n  Distribution of null weights:")
pure_null = np.sum(null_weight > 0.99)
pure_live = np.sum(null_weight < 0.01)
print(f"    Fully null (>0.99): {pure_null}")
print(f"    Fully live (<0.01): {pure_live}")
print(f"    Mixed: {bs - pure_null - pure_live}")
print(f"    Mean null weight: {np.mean(null_weight):.4f} (expected: {n_null/bs:.4f})")

# Block-level summary
print(f"\n  Mean null weight by block:")
for bi, vb in enumerate(vac_blocks):
    s, e = block_ranges[vb]
    print(f"    Block {vb}: {np.mean(null_weight[s:e]):.4f}")
