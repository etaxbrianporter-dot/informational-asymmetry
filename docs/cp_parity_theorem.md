# CP Parity Selection Theorem

## Statement

**Theorem (CP Parity Selection Rule).** Let O be the vacuum overlap matrix of K₂ₙ (n ≥ 4), restricted to the vacuum orbit of 2kp matchings (k = torus orbit length, p = 2n−1). Let the symmetry group G act on the vacuum sector, where G = C₂ₖ when 2^k ≡ −1 mod p, and G = C₂ × Cₖ otherwise. Then:

1. O has rank kp (half of 2kp) — the Half-Rank Law.
2. Each irrep of G has dimension p (odd).
3. In every CP-even irrep: null = (p−1)/2, live = (p+1)/2.
4. In every CP-odd irrep: null = (p+1)/2, live = (p−1)/2.

The CP-odd null excess per irrep is exactly 1. Total CP-odd excess = k.

## Proof

### Step 1: Conjugation preserves direction blocks

The conjugation map Q: v → (p−v) mod p sends direction class d = min(δ, p−δ) to itself, since min(−δ mod p, p−(−δ mod p)) = min(p−δ, δ) = d. Therefore Q permutes matchings *within* each direction block.

**Verified computationally:** At K₈, K₁₀, K₁₂, the overlap ⟨**1**_block_i | Q | **1**_block_j⟩ = δ_{ij} · bs_block. No inter-block mixing.

### Step 2: Block-Fourier vectors are CP-even

Define k vectors:

  v_j = Σ_{i=0}^{k−1} ω^{ij} · **1**_{block_i},  j = 0, ..., k−1

where ω = exp(2πi/k). Since Q preserves each block and Q · **1**_{block_i} = **1**_{block_i}, we have Q · v_j = v_j for all j. All block-Fourier vectors are CP-even.

### Step 3: Block-Fourier vectors are live eigenvectors of O

The overlap matrix acts on block indicator vectors through the inter-block overlap structure. Since all direction blocks in the vacuum orbit have identical internal structure (they are related by the torus automorphism), the matrix ⟨**1**_{block_i} | O | **1**_{block_j}⟩ depends only on the relative position of blocks i and j in the orbit. The Fourier modes diagonalize this circulant structure.

**Computed eigenvalues:**

| Level | v₀ | v₁ | v₂ | v₃ | v₄ |
|-------|-----|-----|-----|-----|-----|
| K₈  | 48 | 12 | 12 | — | — |
| K₁₀ | 76 | 28 | 28 | — | — |
| K₁₂ | 120 | 42.1 | 77.9 | 77.9 | 42.1 |

All nonzero. These are precisely the non-generic eigenvalues of the full spectrum (those with multiplicity ≠ 2k).

### Step 4: One block-Fourier vector per CP-even irrep

**Key computation.** The irrep decomposition of each v_j was measured by projecting onto the eigenbasis of the symmetry generator P:

| Level | Vector | Eigenvalue | Irrep | CP |
|-------|--------|-----------|-------|-----|
| K₁₂ | v₀ | 120.0 | j=0 | + |
| K₁₂ | v₁ | 42.1 | j=8 | + |
| K₁₂ | v₂ | 77.9 | j=6 | + |
| K₁₂ | v₃ | 77.9 | j=4 | + |
| K₁₂ | v₄ | 42.1 | j=2 | + |

Each v_j sits in a **single** irrep with weight 1.0000. The k vectors distribute one-to-one across the k CP-even irreps (j = 0, 2, 4, 6, 8 in C₁₀). Same pattern verified at K₈ and K₁₀.

**Why one-to-one:** In C₂ₖ, the P-eigenvalue of v_j is determined by how P permutes direction blocks. Since P (the ×2 map) generates the orbit of length k, the P-eigenvalue of v_j is a k-th root of unity raised to a power that depends on the block permutation. The map j → irrep label is a bijection onto the even-index irreps.

### Step 5: The ±1 split follows

Within each p-dimensional irrep of G:
- **CP-even:** The block-Fourier vector provides one guaranteed live mode (nonzero O-eigenvalue). The remaining p−1 dimensions partition into null and live.
- **CP-odd:** No such guaranteed mode exists. All p dimensions compete for null/live status.

The total null count is kp (Half-Rank Law). With 2k irreps each of dimension p, and one guaranteed live mode per CP-even irrep:

  null(CP+) + null(CP−) = kp
  null(CP+) ≤ k(p−1)  [since one mode per CP+ irrep is live]
  null(CP−) ≤ kp       [no constraint]

The unique solution with uniform per-irrep distribution:

  **null per CP+ irrep = (p−1)/2**
  **null per CP− irrep = (p+1)/2**

Check: k·(p−1)/2 + k·(p+1)/2 = kp ✓

### Step 6: Uniformity across Cₖ sectors

The per-irrep spectrum confirms the split is exactly (p±1)/2, not some other partition. The eigenvalues of O restricted to each irrep were computed explicitly:

**K₁₂, C₁₀ irreps (each 11-dimensional):**

| Irrep | CP | Live eigenvalues | null | live |
|-------|-----|-----------------|------|------|
| j=0 | + | 120, 44.4, 21.4, 15.3, 11.2, 3.7 | 5 | 6 |
| j=1 | − | 44.4, 21.4, 15.3, 11.2, 3.7 | 6 | 5 |
| j=2 | + | 44.4, 42.1, 21.4, 15.3, 11.2, 3.7 | 5 | 6 |
| j=3 | − | 44.4, 21.4, 15.3, 11.2, 3.7 | 6 | 5 |
| j=4 | + | 77.9, 44.4, 21.4, 15.3, 11.2, 3.7 | 5 | 6 |
| j=5 | − | 44.4, 21.4, 15.3, 11.2, 3.7 | 6 | 5 |
| j=6 | + | 77.9, 44.4, 21.4, 15.3, 11.2, 3.7 | 5 | 6 |
| j=7 | − | 44.4, 21.4, 15.3, 11.2, 3.7 | 6 | 5 |
| j=8 | + | 44.4, 42.1, 21.4, 15.3, 11.2, 3.7 | 5 | 6 |
| j=9 | − | 44.4, 21.4, 15.3, 11.2, 3.7 | 6 | 5 |

**Structure:** All CP-odd irreps share the identical "common core" of (p−1)/2 = 5 eigenvalues: {44.4, 21.4, 15.3, 11.2, 3.7}. Each CP-even irrep has this same core plus one additional eigenvalue from the block-Fourier set {120, 42.1, 77.9, 77.9, 42.1}.

## Physical Interpretation

The ±1 split is a **combinatorial Goldstone mechanism**:

- The k block-Fourier vectors are the matching-space analogues of "would-be Goldstone bosons" — collective modes of the vacuum orbit that couple to the symmetry breaking direction (CP-even, block-coherent).
- These modes are absorbed into the live sector of CP-even irreps, leaving CP-odd with one fewer live mode.
- The asymmetry is exactly 1 per irrep because there is exactly one block-Fourier vector per CP-even irrep — the minimum needed to break the null/live symmetry when p is odd.

## Verified At

| Level | p | k | Group | null(CP+) per irrep | null(CP−) per irrep | Status |
|-------|---|---|-------|-------|-------|--------|
| K₈  | 7 | 3 | C₂ × C₃ | 3 = (p−1)/2 | 4 = (p+1)/2 | ✓ |
| K₁₀ | 9 | 3 | C₆ | 4 = (p−1)/2 | 5 = (p+1)/2 | ✓ |
| K₁₂ | 11 | 5 | C₁₀ | 5 = (p−1)/2 | 6 = (p+1)/2 | ✓ |

## Open

The uniformity across Cₖ sectors (same rank in every CP+ irrep, same in every CP−) is verified computationally but the algebraic proof of uniformity requires showing that O restricted to the complement of the block-Fourier subspace has equal rank in every Cₖ sector. The "common core" eigenvalue set being identical across all irreps (at K₁₀ and K₁₂) suggests this follows from the block-circulant structure of O modulo the block-Fourier subspace, but the full argument remains to be formalized.
