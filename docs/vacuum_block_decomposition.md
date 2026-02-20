# Vacuum Block Decomposition: Half-Rank is Global

## The Question

The vacuum orbit overlap matrix O_vac has the half-rank law (Law 1): exactly half its eigenvalues are zero. Each vacuum orbit decomposes into k direction blocks V_0, ..., V_{k-1} under the torus action ×2 mod p.

**Does each block already have half-rank individually, or does the half-rank only emerge from interference between blocks?**

## Answer: Half-Rank is GLOBAL

Individual direction blocks are nearly full rank. Only 1 null mode per block — nowhere near bs/2.

| Level | Block size | Block rank | Block null | bs/2 |
|-------|-----------|-----------|-----------|------|
| K₈    | 14        | 13        | 1         | 7    |
| K₁₀   | 18        | 17        | 1         | 9    |
| K₁₂   | 22        | 21        | 1         | 11   |

**Pattern: block null = 1, universally.** Each direction sector is (2n+2)-dimensional with exactly one null mode (the "constant mode" with null projector diagonal = 1/(2n-1)). The remaining half-rank deficit is created entirely by destructive interference when blocks combine.

## Block-Circulant Structure

The vacuum overlap matrix is a k×k block-circulant matrix under the correct generator g = 2^(ord(2 mod p)/k) mod p. This generator has order exactly k on vertices, resolving the stabilizer twist from ×(-1).

| Level | p  | ord(2) | k | Generator | Block sequence |
|-------|----|--------|---|-----------|----------------|
| K₈    | 7  | 3      | 3 | g=2       | [0, 1, 2]      |
| K₁₀   | 9  | 6      | 3 | g=4       | [0, 2, 1]      |
| K₁₂   | 11 | 10     | 5 | g=4       | [0, 2, 4, 1, 3] |

**Why g=4, not g=2?** The direction permutation ×2 has order 2k on vertices but only k on matchings (each matching is stabilized by ×(-1) = ×2^k). The naive circulant with g=2 fails because π^{-1}_vertices ≠ π^{k-1}_vertices. Using g = 2^(ord/k) gives a vertex map of order exactly k, making the block-circulant exact.

### Reconstruction: Exact at All Levels

```
K₈:  reconstruction error = 2.84e-14  ✓
K₁₀: reconstruction error = 2.84e-14  ✓
K₁₂: reconstruction error = 7.11e-14  ✓
```

Fourier decomposition: O_vac → {M_ℓ = Σ_j B_j ω^{jℓ}}_ℓ=0,...,k-1 where ω = e^{2πi/k}.

## The Democratic Half-Rank Theorem

**Every Fourier channel independently has exactly bs_blk/2 null modes.**

| Level | bs_blk | Null/channel | Channels | Total null |
|-------|--------|-------------|----------|-----------|
| K₈    | 14     | 7           | 3        | 21        |
| K₁₀   | 18     | 9           | 3        | 27        |
| K₁₂   | 22     | 11          | 5        | 55        |

**Formula:** null_total = k × (bs_blk/2) = k × (n+1) = (2n+2)k/2 = bs_total/2. ✓

This is not a trivial consequence of the block structure. Individual blocks have rank bs-1 (only 1 null mode), but the Fourier transform of the off-diagonal coupling creates exactly (bs/2 - 1) additional null modes in every single channel. The null space is distributed with perfect democracy across the representation theory of ℤ_k.

## Channel Spectra: Universal Bulk with Distinguished Top

The non-zero eigenvalues show remarkable structure:

### K₈ (k=3 channels, bs_blk=14)
```
ℓ=0: 48.00, {24.49}×2, {12.00}×2, {7.51}×2, {0}×7
ℓ=1: {24.49}×2, {12.00}×3, {7.51}×2, {0}×7
ℓ=2: {24.49}×2, {12.00}×3, {7.51}×2, {0}×7
```

### K₁₀ (k=3 channels, bs_blk=18)
```
ℓ=0: 76.00, {29.00}×2, {16.00}×2, {11.80}×2, {7.20}×2, {0}×9
ℓ=1: {29.00}×2, 28.00, {22.00}×2, {11.80}×2, {7.20}×2, {0}×9
ℓ=2: {29.00}×2, 28.00, {22.00}×2, {11.80}×2, {7.20}×2, {0}×9
```

### K₁₂ (k=5 channels, bs_blk=22)
```
ℓ=0: 120.00, {44.36}×2, {21.41}×2, {15.32}×2, {11.21}×2, {3.71}×2, {0}×11
ℓ=1: 77.89, {44.36}×2, {21.41}×2, {15.32}×2, {11.21}×2, {3.71}×2, {0}×11
ℓ=2: 42.11, {44.36}×2, {21.41}×2, {15.32}×2, {11.21}×2, {3.71}×2, {0}×11
ℓ=3: 42.11, {44.36}×2, {21.41}×2, {15.32}×2, {11.21}×2, {3.71}×2, {0}×11
ℓ=4: 77.89, {44.36}×2, {21.41}×2, {15.32}×2, {11.21}×2, {3.71}×2, {0}×11
```

**Pattern:** Each channel has:
- **1 distinguished eigenvalue** (the "top mode") that varies by channel
- **n-1 universal doublet pairs** shared identically across ALL channels
- **bs/2 = n+1 null modes**

The distinguished top eigenvalues correspond exactly to the full orbit eigenvalues that have non-2k multiplicities (the "special" eigenvalues from Law 3).

At K₁₂: the five top eigenvalues are {120, 77.89, 42.11, 42.11, 77.89}, which are the orbit eigenvalues with multiplicities {1, 2, 2}. The 10 universal doublet eigenvalues {44.36, 21.41, 15.32, 11.21, 3.71}×2 fill out the multiplicity-10 eigenvalues of the full orbit.

**Interpretation:** The multiplicity-2k eigenvalues of the full orbit arise from the diagonal block spectrum — they're the same in every Fourier channel. The special eigenvalues (mult < 2k) arise from the off-diagonal coupling between direction sectors. The "tuning fork" analogy: the universal doublets are the prong frequencies (direction-sector physics), while the distinguished top mode is the coupling resonance that varies by channel.

## Vacuum Internal Coupling

Off-diagonal blocks (cross-direction coupling within the vacuum) also show the flat/structured dichotomy:

### K₈ (k=3)
```
shift=1: rank=13, hier=15.80, ||F||²=448  (structured)
shift=2: rank=13, hier=15.80, ||F||²=448  (structured, = shift 1 by symmetry)
```
All cross-direction couplings are structured. No internal flat couplings.

### K₁₀ (k=3)
```
shift=1: rank=17, hier=44.50, ||F||²=720  (structured)
shift=2: rank=17, hier=44.50, ||F||²=720  (structured)
```
Same pattern. All structured.

### K₁₂ (k=5)
```
shift=1: rank=11, hier=1.00, ||F||²=176   ** FLAT **
shift=2: rank=21, hier=36.11, ||F||²=1056  (structured)
shift=3: rank=21, hier=36.11, ||F||²=1056  (structured)
shift=4: rank=11, hier=1.00, ||F||²=176   ** FLAT **
```

K₁₂ has internal flat couplings! Shifts 1 and 4 (= k-1) are direction-disjoint within the vacuum. Shifts 2 and 3 carry rich hierarchical structure. This is the first level where the vacuum itself has informationally decoupled direction sectors.

### Frobenius Budget
| Level | Diagonal fraction | Off-diagonal fraction |
|-------|------------------|-----------------------|
| K₈    | 63.6%            | 36.4%                 |
| K₁₀   | 73.0%            | 27.0%                 |
| K₁₂   | 78.8%            | 21.2%                 |

Diagonal (same-direction) fraction grows with level. The vacuum becomes increasingly dominated by intra-direction self-coupling at higher levels.

## Physical Interpretation

### The Half-Rank Mechanism

1. **Each direction block is nearly full rank** (bs-1 out of bs). The one null mode per block is a "gauge" constraint.

2. **The block-circulant Fourier transform decomposes O_vac** into k independent channels M_ℓ, one per irrep of ℤ_k.

3. **Cross-direction coupling creates destructive interference** that independently kills exactly half the modes in every channel.

4. **The half-rank law is a spectral democracy**: it's not imposed by any single direction sector or any single Fourier channel. Every channel carries its own copy of the null space.

### Implications for Physics

The "tuning fork" picture refines: the vacuum resonance has k overtones (Fourier channels), each with its own frequency (distinguished eigenvalue) but the same quality factor (universal doublets determine the spectral shape). The null space — which controls what physics can emerge — is not localized in any particular overtone.

For the K₆ Higgs: the spectral ratio a₄/a₂ that gives 125 GeV involves eigenvalues from the full orbit. The block decomposition shows these eigenvalues split into universal bulk (from intra-direction physics) and distinguished modes (from inter-direction coupling). The physical prediction emerges from both contributions.

For the K₈ Yukawa hierarchy: the 415:135:1 ratio comes from the orbit-level coupling. Sense 3 showed this hierarchy is compressed by orbit averaging. The block-circulant Fourier decomposition shows the compression is democratic — every Fourier channel contributes equally to the null space that constrains the Yukawa couplings.

## Technical: The Stabilizer Twist

For K₂ₙ with n ≥ 4, the direction permutation π = ×2 mod p has order 2k on vertices (where k = orbit size). The stabilizer is π^k = ×(-1) mod p, which fixes all matchings but permutes vertices. This creates a ℤ₂ extension:

```
1 → ℤ₂ → ℤ_{2k} → ℤ_k → 1
```

The naive block-circulant using g=2 fails because it doesn't account for this extension. The correct generator g = 2^(ord/k) = 4 (for most levels) quotients out the stabilizer, giving a true circulant with exact spectral reconstruction.

The block visiting order [0, 2, 4, 1, 3] at K₁₂ is the Galois automorphism σ: ζ_5 → ζ_5^4 acting on the cyclotomic field ℚ(ζ_k). This connects the block-circulant decomposition to the Galois structure that governs discriminant factorizations at higher levels.

## Summary Table

| Property | K₈ | K₁₀ | K₁₂ |
|----------|-----|------|------|
| k (blocks) | 3 | 3 | 5 |
| bs_blk | 14 | 18 | 22 |
| Block rank | 13 | 17 | 21 |
| Block null | **1** | **1** | **1** |
| Null/channel | **7** | **9** | **11** |
| = bs_blk/2 | ✓ | ✓ | ✓ |
| Reconstruction | exact | exact | exact |
| Internal flat | 0 | 0 | 2/4 |
| Diagonal frac | 63.6% | 73.0% | 78.8% |
| Generator | ×2 | ×4 | ×4 |
| Block seq | 0,1,2 | 0,2,1 | 0,2,4,1,3 |

## Files

- `coupling_probe_unified.py`: Updated with `sense3_vacuum_blocks()` method
- Method includes: diagonal blocks, off-diagonal SVD, Frobenius budget, block-circulant Fourier decomposition with correct generator, null space by channel
