# Sense 3: Block-Level Coupling Results

## Question

Does the Yukawa hierarchy live at the individual direction-block level, or does it only emerge after orbit averaging (summing over the k direction blocks within each orbit)?

## Method

Decompose each orbit-level coupling matrix O(vacuum, other) into a k×k grid of block-block matrices O(Vᵢ, Bⱼ), where Vᵢ is vacuum direction block i and Bⱼ is the other orbit's direction block j. SVD each block pair independently.

The direction permutation (×2 mod p) cyclically permutes blocks within each orbit, so the k×k grid has circulant structure: block pair (Vᵢ, Bⱼ) is spectrally equivalent to (Vᵢ₊₁, Bⱼ₊₁). This reduces k² pairs to k distinct "shift classes."

## Key Findings

### 1. Hierarchy Exists at Block Level — Orbit Averaging COMPRESSES It

**The hierarchy is not an orbit-averaging artifact.** Individual same-direction (shift=0) block pairs often have LARGER σ_max/σ_min ratios than the orbit-level coupling.

| Level | Orbit | Orbit hier | Block hier (shift=0) | Ratio |
|-------|-------|-----------|----------------------|-------|
| K₈    | 0     | 10.79     | 5.70                 | 0.53× |
| K₁₀   | 7     | 5.50      | 167.70               | 30.5× |
| K₁₀   | 2     | 17.28     | 219.67               | 12.7× |
| K₁₀   | 5     | 8.60      | 141.33               | 16.4× |
| K₁₂   | *     | ~3-24     | up to ~50+           | varies |

At K₁₀ and K₁₂, block-level hierarchies regularly exceed orbit-level by 10-30×. The torus symmetry acts as a **regulator**: orbit averaging mixes the structured same-direction coupling with the flatter cross-direction couplings, compressing the range.

### 2. Flat Cross-Direction Couplings Are Universal

Certain shift classes produce perfectly flat SVD spectra — every singular value is identical (hierarchy = 1.000).

| Level | Flat block pairs | Total block pairs | Fraction |
|-------|-----------------|-------------------|----------|
| K₈    | 3               | 27                | 11%      |
| K₁₀   | 12              | 99                | 12%      |
| K₁₂   | 95              | 605               | 16%      |

Flat couplings occur when a non-vacuum orbit's direction multiset shares NO direction classes with the vacuum. In these cross-direction pairs, every vacuum matching overlaps every non-vacuum matching with exactly the same count — the coupling is democratic, carrying no information.

These are the block-level analog of the direction-disjoint sectors discovered at K₁₆.

### 3. Rank Cancellation Is Massive and Universal

Block ranks always sum to far more than the orbit rank:

| Level | Typical Σ block ranks | Orbit rank | Cancelled | Ratio |
|-------|-----------------------|------------|-----------|-------|
| K₈    | 63                    | 21         | 42        | 3:1   |
| K₁₀   | 117-153               | 27         | 90-126    | 4-6:1 |
| K₁₂   | 475-625               | 55         | 420-570   | 9-11:1|

The cancellation ratio grows with level. This is **destructive interference between direction blocks**: modes that are live in one block pair are killed when combined with other blocks. The orbit average selects only modes that survive interference across all k shift classes simultaneously.

### 4. The ||F||² Grid Is Circulant

For every non-FP orbit pair, the ||F||² grid has exact circulant structure:

```
K₈, Orbit 0:  [616, 56, 336]    → ratios 11:1:6
K₁₀, Orbit 0: [1944, 72, 360]   → ratios 27:1:5
K₁₀, Orbit 7: [1440, 144, 576]  → ratios 10:1:4
```

Same-direction (shift=0) always carries the largest ||F||², confirming that direction-matched blocks dominate the coupling. But Frobenius norm is exactly additive (Σ blocks = orbit total), so no energy is lost — just redistributed.

### 5. FP Blocks Are All Equivalent

The fixed-point orbit has a single block (all directions present). All k vacuum blocks couple identically to it:

| Level | FP block rank | FP orbit rank | Per-block cancelled |
|-------|--------------|---------------|---------------------|
| K₈    | 13           | 19            | (39-19)/3 ≈ 7      |
| K₁₀   | 17           | 25            | (51-25)/3 ≈ 9      |
| K₁₂   | 21           | 51            | (105-51)/5 ≈ 11    |

FP block rank follows the pattern: (block size)−(number of directions)+1 at K₈ and K₁₀.

### 6. Block Rank Depends on Direction Overlap

Same-direction blocks have LOWER rank than cross-direction blocks in some cases:

```
K₁₀, Orbit 0:  shift=0 → rank=7,  shift=1 → rank=9
K₁₀, Orbit 11: shift=0 → rank=13, shift=1 → rank=9 (flat)
```

This is counterintuitive: blocks that share direction classes have more constraints (lower rank) but richer structure (larger hierarchy). Blocks with no shared directions have full rank but carry no information (flat).

## Physical Interpretation

The block decomposition reveals the **mechanism of torus regularization**:

1. **Same-direction coupling carries the physics** — hierarchy, rank deficit, structure
2. **Cross-direction coupling carries noise** — flat, democratic, uninformative
3. **Orbit averaging = interference** — cancels ~80-90% of modes, selecting only those that survive across all direction classes
4. **The surviving modes are the physical ones** — they respect the full torus symmetry

This is analogous to how Fourier transforms project out specific momentum modes: the orbit average is a momentum-space projector, and the block coupling is position-space. The hierarchy in position space (blocks) is compressed by the projection into momentum space (orbits), but the essential structure — which modes couple to which — is preserved.

## Integration

Added to `coupling_probe_unified.py` as `sense3_block_coupling()`, controlled by `--no-blocks` flag.
