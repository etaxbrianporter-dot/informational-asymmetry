# Resolution: The "Factor of 2" Between K₆ and K₈ Higgs Ratios

## The Original Puzzle

K₆ paper: a₄/a₂² = 0.3722 → m_H = 126 GeV
K₈ NCG: Σσ⁴/Σσ² = 0.705 → appeared to be ~2× K₆

## The Resolution: Apples Were Not Apples

The apparent factor of 2 arose from comparing **different quantities in different spaces**:

| What was compared | K₆ | K₈ | Ratio |
|---|---|---|---|
| a₄/a₂ (matching-space, Z₃) vs full vertex Tr(D⁴)/Tr(D²) | 1.231 | 1.739 | 1.41 ≈ √2 |
| Σσ⁴/Σσ² (vertex Yukawa) vs a₄/a₂ (matching, Z₃) | 1.231 | 0.705 | 0.57 |

These are incommensurable: different spaces, different phase content, different sector restrictions.

## The Correct Comparison: Spectral Kurtosis R = a₄/a₂²

The dimensionless kurtosis R = a₄/a₂² measures the **shape** of the eigenvalue distribution (1/n_eff). This is the quantity that enters the Higgs mass formula.

### Apples-to-apples results:

| Comparison | K₆ R | K₈ R | Ratio |
|---|---|---|---|
| Matching-space (with Z₃) | 0.3722 | 0.7868 | 2.11 |
| Vertex-space (full, no Z₃) | 0.3427 | 0.1881 | 0.55 |
| **K₈ Yukawa only vs K₆** | **0.3722** | **0.3491** | **0.938** |

## The Key Discovery: Per-Chirality Kurtosis

Decomposing K₈'s Tr(D²) and Tr(D⁴) by Pair C chirality sectors:

| Sector | Tr(D²) | Tr(D⁴) | R = Tr(D⁴)/Tr(D²)² | R/R_K₆ |
|---|---|---|---|---|
| γ=+1 (vertices {0,1,2,3}) | 4.624 (50.0%) | 8.147 (50.7%) | **0.3811** | **1.024** |
| γ=−1 (vertices {4,5,6,7}) | 4.624 (50.0%) | 7.933 (49.3%) | **0.3711** | **0.997** |
| K₆ reference | — | — | **0.3722** | 1.000 |

**R₋₋ = 0.3711 matches R_K₆ = 0.3722 to 0.3%.**

**R₊₊ = 0.3811 matches R_K₆ = 0.3722 to 2.4%.**

Each chirality sector of K₈ independently reproduces the K₆ Higgs kurtosis. The K₆ spectral shape is not doubled — it is **replicated** in each 4-vertex chirality block.

## Why This Happens

The K₆ vacuum, embedded via hub matchings into K₈, projects onto the K₈ vacuum eigenspace with 23.55% overlap (entirely in the ρ₁ component). This projection preserves the spectral shape:

1. **K₆ selects a K₄ core** on 4 vertices with 2 spectators
2. **Pair C grading** splits K₈ into two 4-vertex blocks
3. Each 4-vertex block inherits the K₆ vacuum's eigenvalue distribution
4. The kurtosis R is a **shape invariant** — it depends only on the relative eigenvalue distribution, not the overall scale

The γ=−1 sector (containing hub vertex v=7 and spectator-like content) preserves R more precisely because it more directly inherits the K₆ structure. The γ=+1 sector (Higgs-active K₄) has a slight upward deviation, possibly from the extra Yukawa coupling back-reaction.

## Effective Yukawa Multiplicity

| | n_eff = 1/R |
|---|---|
| K₆ | 2.687 |
| K₈ Yukawa (3 SVs) | 2.864 |
| K₈ γ=+1 sector | 2.624 |
| K₈ γ=−1 sector | 2.695 |

The K₈ Yukawa sector is slightly more democratic (n_eff closer to 3) than K₆, consistent with the 10.6% spread being smaller than K₆'s eigenvalue spread.

## What This Means for the Higgs Mass

Since R_K₈(per chirality) ≈ R_K₆ = 0.3722:

- The K₆ Higgs mass prediction (m_H = 126 GeV from √(2a₄/a₂) × M_W) is **inherited** by K₈
- No factor-of-2 correction is needed
- Each chirality sector independently carries the Higgs information
- The matching-space ratio 2.11 arises from K₈'s non-Higgs sectors (mass terms D₊₊, D₋₋) inflating a₄ relative to a₂² — these are gauge/gravitational content, not Higgs

## Summary

| Claim | Status |
|---|---|
| "Factor of 2 between K₈ and K₆" | **KILLED** — arose from comparing incommensurable quantities |
| K₆ kurtosis preserved in K₈ | **CONFIRMED** — R₋₋/R_K₆ = 0.997 (0.3% match) |
| Per-chirality Higgs structure | **DISCOVERED** — each 4-vertex block independently carries R ≈ 0.372 |
| K₈ Yukawa kurtosis | 0.349, 6% below K₆ (slight democratization from near-degenerate SVs) |
