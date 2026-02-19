# Overlap Moment Results — Banked Structure
## February 19, 2026

### Verified Exact Formulas

**Overlap matrix eigenvalues** for K₂ₙ (3 distinct, corresponding to S₂ₙ irreps):

| Level | λ_max (trivial, ×1) | λ_mid (physical, ×dim) | λ_min (kernel) |
|-------|---------------------|----------------------|----------------|
| K₄   | 4 (×3, all degenerate) | — | — |
| K₆   | 18 (×1) | 8 (×9) | 0 (×5) |
| K₈   | 120 (×1) | 36 (×20) | 0 (×84) |

General: λ_max = 2n(2n-3)!!, λ_mid = 4(n-1)(2n-5)!!, kernel dim = (2n-1)!! - dim(V_{(2n-2,2)}) - 1

**K₄ is special:** O₄ = 4I₃. The 3 matchings partition K₄'s 6 edges into 3 disjoint pairs. No overlap at all. This is unique to K₄.

**Vertex-resolved block structure:**
- At each vertex v, matchings partition into (2n-1) blocks of size (2n-3)!!
- Each O_v is block-diagonal with all-ones blocks
- Exact: Tr(O_v^k) = (2n-1) × ((2n-3)!!)^k
- Summed: Σ_v Tr(O_v^k) = 2n(2n-1)((2n-3)!!)^k

**Vertex-space average adjacency:** D = (1/(2n-1))(J-I) for all K₂ₙ. Universal, no interesting variation. The matching algebra's projection to vertex space is the complete graph Laplacian.

### What These Are NOT

These are **intra-level** spectral geometry — structure within a single K₂ₙ matching algebra. They are NOT the Seeley-DeWitt f₀ correction to gauge couplings, which involves D_F² = Y†Y connecting BETWEEN levels.

### f₀ Investigation Status: CLOSED (mechanism identified, not yet quantitative)

**Finding:** The f₀ boundary correction uses the same aᵢ = {17/10, 3/2, 2} coefficients as 2-loop Yukawa running. The correction has the RIGHT SIGN to reduce c₁/c₃ (because a₁/c₁ = 0.2125 < a₃/c₃ = 0.500).

**Problem 1:** a₂/c₂ - a₃/c₃ = -0.25 ≠ 0, so f₀ also disturbs c₂/c₃ = 3/2. Not a clean single-parameter fix.

**Problem 2:** Required ρ = f₀/(f₂Λ²) ≈ 0.41. This is a 40% correction — the Seeley-DeWitt expansion is not converging at Λ ~ 10⁸ GeV. Either the cutoff is too low for asymptotic expansion, or the full self-consistent system (boundary + running + thresholds) is needed.

**Honest status of c₁ = 8:** 2.5% from integer at 2-loop. The mechanism (spectator counting) is confirmed. The residual points to the boundary/running interplay, which requires a self-consistent treatment beyond perturbative Seeley-DeWitt. The direction of the correction is right; the magnitude needs the full calculation.

### Where Overlap Moments May Matter

- Higgs potential: a₂, a₄ from K₆ involve the physical subspace (λ_mid eigenspace)
- Spectral dimension of matching algebras as noncommutative spaces  
- Heat kernel / zeta function of matching geometry
- Continuum limit: how discrete matching structure approaches continuous manifold
- The 84-dimensional kernel of K₈'s overlap → gauge redundancy in fermion sector
