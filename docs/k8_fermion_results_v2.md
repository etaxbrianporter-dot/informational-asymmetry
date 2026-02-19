---
title: K₈ Fermion Sector: Structural Results (v2)
section: K₈ Fermion Sector
status: active
---

# K₈ Fermion Sector: Structural Results (v2)

## Framework

K₈ (complete graph on 8 vertices, 28 edges, 105 perfect matchings) embeds on
genus 2 (double torus). Edge directions from K₇ Heawood classification:
- d₀ = ±1 mod 7 (7 edges, torus direction 1)
- d₁ = ±2 mod 7 (7 edges, torus direction 2)
- d₂ = ±3 mod 7 (7 edges, torus direction 3)
- d₃ = handle edges (i,7) (7 edges, second handle)

Z₃ phases on torus directions; handle phase is gauge-redundant (confirmed).

Generation pairing: (0,1) = Higgs doublet, (2,3)/(4,5)/(6,7) = 3 generations.

## Vacuum Structure

The 105×105 Gram matrix (direction-weighted, Z₃-phased) has:
- 5 null modes (momentum conservation)
- Physical vacuum: λ = 1.9595, **6-fold degenerate**
- 6D = 2D(ρ₁) + 2D(ρ₂) + 2D(ρ₃) under Z₇ symmetry
- Zero trivial component → vacuum MUST break Z₇

## Main Result: K₆ Vacuum → K₈ Yukawa (Zero Parameters)

The physical K₆ vacuum (the one giving mH = 125 GeV) embeds as hub matchings
in K₈. Its projection onto the 6D vacuum eigenspace is:

**Purely ρ₁ (100%)** — single Z₇ irrep, zero free parameters.

### Zero-Parameter Yukawa Hierarchy

**415 : 135 : 1** (rank 3, three massive generations)

## Selection Rules (Proven)

| Z₇ irrep | Accessible from K₆? | Rank | Hierarchy |
|-----------|---------------------|------|-----------|
| ρ₁ (k=1,6) | ✓ (physical vacuum) | 3 | 415:135:1 |
| ρ₂ (k=2,5) | ✗ (never, from ANY K₆ eigenvector) | — | decoupled |
| ρ₃ (k=3,4) | ✓ (excited K₆ states only) | **2** | ~77:1 |

**ρ₂ is structurally inaccessible from the Higgs sector.**

### Physical Interpretation

- ρ₁: coupled to Higgs → charged fermion masses (3 generations, hierarchical)
- ρ₃: coupled to excited Higgs states → rank 2 (one massless generation)
- ρ₂: **decoupled from Higgs** → candidate neutrino sector (mass from non-Higgs mechanism)

## Z₇ Doublet Shifts

The 7 doublet positions within K₇ give different hierarchies (all purely ρ₁):

| Doublet | Hierarchy | Rank |
|---------|-----------|------|
| {0,1} | 415:135:1 | 3 |
| {1,2} | 544:74:1 | 3 |
| {2,3} | 25:2:1 | 3 |
| {3,4} | 7:1 | 2 |
| {4,5} | 420:109:1 | 3 |
| {5,6} | 18:2:1 | 3 |
| {6,0} | 54:1 | 2 |

Doublets far from vertex 7 (handle vertex) → rank 3, large hierarchy.
Doublets adjacent to vertex 7 → rank 2, small hierarchy.

## Key Algebraic Facts

1. dim(vacuum) = 6 = dim(so(4)) = dim(SU(2)_L × SU(2)_R) (custodial symmetry)
2. 6 = 2+2+2 under Z₇ (three conjugate irrep pairs)
3. K₆ hub → ρ₁ ⊕ ρ₃ only (ρ₂ selection rule)
4. Handle phase is gauge-redundant (every matching has exactly 1 handle edge)
5. 105 = 7 × 15 = |Z₇| × |K₆ matchings| (K₈ = K₆ ⊗ Z₇)
6. Z₇ trivial irrep absent from vacuum (symmetry must break)

## Connection to K₆

| Property | K₆ on T² | K₈ on Σ₂ |
|----------|----------|----------|
| Observable | Higgs mass (mH) | Fermion masses (Yukawa) |
| Action sector | Bosonic (a₄/a₂²) | Fermionic (vacuum eigenvectors) |
| Vacuum | Unique (λ=3.306) | 6D degenerate (λ=1.960) |
| Free parameters | 0 | 0 (K₆ vacuum fixes ρ₁ direction) |

**K₆ → mH = 125 GeV (bosonic). K₈ → Yukawa 415:135:1 (fermionic).**
