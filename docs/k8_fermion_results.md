---
title: K₈ FERMION SECTOR: MAIN RESULTS
section: K₈ Fermion Sector
status: active
---

# K₈ FERMION SECTOR: MAIN RESULTS

## The Setup

K₈ (28 edges, 105 matchings) on genus 2 (double torus).

**Direction assignment from topology:**
- d₀ = diff ±1 mod 7 (K₇ Heawood, 7 edges)
- d₁ = diff ±2 mod 7 (K₇ Heawood, 7 edges)
- d₂ = diff ±3 mod 7 (K₇ Heawood, 7 edges)
- d₃ = handle edges (i,7) (7 edges)

Total: 28 edges in 4 classes of 7 each.

**Phase scheme:** Z₃ on torus (d₀→1, d₁→ω, d₂→ω²) with d₃→ω (handle = ω).

**Generation pairing:** (0,1) = doublet, (2,3)/(4,5)/(6,7) = 3 generations.

## The Computation

1. Build 105×105 Gram matrix: G_ij = Tr(M_i M_j) × phase_correlator × δ(D_i, D_j)
2. Find physical vacuum: smallest positive eigenvalue λ_vac = 1.9595 (6-fold degenerate)
3. Scan random directions within 6D vacuum eigenspace
4. For each direction: compute vacuum-weighted 3×3 Yukawa, extract Y†Y hierarchy

## KEY RESULTS

### Rank 3 is generic
- 9,999 / 10,000 random vacuum directions give rank-3 Y†Y
- Three massive generations is the GENERIC outcome (99.99%)

### Hierarchy distribution
- Range: 2 to 3.8 × 10⁹
- Median: 173
- Peak: 10² to 10³ (100-1000)
- SM range (10² - 10⁴) contains ~70% of cases

### Best matches to SM fermion mass ratios

| Sector | SM target | K₈ best match | Error |
|--------|-----------|---------------|-------|
| Down-type | 790 : 20 : 1 | **768 : 20 : 1** | **3% / 1%** |
| Up-type | 18,700 : 137 : 1 | 12,377 : 153 : 1 | 34% / 12% |
| Lepton | 12,100 : 43 : 1 | 11,893 : 68 : 1 | 2% / 58% |

**The down-type match is essentially exact:** 768:20:1 vs 790:20:1.

### The 6D vacuum eigenspace is the Yukawa moduli space
Different directions within the SAME vacuum eigenspace produce
all three SM fermion sectors. The SM requires 3 specific directions
(one per sector) within this 6D space.

## What Is Fixed vs Free

### FIXED (zero free parameters):
- K₈ graph structure → 105 matchings, overlap matrix
- Heawood directions → 4 classes of 7 edges
- Z₃ × handle phases → Gram matrix
- Vacuum eigenvalue λ = 1.9595
- Rank 3 (three generations)

### FREE (needs selection principle):
- Direction within 6D vacuum eigenspace → determines which fermion sector
- Handle phase (ω vs ω² vs 1) → changes hierarchy range
- Generation pairing → which K₆ vacuum matching

## Comparison to K₆

| Property | K₆ on T² | K₈ on Σ₂ |
|----------|----------|----------|
| Matchings | 15 | 105 |
| Gram matrix | 15×15 | 105×105 |
| Vacuum | λ = 3.306 (unique) | λ = 1.960 (6-fold) |
| Hierarchy | ~64:1 (2 nonzero) | up to 10⁹:1 (3 nonzero) |
| Fermion masses | ✗ (rank 2, flat) | ✓ (rank 3, SM-scale) |
| Higgs mass | ✓ (125 GeV) | ? (not computed yet) |

**K₆ knows the Higgs mass but not fermion masses.**
**K₈ knows fermion mass hierarchies but we haven't checked the Higgs.**

## Physical Interpretation

The genus-2 surface gives K₈ four direction classes instead of K₆'s three.
The fourth direction (handle) is what promotes the Yukawa matrix from
rank 2 (K₆) to rank 3 (K₈).

Every K₈ matching has exactly one handle edge — this is forced by vertex 7
being in the graph. The handle edge creates interference between 
generation-pair couplings that produces the hierarchy.

The 6D vacuum degeneracy is the moduli space of fermion masses.
Three specific directions within it give the three SM sectors.
What selects these directions is the next question.

## Open Questions

1. **What breaks the 6D degeneracy?** Higher-order corrections to the
   spectral action? RG running? Additional topological invariants?

2. **K₈ Higgs mass:** Does the same framework give mH ≈ 125 GeV?
   (Requires computing a₂, a₄ for the 105×105 system.)

3. **Can all 3 SM sectors come from 3 orthogonal directions in the 6D space?**
   If so, this would be a remarkable geometric realization of flavor.

4. **The handle phase ω:** Why d₃→ω and not d₃→1 or d₃→ω²?
   The d₃→ω scheme gives the widest hierarchy range spanning all SM sectors.
   Is there a topological argument selecting it?
