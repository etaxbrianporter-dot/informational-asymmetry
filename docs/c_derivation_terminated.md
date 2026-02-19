---
title: c = a₄/a₂: Computation to Termination
section: K₆ Higgs Sector
status: active
---

# c = a₄/a₂: Computation to Termination

## The Question
Does the spectral action on K₄(ST) × K₆ × K₄(EW) determine the normalization constant c, and if so, does it equal a₄/a₂?

## What Was Computed

### 1. Sign Convention Resolution
The paper's K₆ matching matrices are **symmetric** (M[a,b] = M[b,a] = 1), not antisymmetric. This was identified by testing both conventions against the paper's known values:

| Convention | a₂ (Gram) | a₄ (quartic) | Match paper? |
|---|---|---|---|
| Antisymmetric | 3.3060 ✓ | 3.5366 ✗ | Partial |
| **Symmetric** | **3.3060 ✓** | **4.0682 ✓** | **Both** |

With symmetric matching matrices, Tr(D†D) = a₂ and Tr((D†D)²) = a₄ **exactly**, because the quartic form H(v,v,v,v) = Σ v_i v_j v_k v_l Tr(M_i M_j M_k M_l) × phase coincides with the matrix trace when M_i are symmetric and real.

### 2. The Numerical Result
With verified K₆ spectral data:
- a₂ = 3.3060051829
- a₄ = 4.0681621736  
- a₄/a₂ = 1.2305371433

Setting c = a₄/a₂ gives:
- mH = mW × √(8/a₂) = 80.379 × √(8/3.306) = **125.04 GeV**
- Experiment: 125.09 ± 0.11 GeV
- Deviation: −0.05 GeV (0.5σ)

This is closer than any other candidate:

| c | Value | mH (GeV) | Δ from exp |
|---|---|---|---|
| **a₄/a₂** | **1.2305** | **125.04** | **−0.05** |
| π²/8 | 1.2337 | 124.88 | −0.21 |
| 4/a₂ | 1.2099 | 126.10 | +1.01 |
| c_exact | 1.2295 | 125.09 | 0.00 |

### 3. Tests for Algebraic Identity

**Test: Is a₄/a₂ constant for all eigenvectors of G?**
No. It ranges from 1.00 (K₄-like modes) to 2.19 across the 15 eigenvectors. It is specific to the vacuum (minimum eigenvector).

**Test: Is H(v,·,v,·) proportional to G?**
No. The ratio matrix H(v,·,v,·)/G has eigenvalue ratios ranging from −0.62 to +0.47. The quartic tensor is structurally distinct from G⊗G.

**Test: Is it universal for random Gram/quartic pairs?**
No. For random symmetric matrices, a₄/a₂ at the minimum eigenvector has no preferred value (mean ≈ 0, std ≈ 6).

**Test: Does it hold for other complete graph geometries?**
- K₄: a₄/a₂ = 1.000 (trivially — all eigenvalues equal)
- K₆: a₄/a₂ = 1.231
- K₈: vacuum is degenerate (a₂ → 0), ratio undefined

**Test: Is a₄/a₂ at the vacuum an extremum?**
No. Over 1000 random unit vectors in the dir=0 subspace: a₄/a₂ ranges from 0.91 to 2.19 with mean 1.48. The vacuum value 1.23 falls within the distribution, not at any boundary.

### 4. Product Geometry Extraction

The naive product D_F = D₆ ⊗ I₄ + I₆ ⊗ D₄ gives, for the K₄(EW) Ψ₃ doublet decomposition:
- Yukawa R_Y = 0.794
- With N_c = 3: λ/g² = 0.132, c = 2.81, mH = 82.7 GeV

This is **wrong by a factor of ~2.3**. The product geometry with standard color multiplicity gives c ≈ 2.8, not 1.2.

The value c = a₄/a₂ ≈ 1.23 would correspond to an effective N_c ≈ 1.3 — between "no color" (N_c = 1) and "full QCD" (N_c = 3).

### 5. The K₄ Algebra Issue

Symmetric K₄ matching matrices commute: Ψ_a Ψ_b = Ψ_b Ψ_a. They form the Klein four-group V₄ = Z₂ × Z₂, not SU(2). This is true for BOTH symmetric and antisymmetric representations — the matchings are commuting involutions in either case.

This means: the K₄(EW) matchings do NOT generate the electroweak SU(2). They provide a Z₂ × Z₂ grading that defines doublet structure, but the non-abelian gauge symmetry requires additional structure beyond perfect matchings.

## Diagnosis

**Why c = a₄/a₂ cannot be derived from the spectral action:**

1. **c is not a spectral invariant of K₆ alone.** It depends on the product structure K₄(EW) × K₆, specifically on how the gauge and Higgs sectors are separated. Different embeddings of K₄ into K₆ give different values of c.

2. **The spectral action determines mH² = 8(λ/g²)mW², but λ/g² requires decomposing D_F into gauge and Higgs parts.** This decomposition is not intrinsic to K₆ — it comes from the real structure J and the product geometry.

3. **For c = a₄/a₂, the required identity is λ/g² = 1/Tr(D†D).** This would mean [Higgs quartic] × [gauge quadratic] = 1/a₂. No spectral action theorem gives this — it would relate quantities from different sectors of the spectral triple.

4. **The product geometry gives c ≈ 2.8 (with color), not 1.2.** Getting c ≈ 1.2 requires either dropping color entirely (but then c = 1, not 1.2) or having a partial color factor (N_c ≈ 1.3) that has no standard interpretation.

## The Honest Conclusion

**c = a₄/a₂ is a numerical observation, not a derivation.**

It gives mH = 125.04 GeV (within 0.5σ of experiment), which is better than any other candidate formula. But:

- It is not an algebraic identity of the K₆ Gram/quartic structure
- It is not a consequence of any known spectral action theorem
- It is not reproduced by the naive product geometry computation
- The 0.05 GeV deviation is within experimental uncertainty but is also consistent with the relation being approximate

**What would be needed to promote it to a theorem:**

Either (a) a structural identity relating Tr((D†D)²)/Tr(D†D) to the gauge-Higgs decomposition ratio for the specific embedding K₄ ↪ K₆, or (b) a new normalization principle in the spectral action that selects c = Σλᵢ²/Σλᵢ as the correct "spectral weight" for the Higgs quartic.

Neither exists in the current literature or follows from the computations performed here.

**Status: TERMINATED.** The question has a definite negative answer: no spectral action argument produces c = a₄/a₂. The numerical proximity (0.086%) remains unexplained but cannot be elevated to a derivation without new mathematical input.
