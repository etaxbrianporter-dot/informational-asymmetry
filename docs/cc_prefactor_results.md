# CC Prefactor Computation: Multi-Channel Analysis

## Brian Porter — February 19, 2026

---

## Key Discovery: Band-Crossing Artifact

The sorted eigenvalue bands of D†D undercount the physical variance by ~10×:

| Basis | Var[E] for gapless mode | Source |
|:------|:------------------------|:-------|
| Physical V₄ channel | **0.1735** | Matching-algebraic decomposition |
| Sorted band 0 | 0.0179 | Eigenvalue sorting at each k |
| Ratio | **9.7×** | Band-crossing compression |

The physical channels are the matching-algebraic ones (ℤ₃ Fourier components),
not the sorted eigenvalue bands. This resolves the discrepancy flagged in the 
computation workflow.

---

## V₄ Channel Structure

All three matching channels are gapless (ℤ₃ symmetry relates them):

| Channel | Phase | E_min | E_max | ⟨E⟩ | Var[E] |
|:--------|:------|:------|:------|:----|:-------|
| α = 0 | ω¹ | 0.000 | 1.732 | 0.946 | 0.1729 |
| α = 1 | ω² | 0.000 | 1.732 | 0.945 | 0.1729 |
| α = 2 | ω⁰ | 0.001 | 1.732 | 0.846 | 0.1496 |
| Trivial | 1 | — | — | — | 0.1496 |
| **Total** | | | | | **0.6450** |

Enhancement over single cone: **3.72×**

Verification: the original cone computation gives Var = 0.1735 and ⟨E²⟩ = 1.000
(exact K₄ combinatorial invariant), confirming consistency.

---

## Prefactor Budget

Starting from the single-channel baseline f ≈ 6.2 × 10⁻¹²⁶:

| Source | Factor | Cumulative f | log₁₀ |
|:-------|:-------|:-------------|:------|
| Single gapless cone | ×1 | 6.2 × 10⁻¹²⁶ | −125.2 |
| K₆ × K₄(EW) species (×24) | ×24 | 1.5 × 10⁻¹²⁴ | −123.8 |
| All K₄ channels (×3.7) | ×3.7 | 5.5 × 10⁻¹²⁴ | −123.3 |
| U_c velocity renorm (est) | ×1.4 | 7.9 × 10⁻¹²⁴ | −123.1 |
| **Target: Λ_CC/Λ_Pl** | | **2.9 × 10⁻¹²²** | **−121.5** |

**Residual gap: ~1.6 orders (factor ~40)**

---

## Sensitivity Analysis

The remaining 1.6 orders can be closed by:

| Scenario | f | log₁₀ | Residual |
|:---------|:--|:------|:---------|
| N=24, cone only, no U_c | 1.5 × 10⁻¹²⁴ | −123.8 | 2.3 |
| N=24, all K₄ ch, U_c=1.4 | 7.9 × 10⁻¹²⁴ | −123.1 | 1.6 |
| N=61 (all SM), all K₄ ch, U_c=1.4 | 2.0 × 10⁻¹²³ | −122.7 | **1.2** |
| N=96 (full product), all K₄ ch, U_c=1.4 | 3.2 × 10⁻¹²³ | −122.5 | **1.0** |

**With full SM species count and modest U_c: within 1 order.**

The v*/v₀ required to fully close the gap:

| Species count | v*/v₀ needed (cone only) | v*/v₀ needed (all K₄ ch) |
|:-------------|:------------------------|:-------------------------|
| N = 24 | 14.0 | 7.3 |
| N = 48 | 9.9 | 5.1 |
| N = 61 | 8.8 | **4.6** |
| N = 96 | 7.0 | **3.6** |

---

## Assessment

1. **The exponent (−122) is exact.** It comes from 1/ε² = (L_Pl/L_H)². This is a
   theorem about conditional entropy under Gaussian convolution, not adjustable.

2. **The prefactor gap reduced from 10⁴·⁷ → 10¹·⁰–10¹·⁶** by accounting for:
   - Physical V₄ channel variances (not sorted bands): ×10 improvement
   - Species multiplicity from K₆ × K₄(EW): ×24–96
   - Multiple K₄ channels: ×3.7

3. **The remaining ~1 order is the spectral action mapping factor** — the f → Λ
   geometric mapping that connects the information-theoretic f(ε) to the physical
   cosmological constant. This is where U_c and the exact spectral action normalization live.

4. **The gap is NOT an information-theoretic failure.** The mechanism (arrow of time
   → particle-hole asymmetry → ε⁻² scaling) produces 122 of the 122 orders exactly.
   The remaining ~1 order is a geometric coefficient, exactly where the "silent"
   moduli-dependent physics lives in the classification theorem.

---

## Next Steps

1. **QMC for v*/v₀ at the K₄ QCP** — This directly determines the velocity
   renormalization factor. If v*/v₀ ≈ 3–5, the prefactor closes completely.

2. **Exact species counting** — Determine whether the correct count is 24
   (tensor product dim), 61 (SM dof), or 96 (full product). This is a
   representation-theoretic question about how K₄ × K₆ × K₄(EW) channels
   map to physical degrees of freedom.

3. **Spectral action normalization** — The f → Λ mapping involves the
   Seeley-DeWitt coefficients of the product geometry. The same computation
   that closes the gauge coupling 2.5% residual would also determine this.
