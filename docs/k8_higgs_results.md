---
title: Kâ‚ˆ Spectral Action: Higgs Mass Check â€” Results
section: K₆ Higgs Sector
status: active
---

# Kâ‚ˆ Spectral Action: Higgs Mass Check â€” Results

**Date:** February 18, 2026
**Status:** COMPUTATION COMPLETE â€” verdict clear

---

## 1. What Was Computed

The spectral kurtosis Râ‚ˆ = aâ‚„/aâ‚‚Â² for the Kâ‚ˆ system (105 matchings, 105Ã—105 Gram matrix) on the genus-2 surface, using the Heawood + handle direction assignment and â„¤â‚ƒ phase scheme from Paper III.

**Method:** BZ-averaged Tr(Dâ€ D) and Tr((Dâ€ D)Â²) via 3D momentum grid (Nk=12, converged at Nk=10). Evaluated at all 6 vacuum eigenvectors and 10 random directions within the 6D vacuum eigenspace.

---

## 2. Gram Matrix Verification

All Paper III values reproduced exactly:

| Quantity | Paper III | Computed | Match |
|----------|-----------|----------|-------|
| # matchings | 105 | 105 | âœ“ |
| Overlap distribution | 3150:1680:630 | 3150:1680:630 | âœ“ |
| Zero eigenvalues | 5 | 5 | âœ“ |
| Î»_vac | 1.9595 | 1.9595 | âœ“ |
| Vacuum multiplicity | 6 | 6 | âœ“ |
| Î»_max | 24.0 | 24.0 | âœ“ |

---

## 3. Main Result: Râ‚ˆ â‰  Râ‚†

**Kâ‚ˆ does not reproduce the Kâ‚† Higgs mass.**

| System | aâ‚‚ (= Î»_vac) | R = aâ‚„/aâ‚‚Â² | m_H (c=3) | m_H (câ‰ˆ1.23) |
|--------|---------------|-------------|-----------|---------------|
| Kâ‚† (standalone) | 3.306 | 0.3722 | 80.1 GeV | 124.9 GeV |
| Kâ‚ˆ vacuum (mean) | 1.960 | 0.234 | 63.5 GeV | 99.0 GeV |
| Kâ‚ˆ vacuum (range) | 1.960 | 0.197â€“0.249 | 58â€“66 GeV | 91â€“102 GeV |
| Kâ‚ˆ at Kâ‚†-projection | 1.960 | 0.197 | 58.3 GeV | 90.8 GeV |

The Kâ‚ˆ spectral kurtosis is ~63% of Kâ‚†'s value. No direction within the Kâ‚ˆ vacuum eigenspace reaches Râ‚† = 0.3722.

---

## 4. R Varies Across the Vacuum Eigenspace

This is a new finding not anticipated in Paper III.

The individual vacuum eigenvector R values:

| Eigenvector | R = aâ‚„/aâ‚‚Â² |
|-------------|-------------|
| vâ‚€ | 0.2225 |
| vâ‚ | 0.1970 |
| vâ‚‚ | 0.2334 |
| vâ‚ƒ | 0.2171 |
| vâ‚„ | 0.2390 |
| vâ‚… | 0.2380 |

Random directions: R âˆˆ [0.206, 0.249], mean = 0.234 Â± 0.013.

**The Higgs mass depends on the Yukawa direction.** If the fermion hierarchy selects a specific direction in the 6D vacuum eigenspace (as Paper III argues), that direction simultaneously determines Râ‚ˆ and hence predicts a Higgs mass. This is a coupling between the Higgs sector and the Yukawa sector that doesn't exist in Kâ‚†.

---

## 5. The Hub Sector Anomaly

The 15 hub matchings (containing edge (0,1)) form a sub-Gram matrix whose ground eigenvalue is Î»_hub = 2.764, not the standalone Kâ‚† value of 3.306. The hub R = 0.395, not 0.372.

**Why the difference:** The Kâ‚ˆ hub Gram matrix uses 4-component direction vectors (including the handle), while standalone Kâ‚† uses 3-component vectors. This changes the momentum conservation constraint Î´(D_i, D_j). The handle edge in every hub matching adds a common dâ‚ƒ = (0,0,1,0) component, but this changes which pairs of matchings satisfy D_i = D_j.

**This means Kâ‚† and Kâ‚ˆ are genuinely different systems, not related by simple embedding.**

---

## 6. Kâ‚† Vacuum Projection Confirmed

|Î _vac v_Kâ‚†| = 0.2355 (5.55% of norm squared), confirming Paper III Theorem 5.

The projection is concentrated in specific â„¤â‚‡ components:
- Components (0, ~0, 0.118, -0.078, 0.125, 0.125) in the 6D basis
- vâ‚ component essentially zero â†’ confirms Ïâ‚‚ inaccessibility theorem

R at this projected direction = 0.197, the LOWEST value in the vacuum eigenspace.

---

## 7. Interpretation

### Kâ‚ˆ does NOT supersede Kâ‚† for the Higgs mass.

Three readings:

**(A) Sector separation (most likely):**
Kâ‚† determines the Higgs potential; Kâ‚ˆ determines the Yukawa couplings. They're different levels of the Kâ‚„ â†’ Kâ‚† â†’ Kâ‚ˆ hierarchy, each computing different SM parameters. The Higgs mass comes from Kâ‚† alone (Râ‚† = 0.3722), and Kâ‚ˆ's different Râ‚ˆ is irrelevant to the Higgs â€” it characterizes the quartic structure of the Yukawa sector instead.

This is consistent with the existing framework: Paper I establishes Kâ‚„ (spacetime), Paper II adds Kâ‚† (gauge + Higgs), Paper III adds Kâ‚ˆ (fermion). Each level adds structure without overwriting the previous level.

**(B) Coupled prediction (speculative but testable):**
If Kâ‚ˆ IS the fundamental object, then Râ‚ˆ depends on the Yukawa direction, creating a Higgs-Yukawa coupling. The down-type direction (which gives the best SM match at 768:20:1) selects a specific Râ‚ˆ and hence a specific m_H. Computing Râ‚ˆ at the exact down-type direction would give a Higgs mass prediction that's testable against experiment.

**(C) Normalization mismatch:**
Kâ‚ˆ requires câ‚ˆ â‰ˆ 0.77 for m_H = 125 GeV, vs Kâ‚†'s câ‚† â‰ˆ 1.23. The ratio câ‚ˆ/câ‚† â‰ˆ 0.63 = Râ‚ˆ/Râ‚†. Whether there's a principled reason for Kâ‚ˆ to have a different normalization (perhaps related to the genus-2 surface vs torus geometry) is an open question.

---

## 8. What This Changes

### Paper III should note:

1. **Râ‚ˆ â‰ˆ 0.234 â‰  Râ‚† = 0.372.** Kâ‚ˆ does not independently reproduce the Higgs mass.

2. **Râ‚ˆ varies across the vacuum eigenspace (0.197â€“0.249).** The quartic structure is not a vacuum invariant â€” it depends on the Yukawa direction.

3. **The Kâ‚† projection direction (Ïâ‚ only) gives the LOWEST Râ‚ˆ = 0.197.** This is interesting: the Higgs-sector projection minimizes the quartic coupling within the Kâ‚ˆ vacuum.

4. **Sector separation is the natural reading:** Kâ‚† = Higgs mass, Kâ‚ˆ = fermion masses. The two computations are complementary, not redundant.

### For the c = aâ‚„/aâ‚‚ termination analysis (Paper I, P5):

The c termination analysis noted that "whether Kâ‚ˆ reproduces the same Higgs mass from its own spectral data remains open." **Answer: it does not.** This confirms that the Higgs mass is a Kâ‚† result, not a Kâ‚ˆ result. Paper III's fermion sector is genuinely new physics beyond Kâ‚†, not a recalculation of Kâ‚† physics.

---

## 9. New Open Questions

**Q1.** What is Râ‚ˆ at the exact down-type best-match direction (768:20:1)? If interpretation (B) is correct, this predicts a specific m_H.

**Q2.** Is Râ‚ˆ = 0.197 at the Ïâ‚ direction significant? The minimum-quartic direction is exactly the Higgs-sector direction. Could this reflect a variational principle selecting the physical Higgs coupling?

**Q3.** Does Râ‚ˆ have a maximum at some special direction? The maximum Râ‚ˆ â‰ˆ 0.249 gives m_H â‰ˆ 102 GeV with c â‰ˆ 1.23 â€” still short of 125 GeV. So no direction in the Kâ‚ˆ vacuum reaches the Kâ‚† value, even with the most favorable normalization.
