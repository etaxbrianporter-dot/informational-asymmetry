---
title: The Normalization Constant: What Is Actually Proved
section: K₆ Higgs Sector
status: active
---

# The Normalization Constant: What Is Actually Proved

## Kâ‚† Spectral Action â†’ Higgs Mass: Honest Audit

**Date:** February 17, 2026  
**Status:** All paper values reproduced. No first-principles derivation of c exists.  
**Supersedes:** normalization_derivation_v2 through v7

---

## 1. Hard Facts (proved, verified to 10 digits)

The Gram matrix formula G_ij = Tr(Máµ¢ Mâ±¼) Ã— Re(Ï‰^(zâ‚ƒâ±¼ âˆ’ zâ‚ƒáµ¢)) Ã— Î´(diráµ¢, dirâ±¼) with the sorted assignment reproduces every quantity in the paper exactly.

| Quantity | Value | Status |
|----------|-------|--------|
| R = aâ‚„/aâ‚‚Â² | 0.3722127085 | Proved: pure Kâ‚† combinatorics |
| aâ‚‚ = Î»_min(G) | 3.3060051829 | Proved: Gram matrix eigenvalue |
| aâ‚„ = H(v,v,v,v) | 4.0681621736 | Proved: quartic form on ground eigenvector |
| mHÂ² = 8(R/c)mWÂ² | verified | Proved: c=3 â†’ 80.08 GeV matches paper |
| n = 3 from saturation | unique solution | Proved: (2nâˆ’1)!! = n(2nâˆ’1) |
| BZ traces k-independent at vacuum | Var(Sâ‚‚) = 0 | Proved: all vacuum matchings share dir |
| Gauge kinetic = aâ‚‚ | ratio = 1.000000 | Proved: |dáµ¢|Â² = 1 on unit lattice |
| Internal gauge curvature = 0 | Fâ‚â‚‚ = 0 | Proved: abelian vacuum |

---

## 2. The Normalization Problem

The Higgs mass formula contains one undetermined constant:

    mHÂ² = 8 Â· (R/c) Â· mWÂ²

R = 0.3722 is proved. mW = 80.379 GeV is experimental. The value of c determines everything:

| c | mH (GeV) | Source |
|---|----------|--------|
| 6.000 | 56.6 | Standard CCM on Mâ´ Ã— F |
| 3.000 | 80.1 | Paper's stated normalization |
| **1.230** | **125.1** | Required by experiment |
| 1.234 (= Ï€Â²/8) | 124.9 | Numerical near-miss |
| 1.210 (= 4/aâ‚‚) | 126.1 | Numerical near-miss |

---

## 3. Why Neither "Derivation" of c Survives

### 3.1 The c = 4/aâ‚‚ claim (mHÂ² = 2Â·aâ‚„/aâ‚‚Â·mWÂ²)

**Stated derivation:** Map CCM trace formula Î»/gÂ² = b/(4a) onto Kâ‚† via a_CCM = aâ‚‚, b_CCM = aâ‚„, giving Î»/gÂ² = aâ‚„/(4aâ‚‚).

**Cross-check that kills it:** At the democratic point, aâ‚„/(4aâ‚‚) = 0.429, but the known Î»/gÂ² = 7/48 = 0.146. The mapping is off by a factor of 2.94. The a_CCM â†’ aâ‚‚ identification is wrong â€” CCM traces are sums of squared Yukawa eigenvalues, not traces of the 6Ã—6 Kâ‚† Dirac matrix.

**Status:** No surviving derivation. The formula mHÂ² = 2(aâ‚„/aâ‚‚)mWÂ² = 126.1 GeV is a numerical observation, not a derived result.

### 3.2 The c = Ï€Â²/8 claim (= 3 Ã— Î¶(2)/4)

**Stated derivation:** Decompose c into N_c Ã— Î¶(2)/dim(coset) where N_c = 3 from Kâ‚† matching algebra, Î¶(2) = Ï€Â²/6 from flat torus spectral theory, and 4 from electroweak coset dimension.

**Three problems:**

1. **The "3" is undetermined.** N_c and N_gen are computationally degenerate â€” both equal 3 from the same saturation equation. The complexifier needed to identify the 3 as specifically "color" is not vacuum-selected: the Pati-Salam complexifier (M07) has vacuum overlap 0.45, while the vacuum-preferred M03 gives so(4)âŠ•so(2), not su(3)âŠ•u(1).

2. **The Î¶(2)/4 factor lacks derivation.** The claim that the torus spectral invariant enters the normalization as c â†’ c Ã— Î¶(2)/4 has not been derived from the spectral action expansion on TÂ². It is a conjecture motivated by the numerical coincidence.

3. **Not universal across moduli space.** At generic random assignments, R drops to 0.26â€“0.33, giving mH = 104â€“118 GeV with c = Ï€Â²/8 fixed. The near-125 result at both democratic and sorted is because R â‰ˆ 0.37 at both, not because c = Ï€Â²/8 is structurally preferred.

**Status:** Algebraically suggestive but physically unjustified. No derivation from the spectral action on TÂ² exists.

### 3.3 The "bracket" claim

The statement "mH âˆˆ [124.9, 126.1]" presents two unproven approximations as bounds. Neither endpoint has a first-principles derivation. That they happen to straddle the experimental value is a numerical observation, not a theorem. This should not be presented as a zero-parameter prediction.

---

## 4. What the 7/8 Puzzle Reveals

The paper's democratic-point result R = 7/8, c = 6, Î»/gÂ² = 7/48 does not emerge from Kâ‚† alone:

| Computation | R value |
|-------------|---------|
| BZ-averaged with Zâ‚ƒ phases | 0.368 |
| Finite-space (no phases, no BZ) | 0.700 |
| Kâ‚† Dâ€ D eigenvalues (3Â±âˆš3 pairs) | 7/16 |
| Yukawa eigenvalue ratio | 1/3 |
| Paper's claim | 7/8 |

The 7/8 requires the full SM fermion content â€” quarks Ã— N_c + leptons with electroweak quantum numbers. This is Kâ‚„(EW) sector physics, not Kâ‚† physics. The normalization c = 6 = 2N_c is specific to the Mâ´ spectral action with the SM finite space.

**Implication:** The relationship between the Kâ‚†-on-TÂ² computation (R_BZ = 0.372) and the SM Higgs mass necessarily involves the Kâ‚„(EW) sector. The normalization c encodes this relationship. Determining c from first principles requires the full product geometry computation.

---

## 5. What IS the Genuine Predictive Content

### 5.1 Graph-theoretic constraint on R

The 120:90:15 overlap distribution of Kâ‚† matchings constrains R to the interval [0.17, 0.50] over all possible assignments. Over 5000 random samples: mean R = 0.372, std = 0.061. The sorted vacuum gives R = 0.372, near the distribution mean â€” not a fine-tuned value.

This constraint IS proved and IS a prediction: for any fixed c, the Higgs mass is bounded by Kâ‚† combinatorics.

### 5.2 The ratio mH/mW

If c is determined by the geometry (whatever its value), then:

    mH/mW = âˆš(8R/c)

R is fixed by Kâ‚†. So mH/mW is a pure prediction once c is known. The framework predicts a SPECIFIC ratio, not a range.

### 5.3 Structural coincidences

Several features are genuinely remarkable and non-trivial:

- c_exact = 1.2295 is within 0.34% of Ï€Â²/8 and 1.6% of 4/aâ‚‚
- Both candidates have clean algebraic forms involving only Kâ‚† data and standard constants
- The reduction factor from c = 3 to c_exact is 4/Î¶(2) â‰ˆ 2.43, which is a standard spectral quantity of flat tori
- R at the sorted vacuum (0.372) is remarkably stable and near the generic mean

These coincidences are suggestive of real structure. But "suggestive" â‰  "derived."

---

## 6. What Would Close the Gap

To derive c from first principles, at least one of these computations must be performed:

**Option A:** Seeley-DeWitt expansion of Tr(f(DÂ²/Î›Â²)) on TÂ² Ã— Kâ‚†, extracting the gauge-kinetic and Higgs-quartic coefficients separately and forming their ratio.

**Option B:** Full Kâ‚„(ST) Ã— Kâ‚† Ã— Kâ‚„(EW) product Dirac operator decomposition. The 3 Kâ‚„(EW) generators act on the 15-dimensional Kâ‚† moduli space. Decompose Tr(D_productâ´) into FÂ² and |H|â´ sectors.

**Option C:** Lattice gauge theory simulation of the Kâ‚† gauge-Higgs system on TÂ², measuring Î» and gÂ² independently.

All three are well-defined finite computations. None has been performed.

---

## 7. Honest Summary

**The framework determines:**
- R = aâ‚„/aâ‚‚Â² = 0.3722 from Kâ‚† combinatorics (proved)
- The formula mHÂ² = 8(R/c)mWÂ² (verified)
- n = 3 uniquely from saturation (proved)
- R âˆˆ [0.17, 0.50] from graph theory (proved)

**The framework does not determine:**
- The normalization constant c
- Whether c = Ï€Â²/8 or c = 4/aâ‚‚ or something else
- The physical identification of the "3" as color vs generation

**The honest claim:**
IF c â‰ˆ 1.23, THEN mH â‰ˆ 125 GeV. The value c â‰ˆ 1.23 has two attractive near-misses but no first-principles derivation. The paper's c = 3 gives 80 GeV; experiment requires c â‰ˆ 1.23; the factor ~2.4 between them is unexplained.

**What is remarkable (but not proved):**
The required c_exact = 1.2295 sits between Ï€Â²/8 = 1.2337 and 4/aâ‚‚ = 1.2099 â€” two expressions with clean algebraic content. The reduction from the paper's c = 3 is by the factor 4/Î¶(2), a standard torus spectral invariant. These numerical coincidences deserve investigation but do not constitute a derivation.

---

## 8. Recommended Presentation

For the paper, the claim should be:

> The Kâ‚† spectral action determines R = aâ‚„/aâ‚‚Â² = 0.3722 from matching combinatorics. Combined with the spectral action formula mHÂ² = 8(R/c)mWÂ², this yields mH = 125 GeV for c â‰ˆ 1.23. The standard CCM normalization c = 3 gives mH = 80 GeV. The resolution requires deriving the spectral action normalization on TÂ² rather than Mâ´. The value c = Ï€Â²/8 = 3Î¶(2)/4, which would give mH = 124.9 GeV, is motivated by the torus spectral zeta function but has not been derived from first principles. Determining c from the TÂ² Ã— Kâ‚† spectral action is the key remaining computation.

This is weaker than "zero-parameter prediction" but stronger than "postdiction." The framework constrains mH to a narrow range *given* c, and the required c has attractive algebraic structure suggestive of the torus geometry.
