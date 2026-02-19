---
title: The Normalization Constant: Honest Assessment
section: K₆ Higgs Sector
status: active
---

# The Normalization Constant: Honest Assessment

## Kâ‚† Spectral Action â†’ Higgs Mass

**Date:** February 17, 2026  
**Status:** R is proved. c is not derived. No near-miss is privileged.  
**Supersedes:** all prior versions (v1â€“v8)

---

## 1. What Is Proved

| Quantity | Value | Status |
|----------|-------|--------|
| R = aâ‚„/aâ‚‚Â² | 0.3722127085 | Proved: Kâ‚† matching combinatorics |
| aâ‚‚ = Î»_min(G) | 3.3060051829 | Proved: Gram matrix eigenvalue |
| aâ‚„ = H(v,v,v,v) | 4.0681621736 | Proved: quartic form |
| mHÂ² = 8(R/c)mWÂ² | verified | c = 3 â†’ 80.08 GeV matches paper |
| n = 3 from saturation | unique | (2nâˆ’1)!! = n(2nâˆ’1) forces Kâ‚† |
| BZ traces k-independent | Var = 0 | At sorted vacuum |
| Gauge kinetic = aâ‚‚ | ratio = 1 | Unit hexagonal lattice |

These are exact, verified to 10-digit precision, and survive all cross-checks.

---

## 2. The One Undetermined Quantity

The Higgs mass formula contains one undetermined constant c:

    mH = âˆš(8R/c) Ã— mW

| c | mH (GeV) | Origin |
|---|----------|--------|
| 3.000 | 80.1 | Paper's normalization |
| 1.230 | 125.1 | Required by experiment |

The paper claims c = 3. Experiment requires c â‰ˆ 1.23. The factor of 2.4 between them is unexplained.

---

## 3. No Derivation of c Exists

### 3.1 The CCM trace formula (c = 4/aâ‚‚)

The mapping Î»/gÂ² = b/(4a) with a = aâ‚‚, b = aâ‚„ gives c = 4/aâ‚‚ = 1.210 and mH = 126.1 GeV. But at the democratic point, aâ‚„/(4aâ‚‚) = 0.429 while the known Î»/gÂ² = 7/48 = 0.146 â€” off by a factor of 2.94. The mapping fails its own cross-check.

### 3.2 The torus correction (c = Ï€Â²/8)

The decomposition Ï€Â²/8 = 3 Ã— Î¶(2)/4 assigns Î¶(2) to a "torus spectral correction." But the derivation_c.py computation showed explicitly that for a flat torus, the BZ-averaged spectral action gives the same normalization as pointwise evaluation. The flat TÂ² correction is trivial â€” Î¶(2) does not emerge from the spectral action on TÂ². Calling 4/Î¶(2) a "torus spectral invariant" is incorrect.

### 3.3 The 7/8 at democratic

The paper's democratic-point result R = 7/8 with c = 6 does not emerge from any Kâ‚†-only computation. BZ-averaged R gives 0.368; finite-space R gives 0.700; Kâ‚† eigenvalue R gives 7/16. The 7/8 requires the full SM fermion content via Kâ‚„(EW). The normalization c necessarily involves physics beyond Kâ‚†.

---

## 4. The Near-Misses Are Not Privileged

A systematic search over algebraic expressions built from {aâ‚‚, aâ‚„, R, Ï€, e, Î¶(2), Î¶(3), Î“(3/4), âˆš2, âˆš3, ln 2, integers} finds 41 expressions within 0.5% of c_exact â‰ˆ 1.2295. The closest:

| Expression | Value | Distance to c_exact |
|-----------|-------|---------------------|
| aâ‚„/aâ‚‚ | 1.2305 | 0.09% |
| Î“(3/4) | 1.2254 | 0.33% |
| Ï€Â²/8 | 1.2337 | 0.34% |
| âˆš(3/2) | 1.2247 | 0.39% |
| e/âˆš5 | 1.2157 | 1.12% |
| 4/aâ‚‚ | 1.2099 | 1.59% |

Ï€Â²/8 is not even the closest clean expression â€” aâ‚„/aâ‚‚ (a pure Kâ‚† ratio) is four times closer. Finding an algebraic expression within 0.5% of a target number near 1.23 is expected, not remarkable.

The "bracket" between 4/aâ‚‚ and Ï€Â²/8 is not a rigorous bound. It is two unproven approximations that happen to straddle experiment. Presenting it as a zero-parameter prediction is not warranted.

---

## 5. What the Framework Actually Determines

The Kâ‚† matching algebra determines one number: **R = 0.3722**. This is exact, combinatorial, and has no free parameters.

Combined with the spectral action formula, this gives:

    mH/mW = âˆš(8R/c) = 1.726/âˆšc

For any fixed c, this is a specific prediction of the mass ratio. The graph-theoretic constraint R âˆˆ [0.17, 0.50] (from 5000 random assignments) bounds mH for any c.

The framework does NOT determine c. The value c = 3 (paper) gives mH = 80 GeV. Either:

(a) The paper's c is wrong for the TÂ² geometry, and the correct c (â‰ˆ 1.23) has not been derived, or

(b) The framework predicts mH = 80 GeV and disagrees with experiment by a factor of 1.56.

Resolving this requires deriving the spectral action normalization for the fibered product geometry on TÂ², which has not been done.

---

## 6. One Numerical Fact Worth Noting

The Gram matrix minimum eigenvalue satisfies:

    aâ‚‚ = 3.3060 â‰ˆ 8mWÂ²/mHÂ² = 3.3032  (0.09% agreement)

This is equivalent to the formula mHÂ² â‰ˆ 2(aâ‚„/aâ‚‚)mWÂ², or c â‰ˆ 4/aâ‚‚. If a derivation of c = 4/aâ‚‚ existed, this would be an exact prediction relating a Kâ‚† graph eigenvalue to a ratio of SM masses. No such derivation exists â€” the one attempted (CCM b/(4a)) fails by a factor of 2.94 at the democratic point.

---

## 7. What Would Constitute Progress

To determine c from first principles, one must compute the Seeley-DeWitt expansion of the spectral action for the fibered product Kâ‚„(ST) Ã— Kâ‚† Ã— Kâ‚„(EW) on TÂ², and separately extract the gauge kinetic and Higgs quartic coefficients. Their ratio determines c. This is a well-defined finite computation that has not been performed.

---

## 8. Recommended Claim

> The Kâ‚† matching algebra determines R = aâ‚„/aâ‚‚Â² = 0.3722 from combinatorics. The spectral action formula mHÂ² = 8(R/c)mWÂ² then gives the Higgs mass for any normalization c. The standard normalization c = 3 gives mH = 80 GeV; the experimental value requires c â‰ˆ 1.23. Computing the spectral action normalization for the TÂ² geometry â€” which differs from the standard Mâ´ geometry â€” is the key open problem. The Gram matrix eigenvalue aâ‚‚ = 3.306 is within 0.1% of 8mWÂ²/mHÂ², a relation that would become exact if c = 4/aâ‚‚; deriving this from the spectral action remains open.
