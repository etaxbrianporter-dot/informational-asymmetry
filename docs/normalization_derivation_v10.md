---
title: The Normalization Constant: Status Report
section: K₆ Higgs Sector
status: active
---

# The Normalization Constant: Status Report

## Kâ‚† Spectral Action and the Higgs Mass

**Date:** February 17, 2026  
**Supersedes:** all prior versions (v1â€“v9)

---

## 1. What Is Proved

The Kâ‚† matching algebra with the sorted assignment and Zâ‚ƒ phase structure determines three spectral invariants exactly:

    aâ‚‚ = Î»_min(G) = 3.3060051829
    aâ‚„ = H(v,v,v,v) = 4.0681621736
    R  = aâ‚„/aâ‚‚Â²    = 0.3722127085

These are pure combinatorics â€” computed from the 15 perfect matchings of Kâ‚†, the 120:90:15 overlap distribution, the hexagonal lattice direction assignment, and the Zâ‚ƒ phase structure. Verified to 10-digit precision against the paper.

Additional proved results at the sorted vacuum:

- BZ traces are k-independent (Var(Sâ‚‚(k)) = 0). All vacuum matchings share direction dâ‚€, so phase factors e^{ikÂ·(dâ±¼âˆ’dáµ¢)} = 1.
- Gauge kinetic trace equals aâ‚‚ exactly. Follows from |dáµ¢|Â² = 1 on the unit hexagonal lattice.
- Internal gauge curvature vanishes. The commutator [âˆ‚D/âˆ‚kâ‚, âˆ‚D/âˆ‚kâ‚‚] = 0 â€” the vacuum is abelian.
- n = 3 is the unique solution to (2nâˆ’1)!! = n(2nâˆ’1). This forces Kâ‚† as the matching algebra.
- Over random assignments, R âˆˆ [0.17, 0.50] with mean 0.372 and std 0.061. The sorted vacuum R = 0.372 is near the mean â€” combinatorially generic, not fine-tuned.

---

## 2. The Higgs Mass Formula

The paper uses:

    mHÂ² = 8 Â· (R/c) Â· mWÂ²

where mW = 80.379 GeV is experimental input and c is a normalization constant from the NCG framework. This formula is verified by the check c = 3 â†’ mH = 80.08 GeV (matching the paper's stated result).

The formula originates from the Chamseddine-Connes-Marcolli (CCM) spectral action on Mâ´ Ã— F, where:

- mHÂ² = 2Î»vÂ² (Higgs mechanism)
- v = 2mW/g (electroweak VEV)
- Î»/gÂ² = R/c (spectral action ratio)
- Therefore mHÂ² = 8(R/c)mWÂ²

This derivation assumes the standard 4D Higgs mechanism with electroweak symmetry breaking pattern, matched to the spectral action via specific trace identities. Whether this formula holds for Kâ‚† fibered over TÂ² â€” rather than a discrete finite space F on Mâ´ â€” is assumed, not derived.

---

## 3. The Normalization Problem

With R = 0.3722 fixed:

| c | mH (GeV) | Source |
|---|----------|--------|
| 6.000 | 56.6 | CCM standard: c = 2N_c |
| 3.000 | 80.1 | Paper's normalization |
| 1.230 | 125.1 | Required by experiment |

The experimental Higgs mass requires c â‰ˆ 1.23. The paper uses c = 3. The factor of 2.4 between them is unexplained.

No first-principles derivation of c exists for Kâ‚† on TÂ².

---

## 4. Failed Derivation Attempts

### 4.1 CCM trace formula: c = 4/aâ‚‚

**Idea:** Map the CCM formula Î»/gÂ² = b/(4a) onto Kâ‚† via a â†’ aâ‚‚, b â†’ aâ‚„, giving c = 4/aâ‚‚ = 1.210 and mH = 126.1 GeV.

**Cross-check:** At the democratic point, aâ‚„/(4aâ‚‚) = 0.429, but the known Î»/gÂ² = 7/48 = 0.146. Off by a factor of 2.94.

**Diagnosis:** The CCM quantities a = Î£|Yáµ¢|Â² and b = Î£|Yáµ¢|â´ are sums over Yukawa eigenvalues. The Kâ‚† quantities aâ‚‚ = Tr(Dâ€ D) and aâ‚„ = Tr((Dâ€ D)Â²) are matrix traces including off-diagonal cross-terms. These are different objects; the mapping fails.

### 4.2 Torus spectral correction: c = Ï€Â²/8

**Idea:** Decompose Ï€Â²/8 = 3 Ã— Î¶(2)/4 and identify the Î¶(2) factor as a spectral correction from the flat torus TÂ².

**Cross-check:** The derivation_c.py computation showed explicitly that for a flat torus, the BZ-averaged spectral action gives the same normalization as pointwise evaluation. The correction from TÂ² is trivial. Î¶(2) does not emerge from the spectral action on TÂ².

**Diagnosis:** Î¶(2) = Ï€Â²/6 is a number theory constant (the Basel sum Î£1/nÂ²). It is not a torus spectral invariant. The heat kernel on TÂ² has finite part 1, not Î¶(2). The Epstein zeta function at s = 1 is lattice-dependent and equals â‰ˆ 5.17 for the hexagonal lattice, not Ï€Â²/6 â‰ˆ 1.64. Calling 4/Î¶(2) a "torus spectral invariant" is incorrect.

### 4.3 The 7/8 at democratic

**Idea:** The paper's democratic-point value R = 7/8 should emerge from Kâ‚† and determine the normalization.

**Cross-check:** No Kâ‚†-only computation gives 7/8. BZ-averaged R = 0.368; finite-space R = 0.700; Kâ‚† eigenvalue R = 7/16; Yukawa eigenvalue R = 1/3. The 7/8 requires the full SM fermion content including the Kâ‚„(EW) sector â€” quarks Ã— N_c + leptons with hypercharge normalization.

**Diagnosis:** The normalization involves physics beyond Kâ‚†. The Kâ‚† computation alone cannot determine c.

---

## 5. The Numerical Near-Misses

Several algebraic expressions are close to c_exact = 1.2295:

| Expression | Value | Distance |
|-----------|-------|----------|
| aâ‚„/aâ‚‚ | 1.2305 | 0.09% |
| Î“(3/4) | 1.2254 | 0.33% |
| Ï€Â²/8 | 1.2337 | 0.34% |
| âˆš(3/2) | 1.2247 | 0.39% |
| 4/aâ‚‚ | 1.2099 | 1.59% |

A systematic search over simple algebraic expressions finds 41 within 0.5% of c_exact. Ï€Â²/8 is not the closest â€” aâ‚„/aâ‚‚ is four times closer. Finding a "clean" expression within 0.5% of any target near 1.23 is expected given the density of such expressions.

The near-miss aâ‚‚ â‰ˆ 8mWÂ²/mHÂ² (0.09% agreement) is equivalent to c â‰ˆ 4/aâ‚‚. This would relate a Kâ‚† graph eigenvalue to a ratio of SM masses â€” if a derivation existed. None does.

---

## 6. The N_c / N_gen Degeneracy

The saturation equation (2nâˆ’1)!! = n(2nâˆ’1) forces n = 3 uniquely, producing N_c = 3 (via complexifier â„â¶ â†’ â„‚Â³), N_gen = 3 (via Dâ‚† eigenvalue pairs), and |V(Kâ‚†)|/2 = 3 â€” all from the same structural source.

The complexifier needed to identify the "3" as specifically "color" (Pati-Salam su(3)âŠ•u(1) decomposition) requires matching M07. But the vacuum preferentially aligns with M03 (overlap 0.77 vs 0.45), which gives the electroweak-type so(4)âŠ•so(2) decomposition. The color structure is a choice, not a consequence of the vacuum.

---

## 7. Open Problems

Determining c from first principles requires resolving four stacked problems, each substantial:

**Problem 1 â€” Formulation.** Specify the product Kâ‚„(EW) Ã— Kâ‚† Dirac operator explicitly. The paper describes this conceptually but does not give the action of Kâ‚„(EW) generators on Kâ‚† moduli. This requires embedding the 3 Kâ‚„ matchings (as SU(2) generators) in the 15-dimensional Kâ‚† matching space with definite grading and real structure.

**Problem 2 â€” Theory.** Extend the Seeley-DeWitt expansion (or find an alternative) for discrete-continuous fibered geometries. The standard theorem applies to smooth manifolds; a discrete graph complement fibered over a continuous torus is a hybrid structure without standard heat kernel expansion. The derivation_c.py computation showed that treating TÂ² as purely "internal" (d = 4) gives a trivial correction, while treating it as physical dimensions (d = 6) conflates the BZ with physical space.

**Problem 3 â€” Structure.** Identify gauge and Higgs sectors within Kâ‚† moduli. This requires both the Kâ‚„(EW) action (Problem 1) and the complexifier choice (not vacuum-selected). The democratic cross-check showed the CCM separation of gauge and Higgs traces does not map onto Kâ‚† BZ-averaged traces.

**Problem 4 â€” Derivation.** Show the formula mHÂ² = 8(R/c)mWÂ² holds for Kâ‚† on TÂ², not just for CCM on Mâ´ Ã— F. The formula is currently imported from the Mâ´ framework and applied to Kâ‚† quantities.

Of these, Problem 1 is a tractable finite algebraic question. Problem 2 is a theoretical challenge. Problems 3 and 4 depend on the first two.

This is a research program, not a remaining calculation.

---

## 8. Summary

**The framework determines:** R = 0.3722 from Kâ‚† matching combinatorics. This is exact, has zero free parameters, and is robust across the moduli space.

**The framework does not determine:** The normalization c. The paper's c = 3 gives mH = 80 GeV. Experiment requires c â‰ˆ 1.23. No derivation of c for the TÂ² geometry exists, and the path to one involves several open problems in noncommutative geometry.

**The honest claim:** The Kâ‚† matching algebra constrains the spectral kurtosis R = 0.3722 from pure combinatorics. Combined with the CCM spectral action formula (assumed to hold for this geometry), this gives mH = 1.726/âˆšc Ã— mW. The paper's normalization c = 3 yields mH = 80 GeV; the experimental value requires c â‰ˆ 1.23. Resolving this discrepancy requires either deriving c for the TÂ² fibered geometry â€” which is an open research problem â€” or accepting that the framework in its current form does not predict the Higgs mass.
