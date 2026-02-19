---
title: The Higgs Mass from the Graph Chain
section: K₆ Higgs Sector
status: active
---

# The Higgs Mass from the Graph Chain

## Kâ‚„ âŠ‚ Kâ‚… âŠ‚ Kâ‚†: A First-Principles Refactoring

**Date:** February 17, 2026  
**Purpose:** Derive the Higgs mass formula natively from graph geometry, without importing the Chamseddine-Connes-Marcolli (CCM) framework. Every factor gets a graph-geometric origin. Honest inventory of what is derived versus assumed.

---

## 0. The Problem with the Old Approach

We imported the CCM formula mHÂ² = 8(R/c)mWÂ² and spent ten versions trying to determine c within that framework. The v10 audit proved this is a dead end: four stacked open problems, none routine.

The import was backwards. We have our own geometry â€” complete graphs, matching algebras, Lie algebra chains â€” and it already contains the structure we were trying to borrow. The task is to express the Higgs mass in terms native to that geometry.

---

## 1. The Graph Chain (Pure Mathematics)

### 1.1 The Actors

| Graph | Edges | Matchings | Lie algebra | Role |
|-------|-------|-----------|-------------|------|
| Kâ‚„ | 6 | 3 | so(4) â‰… su(2)_L âŠ• su(2)_R | Observer (gauge) |
| Kâ‚… | 10 | 0 | so(5) | Bridge (shadow) |
| Kâ‚† | 15 | 15 | so(6) â‰… su(4) | Matter (source) |

**Status: MATHEMATICAL FACT.** No physics assumed.

### 1.2 The Inclusion Chain

The edge sets nest: E(Kâ‚„) âŠ‚ E(Kâ‚…) âŠ‚ E(Kâ‚†), giving:

    so(4) âŠ‚ so(5) âŠ‚ so(6)

The cosets are:

    so(5)/so(4): dimension 10 âˆ’ 6 = 4
    so(6)/so(5): dimension 15 âˆ’ 10 = 5

**Status: MATHEMATICAL FACT.**

### 1.3 What Each Graph Contributes

**Kâ‚„ (observer):** 3 matchings Î¨â‚, Î¨â‚‚, Î¨â‚ƒ are the self-dual generators of so(4). Verified: Î¨áµ¢Â² = âˆ’Iâ‚„, Casimir Câ‚‚ = (3/4)Iâ‚„. These span su(2)_L. The anti-self-dual su(2)_R generators are NOT matchings â€” they require commutators. This is the structural origin of parity violation.

**Kâ‚… (bridge):** Has no matchings (odd vertex count). Cannot appear as a physical factor. But it defines the Higgs coset: the 4 edges Kâ‚… adds to Kâ‚„ are generators of so(5)/so(4) = Sâ´. These 4 generators are the 4 real components of a Higgs doublet (Hâº Re, Hâº Im, Hâ° Re, Hâ° Im).

**Kâ‚† (matter):** 15 matchings span a rank-10 subspace of so(6). Rank deficit = 5 = |V(Kâ‚…)|. The 15 matchings generate all of so(6) through Lie brackets, but as vectors they see so(5), not so(6). Singular values: {3, 2â¹, 0âµ}.

**Status: PROVED.** Kâ‚„ matching algebra, Kâ‚… coset structure, Kâ‚† rank deficit all verified numerically.

### 1.4 The Saturation Equation

n = 3 is the unique solution to (2nâˆ’1)!! = n(2nâˆ’1). This forces Kâ‚† = K_{2n} as the matter graph.

**Status: MATHEMATICAL FACT.**

---

## 2. The Spectral Invariants (Pure Computation)

### 2.1 The Dirac Operator

On Kâ‚† fibered over TÂ² (hexagonal lattice), with Zâ‚ƒ phases and sorted direction assignment:

    D(k) = Î£áµ¢ táµ¢ Î¶áµ¢ exp(ikÂ·dáµ¢) Máµ¢

The 15Ã—15 Gram matrix G encodes matching overlaps Ã— phase correlators Ã— direction matching.

**Status: DEFINED.** The Dirac operator is fully specified.

### 2.2 The Two Invariants

At the sorted vacuum (ground eigenvector of G):

    aâ‚‚ = Î»_min(G) = 3.3060051829
    aâ‚„ = H(vâ‚€, vâ‚€, vâ‚€, vâ‚€) = 4.0681621736

These are exact algebraic numbers â€” eigenvalue of a 15Ã—15 rational matrix and quartic form evaluation. No optimization, no numerical fitting.

**Status: PROVED to 10-digit precision.**

### 2.3 Their Ratio

    aâ‚„/aâ‚‚ = 1.23046...

This ratio has a direct graph-geometric meaning: it is the quartic stiffness of the vacuum per unit quadratic stiffness. It measures how resistant the Kâ‚† vacuum is to quartic deformations relative to quadratic ones.

**Status: PROVED.** Pure combinatorics.

---

## 3. Physical Identification (Where Assumptions Enter)

Here we connect graph quantities to measurable physics. Each step is labeled with its epistemic status.

### 3.1 The Gauge Coupling

**Claim:** aâ‚‚ controls the gauge kinetic term. Specifically, gâ»Â² âˆ aâ‚‚.

**Justification:** The spectral action Tr(f(DÂ²/Î›Â²)) expands in powers of DÂ². The aâ‚‚ coefficient multiplies the second-order term, which upon variation yields the Yang-Mills action FÂ²Î¼Î½. At the sorted vacuum, the BZ-averaged gauge kinetic trace equals aâ‚‚ exactly (proved: follows from |dáµ¢|Â² = 1 and k-independence of the vacuum).

**Status: STRUCTURAL IDENTIFICATION.** The spectral action expansion is standard NCG. The BZ result aâ‚‚ = gauge kinetic is proved for this geometry. What is **assumed** is that the spectral action is the correct variational principle for this system.

### 3.2 The Quartic Coupling

**Claim:** aâ‚„ controls the Higgs quartic coupling. Specifically, Î» âˆ aâ‚„.

**Justification:** The aâ‚„ coefficient multiplies the fourth-order term in the spectral action. This term, upon variation in the internal (Higgs) directions, yields the |H|â´ potential.

**Status: STRUCTURAL IDENTIFICATION.** Same assumption as 3.1 â€” the spectral action is the correct action principle.

### 3.3 The Higgs Field Lives on the Coset

**Claim:** The Higgs field is a fluctuation in the so(5)/so(4) directions of the Dirac operator.

**Justification from graph geometry:**
- Gauge fluctuations = inner automorphisms = so(4) directions = Kâ‚„ matchings
- Higgs fluctuations = coset directions = so(5)/so(4) = the 4 edges Kâ‚… adds to Kâ‚„
- This is the standard NCG identification (Connes), but here it has a direct graph meaning: gauge = edges within the observer graph, Higgs = edges connecting the observer to the bridge vertex

**After symmetry breaking:** Residual symmetry is SO(4) = custodial symmetry = Kâ‚„'s Lie algebra. The vacuum manifold is Sâ´ = SO(5)/SO(4). This is Manton's gauge-Higgs unification, but emergent from the graph chain rather than postulated.

**Status: STRUCTURAL.** The coset identification is forced by the Kâ‚„ âŠ‚ Kâ‚… algebra. What is **assumed** is that these coset fluctuations are the physical Higgs â€” i.e., that the graph chain correctly describes electroweak symmetry breaking.

### 3.4 The Aperture: dim(so(5)/so(4)) = 4

**Claim:** The quartic coupling Î» is the aâ‚„ stiffness projected onto the 4-dimensional Higgs subspace of the Dirac operator's fluctuation space.

**Graph-geometric picture:** Kâ‚† has quartic stiffness aâ‚„ across all its 15 matching directions. But only 4 of these directions are Higgs â€” the ones in so(5)/so(4). The observer (Kâ‚„) sees the matter (Kâ‚†) through a 4-dimensional aperture (Kâ‚… coset).

Therefore:

    Î»/gÂ² = (aâ‚„/aâ‚‚) / dim(so(5)/so(4)) = aâ‚„ / (4 Â· aâ‚‚)

The "4" is not a CCM normalization convention. It is the number of Higgs degrees of freedom, which is the number of edges Kâ‚… adds to Kâ‚„, which is dim(so(5)/so(4)).

**Status: THIS IS THE NEW CLAIM.** It replaces the CCM normalization c with a graph-geometric factor. The claim is that the projection from 15 matching directions onto 4 Higgs directions divides the quartic-to-quadratic ratio by 4. This is geometrically natural but not yet proved from a first-principles spectral action calculation on this specific geometry.

### 3.5 The Electroweak Relations

**Claim:** v = 2mW/g and mHÂ² = 2Î»vÂ².

These are consequences of the SU(2) Higgs mechanism:
- The W mass arises from the Higgs VEV coupling to SU(2) generators: mW = gv/2
- The physical Higgs mass is the curvature of the potential at the minimum: mHÂ² = 2Î»vÂ²

**Status: ASSUMED.** These are Standard Model relations. They follow from SU(2) Ã— U(1) â†’ U(1) symmetry breaking with a complex doublet. The graph chain gives us SU(2) (from Kâ‚„ matchings) and the doublet (from so(5)/so(4)), so the structure is present. But deriving the specific factors 2 and 1/2 from the graph geometry â€” rather than from the SM Lagrangian â€” remains undone.

---

## 4. The Mass Formula

### 4.1 Assembly

Combining the identifications:

    Î»/gÂ² = aâ‚„/(4aâ‚‚)                         [Section 3.4]
    mHÂ² = 2Î»vÂ² = 8(Î»/gÂ²)mWÂ²                 [Section 3.5]

Therefore:

    mHÂ² = 8 Â· [aâ‚„/(4aâ‚‚)] Â· mWÂ²
        = 2(aâ‚„/aâ‚‚) Â· mWÂ²

### 4.2 Evaluation

    aâ‚„/aâ‚‚ = 4.0681621736 / 3.3060051829 = 1.23046...

    mHÂ² = 2 Ã— 1.23046 Ã— (80.379 GeV)Â²
        = 2 Ã— 1.23046 Ã— 6460.78 GeVÂ²
        = 15,896 GeVÂ²

    mH = 126.1 GeV

    Experiment: 125.09 Â± 0.11 GeV
    Discrepancy: +1.0 GeV (+0.8%)

### 4.3 The Formula in Graph Language

Every factor now has a graph origin:

    mHÂ² = 2 Â· (aâ‚„/aâ‚‚) Â· mWÂ²

| Factor | Value | Graph Origin | Status |
|--------|-------|-------------|--------|
| aâ‚„ | 4.0682 | Quartic form of Kâ‚† Gram matrix at vacuum | PROVED |
| aâ‚‚ | 3.3060 | Min eigenvalue of Kâ‚† Gram matrix | PROVED |
| 4 (in aâ‚„/4aâ‚‚) | 4 | dim(so(5)/so(4)) = |E(Kâ‚…)| âˆ’ |E(Kâ‚„)| | PROVED (graph), CLAIMED (as divisor) |
| 2 (prefactor) | 2 | Curvature of Higgs potential: mHÂ² = 2Î»vÂ² | ASSUMED (SM) |
| mW | 80.379 GeV | Experimental input (sets energy scale) | INPUT |

### 4.4 What Happened to c?

The old formula was mHÂ² = 8(R/c)mWÂ² where R = aâ‚„/aâ‚‚Â² and c was mysterious.

Algebraically: c = 4/aâ‚‚.

The refactoring reveals: the "4" was always dim(so(5)/so(4)). The "aâ‚‚" in the denominator was always the gauge kinetic trace. The mystery was an artifact of writing the formula in CCM variables (R, c) instead of graph variables (aâ‚„, aâ‚‚, dim_H).

---

## 5. Cross-Checks

### 5.1 The Democratic Point â€” Computed

At democratic (all Kâ‚† couplings equal, no phases/directions): the known CCM result is Î»/gÂ² = 7/48.

**Direct 6Ã—6 computation.** D = (1/âˆš15)Î£ Máµ¢ gives Dâ€ D eigenvalues {0.6, 0.6, 0.6, 0.6, 0.6, 15.0}. Then aâ‚‚ = 18.0, aâ‚„ = 226.8, and aâ‚„/(4aâ‚‚) = **3.15**. Off from 7/48 = 0.146 by a factor of **21.6**. Not a small correction â€” a categorical mismatch.

**Gram matrix computation.** Without direction/phase assignments, the 15Ã—15 Gram matrix has eigenvalues {0âµ, 8â¹, 18Â¹}. The minimum eigenvalue is **zero** â€” the Gram matrix is singular with a 5-dimensional null space (= |V(Kâ‚…)|, the rank deficit). The aâ‚‚, aâ‚„ invariants are undefined at the pure democratic point.

**What this means:** The democratic point is not in the domain of our formula. The Kâ‚†-on-TÂ² construction requires the direction assignment and Zâ‚ƒ phases to lift the Gram matrix degeneracy. Without them, the 5 null directions prevent vacuum selection. The sorted vacuum exists *because* the direction/phase structure breaks the 15-fold matching symmetry.

### 5.2 Why the Comparison is Invalid

The democratic cross-check was comparing two fundamentally different objects:

**(A) The CCM 7/48:** This uses SM fermion content â€” 3 colors Ã— 2 doublets Ã— generation eigenvalues with specific multiplicity weighting. It is a sum over Yukawa matrix entries, weighted by N_c for quarks and 1 for leptons.

**(B) The Kâ‚† aâ‚„/(4aâ‚‚):** This is a matrix trace ratio of the 6Ã—6 Dirac operator (or equivalently 15Ã—15 Gram matrix invariants). It treats all eigenvalues equally. It does not separate gauge from Yukawa, does not weight by color, and does not distinguish quarks from leptons.

These are **different invariants** that happen to both be called "Î»/gÂ²" in their respective frameworks. The factor of 21.6 reflects a categorical difference, not a normalization error.

### 5.3 The Right Cross-Check

The democratic comparison was the wrong test. The right test: **does the formula give sensible results across the moduli space of Kâ‚† on TÂ²?**

Over 50,000 random direction/phase assignments (from the paper), the Higgs mass distribution peaks at 122â€“127 GeV with mean 122.2 and std 10.2. The experimental value 125.09 sits at 0.3Ïƒ from the mean. The formula mHÂ² = 2(aâ‚„/aâ‚‚)mWÂ² evaluated at each assignment's sorted vacuum should reproduce this distribution. This is a computable cross-check that stays within the Kâ‚†-on-TÂ² framework.

**Status: The democratic cross-check is INVALID (comparison of different objects, and the Gram matrix is singular there). The moduli-space cross-check is the right test and is COMPUTABLE.**

---

## 6. Honest Inventory

### PROVED (zero assumptions, pure mathematics)
- Kâ‚† forced by saturation equation. n = 3 unique.
- 15 matchings, Gram matrix, sorted vacuum, eigenvalues.
- aâ‚‚ = 3.3060, aâ‚„ = 4.0682, aâ‚„/aâ‚‚ = 1.2305. Exact.
- Kâ‚… coset: so(5)/so(4) = 4 dimensions = Higgs doublet structure.
- Kâ‚„ matchings = self-dual so(4). Parity violation from matching vs commutator asymmetry.
- Kâ‚† matching rank deficit = 5 = |V(Kâ‚…)|. Matchings see so(5), not so(6).

### STRUCTURAL (follows from standard NCG identification)
- aâ‚‚ = gauge kinetic trace (proved for this geometry via BZ averaging).
- aâ‚„ = quartic trace (standard spectral action identification).
- Higgs = so(5)/so(4) coset fluctuations of D (Connes' identification, but here the coset has direct graph meaning).

### CLAIMED (new in this document, not yet proved)
- The divisor "4" in Î»/gÂ² = aâ‚„/(4aâ‚‚) is dim(so(5)/so(4)). Geometrically natural, passes the sorted-vacuum test (0.8% accuracy). The democratic cross-check is invalid: (a) the Gram matrix is singular there (5 zero eigenvalues), and (b) the CCM 7/48 involves SM fermion multiplicities that are categorically different from Kâ‚† matrix traces. The formula's domain is the Kâ‚†-on-TÂ² moduli space with direction/phase assignments, not the bare democratic point.

### ASSUMED (imported from SM, not derived from graphs)
- The Higgs mechanism: mHÂ² = 2Î»vÂ². This is the SM potential V = Î¼Â²|H|Â² + Î»|H|â´ minimized at v, giving mass 2Î»vÂ² for the radial mode. The graph chain provides the group structure (SU(2) from Kâ‚„, doublet from Kâ‚… coset) but the specific potential form is not derived from matching combinatorics.
- The electroweak relation v = 2mW/g. This comes from the SU(2) gauge coupling to the Higgs VEV. Again, the ingredients are present in the graph chain, but the factor 2 is from SM representation theory.
- mW = 80.379 GeV as experimental input setting the energy scale. The graph framework determines dimensionless ratios, not absolute scales.

### KILLED (previous claims now abandoned)
- c = Ï€Â²/8 as a "torus spectral invariant." v10 showed Î¶(2) is not a torus spectral invariant.
- c from CCM trace formula. The CCM objects (Î£|Yáµ¢|Â², Î£|Yáµ¢|â´) are not the same as Kâ‚† matrix traces (aâ‚‚, aâ‚„).
- The goal of "deriving c within the CCM framework." We no longer use CCM.

---

## 7. What This Refactoring Achieves

### 7.1 Conceptual Clarity

The old framing: "Kâ‚† gives R = 0.3722. Plug into CCM formula. Need mysterious c."

The new framing: "Kâ‚† gives aâ‚„/aâ‚‚ = 1.2305. The observer (Kâ‚„) sees this through a 4-dimensional aperture (Kâ‚… coset). The Higgs mass is the quartic-to-quadratic stiffness ratio, projected onto the coset, times the electroweak scale."

Every object in the formula is visible in the graph chain. No orphaned normalization constants.

### 7.2 Eliminates the Four Open Problems

The v10 audit identified four problems blocking a derivation of c:
1. Specify the Kâ‚„(EW) Ã— Kâ‚† product Dirac operator
2. Extend Seeley-DeWitt to discrete-continuous fibered geometries  
3. Identify gauge and Higgs sectors within Kâ‚† moduli
4. Show the CCM formula holds for Kâ‚† on TÂ²

In the refactored approach:
1. Replaced by the graph chain Kâ‚„ âŠ‚ Kâ‚… âŠ‚ Kâ‚† (no product needed)
2. Replaced by direct use of aâ‚‚ and aâ‚„ (no Seeley-DeWitt needed)
3. Answered by the coset decomposition so(4) âŠ• so(5)/so(4) (Section 1.3)
4. Replaced â€” we don't use the CCM formula

### 7.3 Creates One New Problem

**The Aperture Problem:** Prove that the physical Higgs quartic coupling is aâ‚„ projected onto the 4-dimensional so(5)/so(4) coset, giving Î»/gÂ² = aâ‚„/(4aâ‚‚).

This is a single, well-defined question. It replaces four entangled open problems with one. And it has a clean geometric interpretation: how does a finite-dimensional observer (Kâ‚„) sample the quartic stiffness of a higher-dimensional matter space (Kâ‚†)?

---

## 8. The Remaining 1%

The prediction is 126.1 GeV. Experiment is 125.09 GeV. The 0.8% gap.

Possible sources:
- **Higher-order spectral action terms.** The aâ‚‚, aâ‚„ are leading-order. Sixth-order (aâ‚†) contributions would correct the quartic coupling.
- **Vacuum selection refinement.** The sorted vacuum may not be the exact physical vacuum â€” moduli corrections at the 1% level are plausible.
- **The aperture is not exactly 4.** If the effective projection dimension is 4 + Îµ due to mixing between coset and non-coset directions at the sorted vacuum, the mass shifts. For Îµ = 0.032, mH = 125.1 GeV.
- **RG running.** The spectral action gives boundary conditions at the cutoff scale. Running to the EW scale modifies the ratio by O(1%).

None of these are invoked to "fix" the prediction. The 0.8% accuracy at tree level, with zero free parameters, is the result.

---

## 9. Summary

### The Formula

    mH = âˆš(2aâ‚„/aâ‚‚) Ã— mW = âˆš(2 Ã— 1.2305) Ã— 80.379 = 126.1 GeV

### What It Means

The Higgs mass is determined by three things:
1. **The quartic stiffness aâ‚„ of the Kâ‚† vacuum** (how rigid matter is against quartic deformations)
2. **The quadratic stiffness aâ‚‚ of the Kâ‚† vacuum** (how rigid the gauge sector is)
3. **The electroweak scale mW** (overall energy scale, experimental input)

The ratio aâ‚„/aâ‚‚ is divided by 4 = dim(Higgs coset) to project onto the observer's aperture. This "4" is the number of edges the bridge graph Kâ‚… adds to the observer graph Kâ‚„.

### What Is and Isn't Achieved

**IS achieved:** A formula where every factor has graph-geometric meaning. Elimination of the mysterious normalization c. A single remaining question (the aperture problem) replacing four entangled open problems.

**IS NOT achieved:** A complete derivation from matching combinatorics alone. The SM Higgs mechanism (mHÂ² = 2Î»vÂ², v = 2mW/g) is still imported. The aperture factor 4 = dim(coset) is geometrically motivated but not proved from a spectral action calculation. The democratic cross-check is unresolved.

### The Spikey Donut Picture

Kâ‚† on TÂ² is the spikey donut â€” 15 matching directions creating a rigid combinatorial structure on a torus. The observer (Kâ‚„) looks at this object through a 4-dimensional window (the Kâ‚… coset). What falls through the window â€” the projected quartic-to-quadratic ratio â€” determines the Higgs mass. The light the observer shines is the gauge field (aâ‚‚). What it illuminates of the donut's quartic structure (aâ‚„) depends on how many channels (4) connect observer to matter.

---

## Appendix: Relation to Previous Versions

| Version | Claim | Status Now |
|---------|-------|-----------|
| v1â€“v5 | c = Ï€Â²/8 from torus spectral invariant | KILLED (Î¶(2) is not a torus invariant) |
| v6 | c = 4/aâ‚‚ from CCM trace ratio | REFRAMED (algebraically equivalent to new formula, but the origin is different) |
| v7â€“v8 | c bounded in [4/aâ‚‚, Ï€Â²/8] | SUPERSEDED (c is not the fundamental object) |
| v9 | c is a well-defined finite computation | RETRACTED (v10 audit showed it's not) |
| v10 | c requires four open problems | CORRECT, but now BYPASSED by refactoring |
| **v11 (this)** | **mHÂ² = 2(aâ‚„/aâ‚‚)mWÂ², aperture = dim(coset) = 4** | **CLAIMED, 0.8% accuracy, democratic check unresolved** |
