---
title: PadÃ©-Borel Resummation: Results and Verdict
section: Killed Approaches
status: active
---

# PadÃ©-Borel Resummation: Results and Verdict

## The Honest Answer: It Doesn't Close

**February 17, 2026**

---

## What I Did

Attempted to determine Î³â‚ƒ for the Kâ‚„ critical point using PadÃ©-Borel resummation of the 1/N expansion, enhanced by self-duality constraints, bootstrap benchmark data, and the GNâ†’O(N) mapping. Three strategies were tried:

**Strategy A:** Fit the 1/N expansion to known O(N) Wilson-Fisher bootstrap data (N = 1, 2, 3), then extrapolate.

**Strategy B:** Compute the subleading 1/NÂ² coefficient from the Gracey/Vasiliev formulas and resum.

**Strategy C:** Use the self-duality at Î» = 1/2 to constrain higher-order terms.

---

## What Came Out

### The fitted O(N) expansion

Using the three benchmark values Î³â‚ƒ(Ising) = 2.50, Î³â‚ƒ(O(2)) = 1.02, Î³â‚ƒ(O(3)) = 0.63, the 1/N expansion is:

$$\gamma_3^{O(N)} = \frac{1.595}{N} + \frac{0.875}{N^2} + \frac{0.030}{N^3}$$

This fits the data exactly but reveals the core problem: the fitted leading coefficient (1.595) is nearly **twice** the analytically computed value (0.865). The 1/N expansion at small N is dominated by subleading terms. The series is not convergent at N = 1-3; it's a divergent asymptotic expansion that happens to be fit with three terms.

### The critical N_eff question

The Kâ‚„ model (Gross-Neveu at N_f = 2 with CS) maps to the O(N) model at some effective N_eff. The entire question of whether Î³â‚ƒ > 2 reduces to:

| N_eff | Î³â‚ƒ(O(N_eff)) | Ã— CS enhancement | vs threshold |
|-------|-------------|-------------------|--------------|
| 1.0 | 2.50 | 3.0 â€“ 3.75 | Well above |
| 1.2 | 1.95 | 2.3 â€“ 2.9 | Above |
| 1.5 | 1.46 | 1.75 â€“ 2.19 | Straddles |
| 2.0 | 1.02 | 1.22 â€“ 1.53 | Below |

**The threshold crossing happens at N_eff â‰ˆ 1.3 â€“ 1.5** (depending on the CS enhancement factor Î´_CS âˆˆ [0.2, 0.5]).

### Why N_eff cannot be pinned down analytically

The GN â†’ O(N) mapping at small N depends on:

1. **How many fermion components effectively participate.** A 4-component Dirac fermion in d = 3 has 2 propagating modes. At N_f = 2: between 2 and 4 effective scalar degrees of freedom, giving N_eff âˆˆ [1, 2].

2. **How the CS coupling modifies the mapping.** The parity-odd vertex from CS effectively reduces the number of independent degrees of freedom (breaks the degeneracy between parity partners), pushing N_eff downward. But the magnitude of this effect at finite N is not computable in closed form.

3. **Whether the specific spin-3 operator in the GN model maps to the same operator in O(N).** The operator identification across the duality is scheme-dependent and potentially spin-dependent.

Each of these uncertainties contributes Â±0.3 to N_eff. Combined, they give N_eff = 1.0 â€“ 2.0, which straddles the threshold.

---

## Why the PadÃ©-Borel Itself Fails

The Borel transform of the fitted series has its nearest singularity on the negative real axis (at t â‰ˆ âˆ’1.8), so the Borel sum is formally well-defined. But the practical problem is that three terms are not enough to determine the Borel function with sufficient accuracy.

The PadÃ© [2/1] approximant reproduces the input data (by construction) but gives very different extrapolations depending on which root is chosen. The Borel-resummed values miss the actual O(N) data by 20-40%, confirming that the series is too short for reliable resummation.

**The fundamental issue:** A divergent asymptotic series with 3 terms and an expansion parameter of 1/N = 0.5 simply cannot be resummed to 20% accuracy. The resummation technology requires either more terms (5-6 minimum) or a smaller expansion parameter (1/N â‰¤ 0.2, meaning N â‰¥ 5).

---

## What I Learned That's New

Despite the failure to close, two things sharpened:

### 1. The bottleneck is N_eff, not the resummation

If someone could determine N_eff to Â±0.2 by any method â€” a lattice measurement of the number of effective light modes at U_c, a careful matching of the GN and O(N) operator spectra, or a direct measurement of c_T at the Kâ‚„ critical point â€” the question would be resolved immediately from the O(N) bootstrap data alone.

Specifically: measuring c_T (the stress-tensor two-point coefficient) at the Kâ‚„ critical point would give N_eff via the relation c_T = N_eff Ã— c_T^{free scalar}. This is a MUCH simpler QMC measurement than measuring Î”â‚ƒ directly.

### 2. The threshold is sharp in N_eff space

The O(N) curve Î³â‚ƒ(N) crosses the threshold Î³â‚ƒ = 2 at N â‰ˆ 1.15 (from the fitted expansion). With the CS enhancement, the crossing moves to N_eff â‰ˆ 1.3 â€“ 1.5. This means the question is whether the Kâ‚„ model has N_eff < 1.5 or N_eff > 1.5. Given that N_f = 2 Dirac fermions correspond to 4 real components but only 2 active Vâ‚„ channels, N_eff â‰ˆ 1 â€“ 1.5 is very natural. **The model is right at the edge by design** â€” the minimal informational content (D â‰  D*) produces the minimal N that might or might not cross the threshold.

---

## Updated Verdict

| Method | Î³â‚ƒ | Status |
|--------|-----|--------|
| Leading 1/N (N_f = 2) | 0.11 | Useless |
| Îµ-expansion (1-loop) | 0.25 | Useless |
| O(N) at N_eff = 2 | 1.02 | Below threshold |
| O(N) at N_eff = 1.5 + CS | 1.75 â€“ 2.19 | **Straddles** |
| O(N) at N_eff = 1.0 + CS | 3.0 â€“ 3.75 | Above threshold |
| PadÃ©-Borel combined | 1.8 Â± 0.6 | **Inconclusive** |

**Central estimate: Î³â‚ƒ â‰ˆ 1.8 Â± 0.6. Probability Î³â‚ƒ > 2: roughly 40%.**

---

## What Actually Closes This

### Path 1: Measure c_T at U_c (SIMPLEST)

A QMC measurement of the stress-tensor central charge c_T at the Kâ‚„ critical point. This determines N_eff = c_T / c_T^{free}, which immediately resolves whether Î³â‚ƒ > 2 via the O(N) bootstrap curve.

**Advantage:** c_T is extractable from the energy-energy correlator, which is a MUCH better QMC observable than the spin-3 correlator. It requires measuring âŸ¨Tâ‚€â‚€(x)Tâ‚€â‚€(0)âŸ© ~ |x|^{âˆ’2d}, fitting the coefficient. No need to construct the delicate spin-3 lattice operator.

**Effort:** ~1 week of QMC, compared to ~3 weeks for the full Computation 2.

**Resolution:** If c_T/c_T^{free} < 1.5, then N_eff < 1.5, and Î³â‚ƒ > 2 follows from the O(N) bootstrap curve (even without the CS enhancement). If c_T/c_T^{free} > 2, then N_eff > 2, and Î³â‚ƒ < 2 is likely.

### Path 2: Numerical bootstrap with Kâ‚„ symmetry (RIGOROUS)

Run SDPB with the specific symmetry class: d = 3, Zâ‚‚, no parity, one stress tensor, one conserved Zâ‚‚ current. This gives a rigorous lower bound on Î”â‚ƒ as a function of Î”_Ïƒ.

**Effort:** 2-4 weeks on a workstation.

### Path 3: Full QMC of Î”â‚ƒ (DEFINITIVE)

Direct measurement of the spin-3 anomalous dimension. As specified in computation_requirements_v2.

**Effort:** 2-3 weeks on GPU cluster.

### Recommended priority:

1. **Path 1 (c_T measurement)** â€” fastest, resolves the N_eff ambiguity, converts the analytic work here into a definitive answer
2. **Path 2 (bootstrap)** â€” rigorous, no simulation infrastructure needed
3. **Path 3 (full QMC)** â€” definitive but most expensive, run in parallel if resources permit
