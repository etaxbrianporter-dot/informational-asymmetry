---
title: Kâ‚„D Theory Part 2: Cosmological Cycling from Boundary Entropy Saturation
section: K₄ Spacetime
status: active
---

# Kâ‚„D Theory Part 2: Cosmological Cycling from Boundary Entropy Saturation

## Working Paper â€” February 2026

---

## 0. Relationship to Part 1

Part 1 derived spacetime (gravity sector) and the Standard Model (gauge sector) from matching combinatorics on Kâ‚„ and Kâ‚†. The cosmological implications were flagged but not developed:

- **Cosmological constant:** d_eff = 3 âˆ’ f(Îµ), status REDUCED. The functional form of f(Îµ) at large Îµ was identified as an open computation.
- **Black hole information paradox:** Status ELIMINATED. "Information cannot be lost inside a bulk that is emergent. Fundamental degrees of freedom live on the 2+1 boundary."
- **Dark matter/dark energy:** Status SILENT. "No mechanism. The framework determines algebraic structure, not cosmological content."

Part 2 takes the structural results of Part 1 as input and asks: **what does the Kâ‚„D boundary theory say about cosmological evolution, entropy saturation, and cycling?**

The guiding principle from Part 1: be honest about status labels. Every claim gets one of: FORCED (follows from proved structure with no additional assumptions), DERIVED (follows from Part 1 results plus a stated assumption), CONJECTURED (plausible but unproved), or OPEN (well-posed question, answer unknown).

---

## 1. What Part 1 Forces About Cosmology

### 1.1 The boundary is fundamental, the bulk is emergent

**Status: FORCED (from Part 1)**

The Kâ‚„ lattice model on a triangular lattice is 2+1-dimensional. The emergent 3+1 spacetime arises when a finite observer (Axiom 4) cannot resolve the spectral structure, and d_eff â†’ 3 at coarse resolution. The bulk is what the boundary looks like to a finite observer. This is not a holographic conjecture â€” it is the construction.

**Implication for cosmology:** Any cosmological question must have an answer stated in boundary language. If it can only be stated in bulk language, the question is ill-posed. "What happens inside a black hole?" is ill-posed. "What is the boundary entanglement configuration corresponding to a black hole?" is well-posed.

### 1.2 Asymmetry is topological and permanent

**Status: FORCED (from Axiom 3, D â‰  D*)**

The â„¤â‚ƒ flux that breaks time-reversal is topological: C = âˆ’2 is a Chern number, an integer invariant. It cannot be continuously deformed to zero. The arrow of time is not emergent from thermodynamics â€” it is algebraically prior. The â„¤â‚ƒ flux exists at every point of the Kâ‚„ boundary regardless of the boundary's entanglement state.

**Implication for cosmology:** No cosmological process â€” not heat death, not entropy saturation, not any phase transition â€” can eliminate the arrow of time. Any cycling mechanism must preserve D â‰  D*. The â„¤â‚ƒ flux survives whatever transition occurs between aeons.

### 1.3 Entanglement entropy is bounded

**Status: DERIVED (from Part 1 + lattice finiteness)**

The boundary is a Kâ‚„ lattice on a triangular torus of size L Ã— L. Each site has 4 orbitals. The maximum entanglement entropy of a boundary region of linear size â„“ is bounded by â„“ Â· log(4) (one bit per orbital per boundary site, area-law saturation). For the full boundary: S_max âˆ LÂ² Â· log(4).

The computed entanglement entropy at the free point is S = (2/3)LÂ·log(L/Îµ), which grows logarithmically â€” far below the area-law maximum. The system has enormous room between its current entanglement and saturation.

**Implication for cosmology:** There exists a maximum entropy state of the boundary. The bulk observer perceives this as a maximum entropy state of the universe. The question is: what does the emergent bulk look like as the boundary approaches this maximum?

### 1.4 The cosmological constant tracks observer resolution

**Status: REDUCED (from Part 1, open computation)**

d_eff = 3 âˆ’ f(Îµ) where f(Îµ) â†’ 0 as Îµ â†’ âˆž. If f(Îµ) ~ 1/ÎµÂ², then Î› ~ (L_Planck/L_Hubble)Â² ~ 10â»Â¹Â²Â². Sign is positive (d_eff < 3).

**Implication for cosmology:** Î› is not a fixed constant of nature. It is a property of the observer-boundary relationship. As the boundary's entanglement state evolves, the effective resolution Îµ changes, and Î› changes with it.

---

## 2. The Derivation Chain

### Step A: Entropy evolution on the boundary

**Status: CONJECTURED â€” needs rigorous derivation from Kâ‚„ dynamics**

**Claim A1:** Under unitary Kâ‚„ boundary dynamics (the Hamiltonian H = Hâ‚€ + H_int from the computation requirements document), the total boundary entanglement entropy increases monotonically for a coarse-grained observer.

**Basis:** This is standard in quantum statistical mechanics â€” entanglement entropy of a subsystem increases under generic unitary evolution when coarse-grained. The Kâ‚„ boundary is a concrete quantum lattice model. The second law should follow from the same mechanism as in any quantum lattice system.

**What needs to be checked:** That the Vâ‚„ Ã— â„¤â‚ƒ symmetry structure doesn't create protected sectors that prevent thermalization. The â„¤â‚ƒ flux breaks time-reversal, which actually *helps* thermalization (prevents many-body localization from time-reversal-protected resonances). The Vâ‚„ symmetry is small (4 elements) and unlikely to create large protected subspaces. But this needs to be verified, not assumed.

**Computation A1:** Exact diagonalization of H at small L (L = 3, 6) to check that the entanglement entropy of half the boundary grows monotonically from a generic initial state. Estimate thermalization time as a function of U/t.

### Step B: Maximum entropy configuration of the boundary

**Status: CONJECTURED â€” definition clear, properties unknown**

**Claim B1:** The maximum entropy configuration of the Kâ‚„ boundary is a state where the two-point correlation function âŸ¨câ€ _i c_jâŸ© decays exponentially at the shortest possible scale (~ 1 lattice spacing) in every Vâ‚„ channel.

**Claim B2:** In this configuration, the Fisher information metric on (kâ‚, kâ‚‚, E) becomes trivial: g_ij â†’ Î´_ij with no structure. A finite observer cannot distinguish any direction or location. The emergent bulk is maximally homogeneous and isotropic.

**Basis for B2:** The Fisher metric is computed from the spectral function of D(k). If all correlations are short-ranged, the spectral function is featureless, and the Fisher metric carries no information about spatial structure. This is what de Sitter space looks like to a finite observer: featureless expansion with no structure to anchor to.

**What needs to be checked:** Explicit computation of the Fisher metric at the maximally mixed state of the boundary. Does it go exactly trivial, or does some residual structure from the â„¤â‚ƒ flux survive?

**Computation B1:** Compute g_ij for the Kâ‚„ lattice at infinite temperature (Î² â†’ 0). The â„¤â‚ƒ flux should still be present (it's topological), so the Fisher metric may retain a vestigial anisotropy even at maximum entropy. If it does, this vestigial structure constrains what can carry over between aeons.

### Step C: What the bulk observer sees at maximum boundary entropy

**Status: DERIVED (from Part 1 machinery + Steps A, B)**

If Step B holds, then at maximum boundary entropy:

- d_eff â†’ 3 with f(Îµ) â†’ 0 (the observer sees pure 3D space, no energy-space mixing)
- Î› â†’ 0 (the cosmological constant vanishes as the boundary homogenizes)
- The emergent geometry is approximately de Sitter with decreasing Î›
- All massive particles have decayed (their boundary entanglement patterns have dissolved into the featureless background)
- Only the topological invariants survive: â„¤â‚ƒ flux, Chern number C = âˆ’2, Vâ‚„ algebra

**This is Penrose's precondition for conformal recycling, derived from boundary physics.**

In Penrose's CCC, the end of an aeon is when only conformally invariant (massless) physics remains and the infinite future becomes conformally equivalent to a big bang. In Kâ‚„D, this translates to: the boundary reaches maximum entanglement, all bulk structure dissolves, but the topological invariants (which cannot dissolve because they're discrete/integer-valued) persist and seed the next epoch.

### Step D: The transition mechanism

**Status: CONJECTURED â€” this is the hardest and most speculative step**

**Question D1:** Is the maximum entropy state of the Kâ‚„ boundary *stable* or *unstable*?

Three possibilities:

**D1a: Stable equilibrium.** The boundary reaches maximum entropy and stays there. Heat death is permanent. No cycling. This is the conservative answer. It would mean Kâ‚„D is compatible with standard cosmology but adds nothing beyond it for late-time evolution.

**D1b: Unstable via boundary phase transition.** The Kâ‚„ boundary at maximum entanglement undergoes a phase transition analogous to deconfinement. At high entanglement density, the effective coupling changes character, and the boundary reorganizes into a new low-entanglement configuration. The bulk observer perceives this as a new big bang.

**Basis for D1b:** The Kâ‚„ model has a known phase structure â€” the semimetal-insulator transition at U = U_c. If the boundary dynamics drives the *effective* coupling through a phase boundary as entanglement accumulates, a first-order transition could catastrophically reorganize the boundary. This is speculative but concrete: the boundary theory is a well-defined quantum lattice model, and its phase diagram is computable.

**D1c: Unstable via topological mechanism.** The â„¤â‚ƒ flux, which is topological and permanent, acts as a "seed" that prevents true equilibrium. Because D â‰  D*, the boundary dynamics is inherently irreversible (Axiom 3). A truly reversible equilibrium is algebraically forbidden â€” the system must always be producing entropy. When it runs out of capacity to absorb new entropy (saturation), it must reorganize.

**Basis for D1c:** This is the most Kâ‚„D-native argument. If D â‰  D* means the system can never reach detailed balance (because time-reversal is broken at the fundamental level), then the maximum entropy state is not a true equilibrium â€” it's a dynamical state that is still evolving, just with nowhere to put its entropy. This creates a pressure for reorganization.

**What would distinguish D1a from D1b/D1c:** The answer depends on whether the Kâ‚„ boundary theory has a unique thermal equilibrium at Î² = 0 (D1a) or whether the topological flux prevents true equilibration (D1b/D1c). This is a well-posed question about a concrete lattice model.

**Computation D1:** Study the Kâ‚„ lattice model at Î² â†’ 0 (high temperature / late universe). Does the â„¤â‚ƒ flux create a residual current or Berry phase that prevents the state from being a true Gibbs equilibrium? If âŸ¨JâŸ© â‰  0 at Î² = 0 (a persistent current from the topological flux), then D1a is eliminated.

### Step E: Memory across aeons

**Status: CONJECTURED â€” depends on Steps C, D**

If cycling occurs (D1b or D1c), the question is: what information survives from one aeon to the next?

**Claim E1: FORCED quantities survive.** By the classification theorem (Part 1), graph-topological invariants are independent of the moduli point and the entanglement state. They are properties of the Kâ‚„ and Kâ‚† graphs themselves. These survive any boundary reorganization:
- Gauge group SU(3) Ã— SU(2) Ã— U(1)
- Three generations
- Parity violation
- Spacetime dimension 3+1
- Chern number C = âˆ’2
- â„¤â‚ƒ flux

**Claim E2: SILENT quantities reset.** Moduli-geometric quantities (masses, couplings, mixing angles) depend on the specific moduli point, which corresponds to the specific entanglement pattern. A boundary reorganization scrambles the entanglement pattern and selects a new moduli point. The laws of physics have the same structure, but different parameters.

**Claim E3: REDUCED quantities are the interesting case.** Moduli-topological quantities (Chern numbers of specific bands, mass hierarchy pattern) take discrete values that depend on which chamber of moduli space the system lands in. A reorganization could land in the same chamber or a different one. The number of distinct aeon-types is the number of chambers in M(Kâ‚„) Ã— M(Kâ‚†).

**Implication:** Each aeon has the *same laws* (same gauge group, same number of generations, same chirality) but potentially *different parameters* (different particle masses, different coupling constants, different cosmological constant). The structural/parametric divide from Part 1's classification theorem becomes the aeon-invariant/aeon-variable divide in Part 2.

---

## 3. Comparison with Penrose CCC

| Feature | Penrose CCC | Kâ‚„D Cycling |
|---------|-------------|-------------|
| **End of aeon** | All massive particles decay; only conformal physics remains | Boundary reaches max entanglement; Fisher metric becomes trivial |
| **Transition mechanism** | Conformal rescaling (mathematical, no dynamics) | Boundary phase transition or topological instability (physical) |
| **What carries over** | Conformal structure, Hawking points | FORCED invariants (gauge group, generations, chirality, C = âˆ’2) |
| **What resets** | Masses, scales | SILENT quantities (masses, couplings, mixing angles) |
| **Observational signature** | Circular temperature anomalies in CMB from prior-aeon BHs | Vâ‚„-structured large-angle CMB anomalies from boundary reorganization pattern |
| **Number of aeons** | Infinite (assumed) | Determined by boundary dynamics (D1a: 1, D1b/D1c: potentially infinite) |
| **Arrow of time** | Resets via conformal rescaling of Weyl tensor | Never resets â€” â„¤â‚ƒ flux is permanent (Axiom 3 is prior to dynamics) |
| **Substrate** | Continuous spacetime | Kâ‚„ lattice on triangular torus |

**Key difference on time's arrow:** In Penrose CCC, the Weyl curvature hypothesis (low initial Weyl curvature at each big bang) resets the gravitational entropy. In Kâ‚„D, the arrow of time (D â‰  D*) never resets because it's topological. What resets is the *entanglement pattern*, not the directionality. This is a stronger position: you don't need an additional hypothesis about why entropy is low at the start of each aeon. The entropy is low because the boundary just reorganized â€” its entanglement was catastrophically rearranged, and the new pattern hasn't yet thermalized.

---

## 4. The Black Hole Role (Revised)

### 4.1 What a black hole IS in Kâ‚„D

**Status: DERIVED (from Part 1 boundary construction)**

A black hole is a region of the 2+1 boundary where entanglement has reached local saturation. The boundary degrees of freedom in that region are maximally entangled with each other, leaving no entanglement capacity for long-range correlations. The bulk observer perceives this as a region from which no information escapes â€” not because it's trapped behind a horizon, but because maximally entangled boundary degrees of freedom look featureless to a finite observer.

The Bekenstein-Hawking entropy S_BH = A/(4l_pÂ²) is the entanglement entropy of the saturated boundary region. The area law follows automatically because boundary entanglement entropy scales with the boundary of the region (perimeter in 2+1), which maps to area in the emergent 3+1 bulk.

### 4.2 Hawking evaporation as boundary disentanglement

**Status: DERIVED**

Hawking radiation is the boundary disentangling. The maximally entangled region is not in a true equilibrium (because D â‰  D* prevents detailed balance even locally). The entanglement slowly leaks from the saturated region into neighboring unsaturated boundary, and the bulk observer perceives this as thermal radiation. The temperature is set by the entanglement gradient at the boundary of the saturated region.

### 4.3 Black holes as cosmological entropy accelerators

**Status: CONJECTURED**

Black holes are local processes that accelerate the boundary's approach to global maximum entropy. Each formation-and-evaporation cycle takes structured (low-entropy) boundary entanglement and redistributes it toward the featureless (high-entropy) maximum. They are not "recycle bins" in the sense of converting matter into new universes. They are more like **mixers** â€” they homogenize the boundary entanglement distribution.

In the cosmological context: a universe with many black holes forming and evaporating will approach boundary entropy saturation faster than one without. The rate of black hole formation controls the timescale of the current aeon.

### 4.4 No interior, no information paradox, no firewall

**Status: FORCED (from Part 1)**

There is no black hole interior in Kâ‚„D. The "interior" is an artifact of taking the emergent 3+1 bulk as fundamental. The fundamental degrees of freedom are on the 2+1 boundary at all times. There is nothing for information to be "lost" into. There is no firewall because there is no horizon in the fundamental description â€” there is only a region of high entanglement density on the boundary.

The Page curve follows automatically: as the boundary disentangles (Hawking evaporation), the entanglement entropy between the "black hole region" and the rest first increases (as new boundary sites become entangled) and then decreases (as the saturated region shrinks). This gives the Page curve without any need for replica wormholes, islands, or other machinery.

---

## 5. Testable Predictions and Observational Signatures

### 5.1 CMB large-angle anomalies (â„“ < 30)

**Status: CONJECTURED â€” prediction is structural, not quantitative**

If cycling occurs, the CMB at the largest angular scales should carry the imprint of the boundary's reorganization pattern. The boundary reorganization respects Vâ‚„ Ã— â„¤â‚ƒ symmetry (these are FORCED and cannot be broken). The dominant large-angle CMB modes should therefore have specific parity properties under Vâ‚„.

**Specific prediction:** The â„“ = 2 (quadrupole) and â„“ = 3 (octupole) modes, which are already known to exhibit anomalous alignment and planarity, should have power concentrated in Vâ‚„-even harmonics rather than Vâ‚„-odd. The hemispherical power asymmetry direction should align with the â„¤â‚ƒ axis of the boundary.

**Computation 5.1:** Decompose the observed CMB low-â„“ multipoles into Vâ‚„ Ã— â„¤â‚ƒ irreps and test for excess power in specific channels. This is a data analysis task, not a simulation.

### 5.2 Non-thermal Hawking radiation structure

**Status: CONJECTURED**

The deviations from thermality in Hawking radiation should carry the Vâ‚„ = â„¤â‚‚ Ã— â„¤â‚‚ fingerprint. The four Kâ‚„ orbitals contribute differently to the radiation depending on their active/inert classification. Active channels (winding âˆ’1) radiate with a spectral correction proportional to C = âˆ’2, while inert channels (winding 0) contribute only at higher order.

**Prediction:** Non-thermal corrections at O(Mâ»Â²_BH) with a specific 2:2 channel structure (two active, two inert contributions with different signs).

**Status: Not testable with current or foreseeable technology.** Included for theoretical completeness.

### 5.3 Time-dependent dark energy

**Status: DERIVED (from Part 1 Î› mechanism + Step A)**

If Î› tracks d_eff = 3 âˆ’ f(Îµ), and Îµ evolves as the boundary entangles, then Î›(t) is time-dependent. The direction of change depends on whether entropy increase causes Îµ to increase or decrease.

**Two scenarios:**
- If increased boundary entropy â†’ coarser effective resolution â†’ Îµ increases â†’ f(Îµ) decreases â†’ Î› decreases: dark energy weakens over cosmological time. (Consistent with asymmetry framework's Î›(t) âˆ âŸ¨dA/dtâŸ©.)
- If increased boundary entropy â†’ more bulk structure to resolve â†’ Îµ decreases â†’ f(Îµ) increases â†’ Î› increases: dark energy strengthens.

**Computation 5.3:** Determine the sign by computing d_eff(Îµ) for the Kâ‚„ lattice at different entanglement densities (different temperatures Î²). This is a well-posed lattice computation.

### 5.4 Aeon-invariant physics

**Status: DERIVED (from classification theorem + Step E)**

If cycling occurs, the following are the same in every aeon: gauge group, generation count, chirality, spacetime dimension. The following may differ: particle masses, coupling constants, mixing angles, cosmological constant.

**Prediction:** There is no "fine-tuning problem" for FORCED quantities (they are topologically protected across aeons). There may be a fine-tuning / anthropic selection problem for SILENT quantities (they are re-randomized each aeon).

---

## 6. Open Computations

| ID | Computation | Difficulty | Blocks |
|----|-------------|------------|--------|
| **A1** | Thermalization of Kâ‚„ boundary (exact diag, L=3,6) | LOW | Step A |
| **B1** | Fisher metric at Î² â†’ 0 (infinite T); â„¤â‚ƒ vestige? | LOW | Step B |
| **D1** | Persistent current from â„¤â‚ƒ flux at Î² â†’ 0 | LOW | Step D (critical) |
| **5.1** | Vâ‚„ Ã— â„¤â‚ƒ decomposition of observed CMB low-â„“ | MEDIUM | Prediction 5.1 |
| **5.3** | d_eff(Îµ) at different entanglement densities | MEDIUM | Prediction 5.3 |
| **Î›** | Large-Îµ asymptotics of f(Îµ) | MEDIUM | Part 1 open question |

**Priority order:** D1 first â€” if âŸ¨JâŸ© = 0 at Î² = 0, scenario D1a (no cycling) is viable and the entire Part 2 edifice is much weaker. If âŸ¨JâŸ© â‰  0, the topological obstruction to true equilibrium is established, and the cycling picture has a foundation.

---

## 7. What Must Be Killed

Following Part 1's discipline of killing claims that don't survive scrutiny:

### 7.1 Individual black holes spawn child universes

**KILLED.** The original asymmetry framework proposed this. Kâ‚„D eliminates it: there is no black hole interior for a child universe to form in. The bulk is emergent. Black holes are boundary entanglement saturation events, not portals.

### 7.2 The arrow of time resets between aeons

**KILLED.** Unlike Penrose CCC, Kâ‚„D's arrow of time is topological (â„¤â‚ƒ flux, C = âˆ’2). It cannot reset. What resets is entanglement structure, not directionality. Each aeon inherits the same arrow of time from the permanent â„¤â‚ƒ flux.

### 7.3 Specific coupling constants can be predicted from cycling

**KILLED.** SILENT quantities are re-randomized each aeon (if cycling occurs). They cannot be predicted from the cycling mechanism â€” they require knowing which moduli point the boundary reorganization selects, which is itself a SILENT question. Cycling explains *why* they have the values they have (random selection from the moduli space at each reorganization), not *what* values they have.

### 7.4 Î›(t) âˆ âŸ¨dA/dtâŸ© as originally formulated

**MODIFIED.** The asymmetry framework's formula had Î› tracking the mean asymmetry production rate. In Kâ‚„D, Î› tracks d_eff = 3 âˆ’ f(Îµ), which is related to but not identical to asymmetry production. The directionality is correct (Î› decreases as the universe approaches heat death), but the functional form needs to be derived from f(Îµ), not assumed from the asymmetry framework.

---

## 8. Status Summary

| Claim | Status | Depends On |
|-------|--------|------------|
| Boundary is fundamental, bulk emergent | **FORCED** | Part 1 |
| â„¤â‚ƒ flux permanent across any transition | **FORCED** | Axiom 3, C = âˆ’2 integer |
| Maximum boundary entropy exists | **DERIVED** | Lattice finiteness |
| Bulk approaches de Sitter at max entropy | **CONJECTURED** | Steps A, B |
| Cycling occurs | **OPEN** | Computation D1 |
| FORCED quantities survive across aeons | **FORCED** | Classification theorem |
| SILENT quantities re-randomize | **CONJECTURED** | Step E |
| Vâ‚„-structured CMB anomalies | **CONJECTURED** | Cycling + symmetry |
| Time-dependent Î› | **DERIVED** | Part 1 Î› mechanism |
| Page curve automatic | **DERIVED** | Boundary construction |
| No BH information paradox | **FORCED** | Part 1 |
| No BH interior | **FORCED** | Part 1 |

**The single most important computation is D1: does the â„¤â‚ƒ flux prevent true equilibrium at Î² = 0?** Everything downstream depends on this.

---

## 9. Relationship to Asymmetry Framework

The asymmetry framework (v2) was the precursor. Its core insight â€” that asymmetry is fundamental and spacetime emerges from it â€” survives and is sharpened by Kâ‚„D:

| Asymmetry Framework | Kâ‚„D Translation | Status |
|---------------------|------------------|--------|
| A_ij = P(iâ†’j) âˆ’ P(jâ†’i) > 0 | D â‰  D* (Axiom 3) | **Absorbed into Part 1** |
| Fisher metric signature flip | Forced by quadratic closure (Axiom 1) + D â‰  D* | **PROVED in Part 1** |
| BH as asymmetry accumulators | BH as boundary entanglement saturation | **Refined** |
| Child spacetime spawning from BH | KILLED â€” no BH interior | **Killed** |
| Parent-child cosmological spawning | Boundary reorganization at entropy saturation | **Refined into cycling** |
| Î›(t) âˆ âŸ¨dA/dtâŸ© | Î›(t) from d_eff = 3 âˆ’ f(Îµ) | **Modified** |
| Fundamental decoherence floor | Not addressed in Part 2; remains open | **Parked** |
| Incomplete measurement collapse | Not addressed; may connect to Axiom 1 + 4 | **Parked** |

The asymmetry framework was the intuition layer. Kâ‚„D Part 1 provided the mathematical scaffolding. Part 2 takes the cosmological intuitions from the asymmetry framework and runs them through the Kâ‚„D machinery to see what survives.

---

## 10. Next Steps

1. **Computation D1** (persistent current at Î² = 0). This is the gate. Small L, exact diag. Can be done immediately.
2. **Computation B1** (Fisher metric at Î² = 0). Determines what the bulk looks like at max entropy. Also small L.
3. **Computation 5.3** (d_eff vs entanglement density). Determines the sign of dÎ›/dt.
4. **Computation A1** (thermalization check). Verifies no unexpected protected sectors.
5. **Computation 5.1** (CMB Vâ‚„ decomposition). The only near-term observational test.
6. **Computation Î›** (large-Îµ asymptotics of f(Îµ)). Part 1 open question, now directly relevant.

Items 1â€“4 are small exact-diagonalization computations. Item 5 is data analysis. Item 6 requires more careful analytical work. None require GPU resources.
