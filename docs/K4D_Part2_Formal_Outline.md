---
title: Kâ‚„D Theory Part 2: Cosmological Cycling from Boundary Entropy Saturation
section: K₄ Spacetime
status: active
---

# Kâ‚„D Theory Part 2: Cosmological Cycling from Boundary Entropy Saturation

## Formal Paper Outline â€” February 2026

---

## Status Label Convention (from Part 1)

Every claim carries one of: **FORCED** (follows from proved structure, no additional assumptions), **DERIVED** (follows from Part 1 results plus a stated assumption), **CONJECTURED** (plausible but unproved), **OPEN** (well-posed question, answer unknown).

---

## Â§0. Relationship to Part 1

- Part 1 derived spacetime (gravity sector) and Standard Model (gauge sector) from matching combinatorics on Kâ‚„ and Kâ‚†
- Cosmological implications flagged but not developed:
  - Cosmological constant: d_eff = 3 âˆ’ f(Îµ), status REDUCED; f(Îµ) at large Îµ left as open computation
  - Black hole information paradox: status ELIMINATED (information cannot be lost inside an emergent bulk)
  - Dark matter/dark energy: status SILENT (no mechanism; framework determines algebraic structure, not cosmological content)
- **Part 2 central question:** What does the Kâ‚„D boundary theory say about cosmological evolution, entropy saturation, and cycling?

---

## Â§1. What Part 1 Forces About Cosmology

### 1.1 The boundary is fundamental, the bulk is emergent [FORCED]

- Kâ‚„ lattice model on triangular lattice is 2+1-dimensional
- Emergent 3+1 spacetime arises when a finite observer (Axiom 4) cannot resolve spectral structure; d_eff â†’ 3 at coarse resolution
- **Cosmological implication:** Any cosmological question must have a boundary-language answer; bulk-only questions are ill-posed
  - Ill-posed: "What happens inside a black hole?"
  - Well-posed: "What is the boundary entanglement configuration corresponding to a black hole?"

### 1.2 Asymmetry is topological and permanent [FORCED]

- â„¤â‚ƒ flux breaking time-reversal is topological: C = âˆ’2 is a Chern number (integer invariant, cannot be continuously deformed to zero)
- Arrow of time is algebraically prior to thermodynamics
- **Cosmological implication:** No cosmological process can eliminate the arrow of time; any cycling mechanism must preserve D â‰  D*; â„¤â‚ƒ flux survives whatever transition occurs between aeons

### 1.3 Entanglement entropy is bounded [DERIVED]

- Boundary is Kâ‚„ lattice on triangular torus L Ã— L, each site has 4 orbitals
- Maximum entanglement entropy: S_max = 4LÂ² Â· ln 2
- **Cosmological implication:** The boundary can reach maximum entropy in finite time â€” this IS heat death in boundary language

### 1.4 The effective dimension encodes cosmological constant [DERIVED]

- d_eff = 3 âˆ’ f(Îµ) where f(Îµ) depends on the observer's spectral resolution
- Î› tracks d_eff; as boundary homogenizes, f(Îµ) â†’ 0 and Î› â†’ 0
- **Cosmological implication:** Î›(t) is time-dependent, tied to boundary entanglement evolution

---

## Â§2. The Cycling Mechanism (Steps Aâ€“E)

### Step A: The boundary thermalizes [CONJECTURED, verifiable]

- Second law follows from same mechanism as any quantum lattice system
- â„¤â‚ƒ flux breaks T-reversal â†’ actually *helps* thermalization (prevents MBL from T-reversal-protected resonances)
- Vâ‚„ symmetry is small (4 elements) â†’ unlikely to create large protected subspaces
- **Computation A1:** Exact diag at L = 3, 6; check monotonic growth of entanglement entropy from generic initial state; estimate thermalization time vs U/t

### Step B: Maximum entropy configuration of the boundary [CONJECTURED]

- **Claim B1:** At max entropy, two-point correlator âŸ¨câ€ _i c_jâŸ© decays exponentially at shortest possible scale in every Vâ‚„ channel
- **Claim B2:** Fisher information metric on (kâ‚, kâ‚‚, E) becomes trivial â†’ emergent bulk is maximally homogeneous and isotropic
- **Key question:** Does residual â„¤â‚ƒ structure survive in the Fisher metric at Î² â†’ 0?
- **Computation B1:** Compute g_ij for Kâ‚„ lattice at infinite temperature; check whether â„¤â‚ƒ vestigial anisotropy persists

### Step C: What the bulk observer sees at maximum boundary entropy [DERIVED from A + B]

- d_eff â†’ 3, f(Îµ) â†’ 0 (pure 3D space, no energy-space mixing)
- Î› â†’ 0 (cosmological constant vanishes as boundary homogenizes)
- Emergent geometry: approximately de Sitter with decreasing Î›
- All massive particles have decayed; only topological invariants survive: â„¤â‚ƒ flux, C = âˆ’2, Vâ‚„ algebra
- **This is Penrose's CCC precondition, derived from boundary physics**

### Step D: The transition mechanism [CONJECTURED â€” hardest and most speculative]

- **Question D1:** Is the maximum entropy state of the Kâ‚„ boundary *stable* or *unstable*?
- Three scenarios:
  - **D1a: Stable equilibrium.** Heat death is permanent. No cycling. Conservative answer.
  - **D1b: Unstable via boundary phase transition.** Semimetal-insulator transition at U = U_c; if effective coupling crosses phase boundary as entanglement accumulates â†’ first-order catastrophic boundary reorganization â†’ bulk observer perceives new big bang
  - **D1c: Unstable via topological mechanism.** D â‰  D* means system can never reach detailed balance; maximum entropy state is not true equilibrium; entropy production with nowhere to absorb â†’ pressure for reorganization
- **Computation D1:** Study Kâ‚„ lattice at Î² â†’ 0; does â„¤â‚ƒ flux create residual current or Berry phase preventing true Gibbs equilibrium? If âŸ¨JâŸ© â‰  0 at Î² = 0 â†’ D1a eliminated

#### D1 Results (executed Feb 17, 2026 â€” single-particle)

- **Spectral asymmetry (clincher):** Mean |Îµ(k) âˆ’ Îµ(âˆ’k)|/|Îµ(k)| = **20.25%** with â„¤â‚ƒ flux vs 0% without; 96% of BZ has asymmetry > 1%
- This is a property of the Hamiltonian, not the quantum state â†’ survives at all temperatures
- **Verdict:** â„¤â‚ƒ flux creates permanent spectral asymmetry that prevents true detailed balance

#### D1-MB Results (3-site many-body, dim = 924)

- â„¤â‚ƒ flux is genuinely complex at many-body level (â€–Im(H)â€–/â€–Re(H)â€– = 1.0 at U = 0)
- Level statistics between Poisson and GOE, NOT near GUE â€” possible anomalous thermalization
- â„¤â‚ƒ flux consistently ENHANCES thermalization (higher âŸ¨râŸ©) compared to no-flux
- **Interpretation layers:** (1) Does system thermalize? (2) Are there protected non-thermalizing sectors? (3) Is there a phase transition in level statistics vs U?
- **4-site run is the decision point** (12,870 states)

### Step E: Memory across aeons [CONJECTURED, depends on C + D]

- **E1: FORCED quantities survive** â€” gauge group, generation count, chirality, d = 3+1, C = âˆ’2, â„¤â‚ƒ flux (graph-topological, independent of moduli/entanglement state)
- **E2: SILENT quantities reset** â€” masses, couplings, mixing angles (moduli-geometric, depend on specific entanglement pattern â†’ scrambled at reorganization)
- **E3: REDUCED quantities are the interesting case** â€” Chern numbers of specific bands, mass hierarchy pattern take discrete values depending on moduli chamber; reorganization may or may not change chamber
- **Implication:** Same laws, different parameters each aeon; structural/parametric divide from Part 1 classification becomes aeon-invariant/aeon-variable divide

---

## Â§3. Comparison with Penrose CCC

| Feature | Penrose CCC | Kâ‚„D Cycling |
|---------|-------------|-------------|
| End of aeon | Massive particles decay; only conformal physics | Boundary max entanglement; Fisher metric trivial |
| Transition | Conformal rescaling (mathematical) | Boundary phase transition or topological instability (physical) |
| What carries over | Conformal structure, Hawking points | FORCED invariants (gauge group, generations, chirality, C = âˆ’2) |
| What resets | Masses, scales | SILENT quantities (masses, couplings, mixing angles) |
| Observational signature | Circular CMB anomalies from prior-aeon BHs | Vâ‚„-structured large-angle CMB anomalies from boundary reorganization |
| Arrow of time | Resets via conformal rescaling of Weyl tensor | **Never resets** â€” â„¤â‚ƒ flux is permanent (Axiom 3) |

**Key advantage:** Kâ‚„D does not need a separate hypothesis for why entropy is low at the start of each aeon â€” the boundary just reorganized, and the new pattern hasn't yet thermalized.

---

## Â§4. The Black Hole Role (Revised)

### 4.1 What a black hole IS in Kâ‚„D [DERIVED]

- Region of 2+1 boundary where entanglement has reached local saturation
- Boundary DOF maximally entangled with each other â†’ no entanglement capacity for long-range correlations
- Bulk observer perceives region from which no information escapes (not trapped behind horizon â€” maximally entangled boundary looks featureless to finite observer)
- Bekenstein-Hawking entropy S_BH = A/(4l_pÂ²) is entanglement entropy of saturated boundary region; area law automatic from 2+1 boundary â†” 3+1 bulk mapping

### 4.2 Hawking evaporation as boundary disentanglement [DERIVED]

- Maximally entangled region not in true equilibrium (D â‰  D* prevents detailed balance)
- Entanglement leaks from saturated region â†’ bulk observer sees thermal radiation
- Temperature set by entanglement gradient at boundary of saturated region

### 4.3 Black holes as cosmological entropy accelerators [CONJECTURED]

- BH formation-and-evaporation cycles homogenize boundary entanglement distribution
- Not "recycle bins" or universe-spawners â€” they are **mixers**
- Rate of BH formation controls timescale of current aeon

### 4.4 No interior, no information paradox, no firewall [FORCED]

- No BH interior in Kâ‚„D (bulk is emergent)
- Page curve automatic: boundary disentanglement during evaporation gives correct entropy trajectory
- No need for replica wormholes, islands, or other machinery

---

## Â§5. Testable Predictions and Observational Signatures

### 5.1 CMB large-angle anomalies (â„“ < 30) [CONJECTURED â€” structural, not quantitative]

- If cycling occurs, CMB at largest angular scales carries imprint of boundary reorganization
- Reorganization respects Vâ‚„ Ã— â„¤â‚ƒ symmetry (FORCED)
- **Prediction:** â„“ = 2 (quadrupole) and â„“ = 3 (octupole) should have power concentrated in Vâ‚„-even harmonics; hemispherical asymmetry direction aligns with â„¤â‚ƒ axis
- **Computation 5.1:** Decompose observed low-â„“ multipoles into Vâ‚„ Ã— â„¤â‚ƒ irreps; test for excess power in specific channels (data analysis, not simulation)

### 5.2 Non-thermal Hawking radiation structure [CONJECTURED â€” not testable]

- Deviations from thermality carry Vâ‚„ = â„¤â‚‚ Ã— â„¤â‚‚ fingerprint
- Two active channels (winding âˆ’1) + two inert channels (winding 0)
- Non-thermal corrections at O(Mâ»Â²_BH) with specific 2:2 channel structure
- Included for theoretical completeness

### 5.3 Time-dependent dark energy [DERIVED from Part 1 Î› mechanism + Step A]

- If Î› tracks d_eff = 3 âˆ’ f(Îµ), and Îµ evolves with boundary entanglement â†’ Î›(t) is time-dependent
- **Two scenarios depending on sign:**
  - Increased entropy â†’ coarser resolution â†’ Îµ increases â†’ Î› decreases (dark energy weakens)
  - Increased entropy â†’ more bulk structure â†’ Îµ decreases â†’ Î› increases (dark energy strengthens)
- **Computation 5.3:** d_eff(Îµ) at different entanglement densities (different Î²) â€” well-posed lattice computation

### 5.4 Aeon-invariant physics [DERIVED from classification theorem + Step E]

- Same in every aeon: gauge group, generation count, chirality, spacetime dimension
- May differ: particle masses, coupling constants, mixing angles, Î›
- No fine-tuning problem for FORCED quantities (topologically protected across aeons)
- Fine-tuning / anthropic selection applies only to SILENT quantities (re-randomized each aeon)

### 5.5 Connection to Paper II predictions (birefringence)

- Paper II derives cosmic birefringence from Kâ‚„ spectral action: Î±â‚€ â‰¤ 0.41Â° (parameter-free)
- Five falsifiable predictions for LiteBIRD/CMB-S4:
  1. Birefringence dipole aligned with TT asymmetry
  2. Dipolar amplitude Î´Î±/Î±â‚€ = (2âˆ’âˆš3)Â²
  3. EB BiPoSH parity alternation (âˆ’1)^{â„“+1}
  4. TT/EB shape split with sign change at â„“ ~ 20
  5. TB/EB ratio = C^{TE}_â„“ / C^{EE}_â„“
- Part 2 provides the cosmological context: the torsion pseudoscalar arises naturally from the boundary's â„¤â‚ƒ flux

---

## Â§6. What Must Be Killed

Following Part 1's discipline of killing claims that don't survive:

1. **Individual BHs spawn child universes** â€” KILLED. No BH interior for child universe to form in.
2. **Arrow of time resets between aeons** â€” KILLED. â„¤â‚ƒ flux is topological (C = âˆ’2). Entanglement resets, directionality doesn't.
3. **Specific coupling constants predictable from cycling** â€” KILLED. SILENT quantities re-randomized each aeon.
4. **Î›(t) âˆ âŸ¨dA/dtâŸ© as originally formulated** â€” MODIFIED. Î› tracks d_eff = 3 âˆ’ f(Îµ), not mean asymmetry rate. Directionality correct, functional form needs derivation from f(Îµ).

---

## Â§7. Relationship to Asymmetry Framework

| Asymmetry Framework | Kâ‚„D Translation | Status |
|---------------------|------------------|--------|
| A_ij = P(iâ†’j) âˆ’ P(jâ†’i) > 0 | D â‰  D* (Axiom 3) | Absorbed into Part 1 |
| Fisher metric signature flip | Forced by quadratic closure + D â‰  D* | PROVED in Part 1 |
| BH as asymmetry accumulators | BH as boundary entanglement saturation | Refined |
| Child spacetime spawning from BH | â€” | **Killed** |
| Parent-child cosmological spawning | Boundary reorganization at entropy saturation | Refined into cycling |
| Î›(t) âˆ âŸ¨dA/dtâŸ© | Î›(t) from d_eff = 3 âˆ’ f(Îµ) | Modified |
| Fundamental decoherence floor | Not addressed in Part 2 | Parked |
| Incomplete measurement collapse | May connect to Axiom 1 + 4 | Parked |

---

## Â§8. Status Summary

| Claim | Status | Depends On |
|-------|--------|------------|
| Boundary fundamental, bulk emergent | **FORCED** | Part 1 |
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

---

## Â§9. Open Computations

| ID | Computation | Difficulty | Blocks | Status |
|----|-------------|------------|--------|--------|
| **D1** | Persistent current from â„¤â‚ƒ flux at Î² â†’ 0 | LOW | Step D (critical) | **DONE** â€” single-particle confirms spectral asymmetry; 3-site MB done, 4-site pending |
| **D1-MB** | Many-body level statistics, symmetry-resolved sectors | MEDIUM | Step D | 3-site done; 4-site is decision point |
| **A1** | Thermalization of Kâ‚„ boundary (exact diag, L=3,6) | LOW | Step A | Not started |
| **B1** | Fisher metric at Î² â†’ 0; â„¤â‚ƒ vestige? | LOW | Step B | Not started |
| **5.1** | Vâ‚„ Ã— â„¤â‚ƒ decomposition of observed CMB low-â„“ | MEDIUM | Prediction 5.1 | Not started |
| **5.3** | d_eff(Îµ) at different entanglement densities | MEDIUM | Prediction 5.3 | Not started |
| **Î›** | Large-Îµ asymptotics of f(Îµ) | MEDIUM | Part 1 open question | Not started |

**Priority:** D1-MB (4-site) â†’ B1 â†’ A1 â†’ 5.3 â†’ 5.1 â†’ Î›

---

## Â§10. Discussion / Outlook

### 10.1 What Part 2 adds to Part 1

- Part 1 derived the *structure* of physics (gauge group, generations, spacetime); Part 2 asks about the *evolution*
- The boundary construction makes cosmological questions well-posed in a way that standard QFT+GR does not
- The classification theorem's FORCED/SILENT/REDUCED distinction acquires cosmological meaning as the aeon-invariant/aeon-variable/aeon-discrete divide

### 10.2 The single gate computation

- Everything in Part 2 stands or falls on whether the â„¤â‚ƒ flux prevents true equilibrium at Î² = 0
- D1 single-particle: 20% spectral asymmetry â†’ **strong positive signal**
- D1-MB 3-site: â„¤â‚ƒ enhances thermalization, statistics between Poisson and GUE â†’ **inconclusive, need 4-site**
- If 4-site confirms anomalous statistics in symmetry-resolved sectors â†’ cycling picture has a foundation
- If 4-site shows standard GUE in all sectors â†’ D1a (no cycling) viable; Part 2 reduces to BH reinterpretation + time-dependent Î›

### 10.3 Items parked for future work

- Fundamental decoherence floor (from asymmetry framework)
- Incomplete measurement collapse (may connect to Axioms 1 + 4)
- Sign problem diagnostic for DQMC at larger system sizes
- Connection to Kâ‚† normalization / Higgs mass (c = Ï€Â²/8 question)
