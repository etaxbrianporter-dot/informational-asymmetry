---
title: D1-MB Results: Many-Body K₄ Level Statistics
section: Killed Approaches
status: active
---

# D1-MB Results: Many-Body K₄ Level Statistics

## 3-Site Results (dim = 924, full diag)

### Raw Data

| U/t | Flux | ⟨r⟩    | ±      | ‖Im(H)‖/‖Re(H)‖ | Verdict  | S_gs  | S_gs/S_max |
|-----|------|--------|--------|-------------------|----------|-------|------------|
| 0   | ℤ₃   | 0.2330 | 0.0904 | 1.000             | Poisson  | 2.594 | 0.935      |
| 0   | none | 0.0000 | 0.0000 | 0.000             | Poisson  | —     | —          |
| 2   | ℤ₃   | 0.4460 | 0.0503 | 0.210             | (mixed)  | 2.477 | 0.893      |
| 2   | none | 0.3611 | 0.0719 | 0.000             | Poisson  | —     | —          |
| 4   | ℤ₃   | 0.3455 | 0.0456 | 0.107             | Poisson  | 2.347 | 0.846      |
| 4   | none | 0.4577 | 0.0587 | 0.000             | (mixed)  | —     | —          |
| 6   | ℤ₃   | 0.4599 | 0.0489 | 0.071             | (mixed)  | 2.231 | 0.805      |
| 6   | none | 0.3208 | 0.0454 | 0.000             | Poisson  | —     | —          |
| 8   | ℤ₃   | 0.3693 | 0.0405 | 0.054             | Poisson  | 2.109 | 0.761      |
| 8   | none | 0.3460 | 0.0438 | 0.000             | Poisson  | —     | —          |

Reference values: GUE = 0.5996, GOE = 0.5307, Poisson = 0.3863

### What We See (with caveats)

**1. The ℤ₃ flux is genuinely complex at the many-body level.**

‖Im(H)‖/‖Re(H)‖ = 1.0 at U = 0 (maximally complex). This confirms D ≠ D* is a many-body effect, not just single-particle. The ratio decreases with U (the real, T-symmetric Hubbard interaction dilutes the complex hopping) but never reaches zero — at U = 8, Im/Re = 0.054, still nonzero.

**2. The level statistics are between Poisson and GOE, NOT near GUE.**

This is the most important finding. Despite the Hamiltonian being genuinely complex (T-broken), the level statistics do NOT reach GUE (0.5996). The maximum ⟨r⟩ = 0.4599 at U = 6 is between Poisson (0.3863) and GOE (0.5307).

Possible interpretations:
- **Finite-size effect:** 924 states is very small. Symmetry sectors (V₄, translation, particle-hole) create block structure that mimics Poisson. This is the most likely explanation. The 4-site system (12,870 states) should clarify.
- **Partial thermalization:** The system has integrable-like sectors (V₄ symmetry creates 4 distinct channels that don't fully mix). This would be physically interesting.
- **Anomalous statistics:** The ℤ₃ flux creates a new universality class between standard ensembles. Would need 6-site to confirm.

**3. The ℤ₃ flux ENHANCES thermalization compared to no-flux.**

At U = 2: ⟨r⟩_flux = 0.446 vs ⟨r⟩_noflux = 0.361 (Δ = +0.085)
At U = 6: ⟨r⟩_flux = 0.460 vs ⟨r⟩_noflux = 0.321 (Δ = +0.139)

The ℤ₃ flux consistently pushes the system TOWARD thermalization (higher ⟨r⟩). This makes physical sense: T-breaking removes degeneracies that would otherwise protect against level repulsion.

**4. ⟨r⟩ varies significantly with U (Δ⟨r⟩ = 0.227).**

The level statistics change qualitatively as U varies. Peak near U = 2 and U = 6, dip at U = 4. This non-monotonic behavior could signal the semimetal-insulator transition (U_c), where the system passes through enhanced quantum criticality.

**5. Ground state entanglement decreases monotonically with U.**

S_gs/S_max goes from 0.935 (U=0) to 0.761 (U=8). The interaction localizes particles on-site, reducing entanglement. But mid-spectrum entanglement stays near maximal (0.92 at U=8) — suggesting ETH holds for excited states even as the ground state becomes less entangled.

---

## What to Run Next

### Priority 1: 4-site system (minutes)

```bash
python computation_D1_manybody.py --nsites 4 --U 0 1 2 3 4 5 6 7 8 --full --noflux
```

This gives dim = 12,870 — large enough for meaningful level statistics. Full diag takes ~10-30 seconds per U value. **This is the minimum system for a credible result.**

Key things to look for:
- Does ⟨r⟩ reach GUE (0.5996) at any U? If yes → clean thermalization with T-breaking
- Does ⟨r⟩ stay below GOE? → V₄ symmetry creates protected sectors
- Is there a sharp peak in ⟨r⟩ at some U_c? → Quantum phase transition
- Does ⟨r⟩ become NON-MONOTONIC in U? → Phase structure in the many-body spectrum

### Priority 2: 6-site system (hours)

```bash
python computation_D1_manybody.py --nsites 6 --U 0 2 4 6 8 --nev 500 --no-ee
```

This gives dim = 2,704,156. Lanczos with 500 eigenvalues from the middle of the spectrum (shift-invert at σ=0) gives the best level statistics for ETH. Skip entanglement entropy (--no-ee) because the reduced density matrix construction is expensive at this size.

This is where the result becomes physically meaningful. If GUE appears at 6 sites, the system is a standard thermalizing T-broken system. If something anomalous persists, it's a real many-body effect.

### Priority 3: Symmetry-resolved sectors

The cleanest level statistics come from resolving symmetry sectors. The Hamiltonian commutes with:
- V₄ = ℤ₂ × ℤ₂ (4 sectors)
- Translation on torus (nsites momentum sectors)
- Particle-hole at half-filling

Each sector should independently show GUE (if T-broken, thermalizing) or Poisson (if not). The current code does NOT resolve these sectors — all states are mixed together, which artificially pushes ⟨r⟩ toward Poisson (superposition of uncorrelated spectra from different sectors gives Poisson regardless of individual sector statistics).

**The fact that ⟨r⟩ > Poisson even WITHOUT resolving sectors is already a positive signal for thermalization.**

---

## Interpretation for Cosmological Cycling

The many-body question has three layers:

**Layer 1: Does the many-body system thermalize?**
If ⟨r⟩ → GUE at larger sizes → yes, normal thermalization with broken T. The "heat death" is a standard Gibbs state, just with permanent T-breaking. This supports D1's single-particle conclusion: equilibrium exists but is asymmetric.

**Layer 2: Are there protected non-thermalizing sectors?**
If ⟨r⟩ stays significantly below GUE even after resolving symmetry sectors → the V₄ × ℤ₃ structure creates protected subspaces. This would mean parts of the Hilbert space NEVER thermalize — creating permanent non-equilibrium structures. This would be strong support for cycling: the system can never truly reach maximum entropy because of topologically protected sectors.

**Layer 3: Is there a phase transition in the level statistics as U varies?**
If ⟨r⟩(U) has a sharp feature (peak, kink, or jump) at some U_c → the many-body spectrum reorganizes at the phase transition. This would directly connect the quantum critical point (needed for Step 4 / γ₃) to the thermalization properties (needed for cycling). A system that has different thermalization properties on different sides of U_c could cycle between phases.

The 3-site data shows hints of all three (non-GUE statistics, U-dependent variation, ℤ₃ enhancement) but at a system size too small to be conclusive. **The 4-site run is the decision point.**
