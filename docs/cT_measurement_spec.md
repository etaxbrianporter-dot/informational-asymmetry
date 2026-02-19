---
title: c_T Measurement Specification for the Kâ‚„ Critical Point
section: Gravitational Sector
status: active
---

# c_T Measurement Specification for the Kâ‚„ Critical Point

## The Shortcut: Measure c_T, Not Î”â‚ƒ

**February 17, 2026**

---

## The Logic Chain

$$c_T \xrightarrow{\text{QMC}} N_{\text{eff}} = \frac{c_T}{c_T^{\text{free}}} \xrightarrow{\text{O(N) bootstrap}} \gamma_3$$

Or better:

$$(\Delta_\sigma, c_T) \xrightarrow{\text{QMC}} \text{SDPB bootstrap} \xrightarrow{\text{rigorous}} \Delta_3 \geq f(\Delta_\sigma, c_T)$$

Instead of measuring the spin-3 anomalous dimension directly (hard, delicate, expensive), measure two simple quantities â€” the stress-tensor central charge c_T and the scalar dimension Î”_Ïƒ â€” then feed them into the conformal bootstrap for a rigorous bound on Î”â‚ƒ.

---

## 1. What c_T Is

In any 3d CFT, the stress-tensor two-point function is fixed by conformal symmetry up to one number:

$$\langle T_{\mu\nu}(x) \, T_{\rho\sigma}(0) \rangle = \frac{c_T}{S_3^2} \cdot \frac{I_{\mu\nu,\rho\sigma}(x)}{|x|^{6}}$$

c_T counts the effective degrees of freedom. For free fields in d = 3 (Osborn-Petkou normalization):

| Field | c_T |
|-------|-----|
| 1 real scalar | 3/(32Ï€Â²) â‰ˆ 0.00950 |
| 1 two-component Dirac | 0.01900 |
| Kâ‚„ free (all 4 channels) | 0.07599 (= 8 scalars) |
| Kâ‚„ free (2 active channels) | 0.03800 (= 4 scalars) |

Define N_eff â‰¡ c_T / c_T^{scalar}. At the free point: N_eff = 8 (all channels) or 4 (active only). At the interacting critical point: N_eff < N_eff^free by the F-theorem.

### Why N_eff determines Î³â‚ƒ

From the O(N) Wilson-Fisher bootstrap data (Kos-Poland-Simmons-Duffin-Vichi, Chester et al.), the spin-3 anomalous dimension as a function of N is precisely known:

| N_eff | Î³â‚ƒ | Î”â‚ƒ | vs threshold |
|-------|-----|-----|-------------|
| 0.8 | 3.42 | 7.42 | Well above |
| 1.0 | 2.50 | 6.50 | Above |
| **1.15** | **2.00** | **6.00** | **Threshold** |
| 1.5 | 1.46 | 5.46 | Below |
| 2.0 | 1.02 | 5.02 | Below |

With CS parity-breaking enhancement (factor 1.3â€“1.5), the threshold shifts to N_eff â‰ˆ 1.5.

**If c_T/c_T^free < 0.30 (relative to active channels): Î³â‚ƒ > 2 confirmed.**

---

## 2. The Lattice Stress-Tensor Operator

### Bond energy operator

For each nearest-neighbor direction Î´áµ¢ on the triangular lattice:

$$B_i(x,\tau) = \zeta_i \cdot c^\dagger(x,\tau) \, M_i \, c(x+\delta_i,\tau) + \text{h.c.}$$

where Mâ‚, Mâ‚‚, Mâ‚ƒ are the Kâ‚„ matching matrices, Î¶â‚ = 1, Î¶â‚‚ = Ï‰, Î¶â‚ƒ = Ï‰Â². This is the kinetic energy along bond direction i â€” a **standard DQMC observable** already computed in every Hubbard model simulation.

### Spin-2 combination (stress tensor)

$$T_+(x,\tau) = B_1(x,\tau) + \tilde\omega \cdot B_2(x,\tau) + \tilde\omega^2 \cdot B_3(x,\tau)$$

where Ï‰Ìƒ = e^{2Ï€i/3}. Under Câ‚ƒ rotation: Tâ‚Š â†’ e^{4Ï€i/3} Tâ‚Š (angular momentum +2), confirming spin-2.

In Cartesian components:
- T_{xx} âˆ’ T_{yy} = (2/3)[2Bâ‚ âˆ’ Bâ‚‚ âˆ’ Bâ‚ƒ]  
- T_{xy} = (1/âˆš3)[Bâ‚‚ âˆ’ Bâ‚ƒ]

### Implementation note

The bond operators Báµ¢ require only nearest-neighbor Green's functions âŸ¨câ€ (x)c(x+Î´)âŸ©, which are standard DQMC outputs. The stress tensor Tâ‚Š is a linear combination of these â€” **no new code needed**, only post-processing of existing observables.

---

## 3. Extraction Methods

### Method A: Ratio to free theory (RECOMMENDED)

1. Measure G_T(Ï„) = Lâ»Â² Î£_x âŸ¨Tâ‚Š(x,Ï„) Tâ‚‹(0,0)âŸ© at U = U_c
2. Compute G_T^free(Ï„) at U = 0 **exactly** from the Kâ‚„ band structure (no QMC needed)
3. Form the ratio R(Ï„, L) = G_T(U_c) / G_T(U=0)
4. At intermediate Ï„, R(Ï„, L) develops a **plateau** â†’ read off c_T/c_T^free directly

The ratio cancels: normalization conventions, leading lattice artifacts, tensor structure prefactors. It has a **plateau** at intermediate Ï„ â€” no fitting required, just read off the value.

### Method B: Exponential fit

At Î² â‰« L (ground state regime), G_T(Ï„) âˆ exp(âˆ’6Ï€Ï„/L) with the exponent **exactly known** (Î”_T = 3, Ward identity). Fit only the coefficient â†’ c_T. This is a **one-parameter fit** (vs. two parameters for the spin-3 measurement).

### Method C: Spatial correlator (cross-check)

G_T(r) = c_T Ã— Aâ‚ƒ / râ¶ with the exponent 6 = 2Î”_T **exact**. Fit coefficient only.

---

## 4. Why This Is Simpler Than Measuring Î”â‚ƒ

|  | c_T measurement | Î”â‚ƒ measurement |
|--|----------------|-----------------|
| **Lattice operator** | Nearest-neighbor bond energy (standard) | 2nd-neighbor traceless rank-3 tensor (delicate) |
| **Operator mixing** | None â€” spin-2 unique in its channel | Severe â€” spin-3 mixes with spin-1 |
| **Unknown exponent** | NO â€” Î”_T = 3 exact (Ward identity) | YES â€” Î”â‚ƒ is the unknown target |
| **Fit parameters** | 1 (coefficient only) | 2 (exponent + coefficient) |
| **Required L** | 6, 9, 12, 15, 18 | 12, 18, 24, 30, 36 |
| **Sweeps per point** | 5,000 â€“ 10,000 | 20,000 â€“ 50,000 |
| **Signal** | Strong (râ»â¶ decay) | Weak (râ»â¸ or worse) |
| **GPU-hours** | ~100â€“300 | ~2,000â€“4,000 |
| **Timeline** | ~1 week | ~2â€“3 weeks |

**The c_T shortcut is 5â€“10Ã— cheaper** and avoids the most technically challenging part of the original Computation 2 (constructing the spin-3 lattice operator with proper trace subtraction).

---

## 5. Computational Parameters

| Parameter | Value |
|-----------|-------|
| Model | Kâ‚„ Hubbard, triangular lattice, Zâ‚ƒ flux |
| Orbitals per site | 4 |
| Filling | Half-filling (n = 2) |
| Coupling | U_c, U_c Â± Î´U |
| Lattice sizes | L = 6, 9, 12, 15, 18 |
| Inverse temperature Î² | 3L |
| Trotter step Î”Ï„ | 0.05 |
| Warmup sweeps | 1,000 |
| Measurement sweeps | 5,000 â€“ 10,000 |
| **Primary observable** | G_T(Ï„) = Lâ»Â² Î£_x âŸ¨Tâ‚Š(x,Ï„) Tâ‚‹(0,0)âŸ© |
| **Secondary observable** | Câ‚€(r) = âŸ¨Oâ‚€(x) Oâ‚€(0)âŸ© with Oâ‚€ = Î£_Î± n_Î± âˆ’ 2 |
| **Free reference** | G_T^free(Ï„) from exact band structure (no QMC) |
| **Extraction** | R(Ï„,L) = G_T(U_c)/G_T(U=0) â†’ plateau â†’ N_eff |

Total runs: 5 sizes Ã— 3 couplings = **15 DQMC runs** + 5 exact free-field computations.

---

## 6. Error Budget

| Source | Î´(c_T/c_T^free) | Controllable? |
|--------|-----------------|---------------|
| QMC statistics | ~3â€“5% | Yes (more sweeps) |
| U_c uncertainty | ~2â€“3% | Yes (3-point bracket) |
| Finite-size scaling | ~3â€“5% | Yes (5 lattice sizes) |
| Trotter error | ~1% | Yes (Î”Ï„ extrapolation) |
| Lattice artifacts | ~1â€“2% | Mostly (ratio cancels) |
| **Total on c_T/c_T^free** | **~5â€“8%** | |
| â†’ **Î´N_eff** | **~Â±0.3** | |

For the Î³â‚ƒ determination:
- Via O(N) curve alone: Î´Î³â‚ƒ â‰ˆ Â±0.5 (includes GNâ†’O(N) mapping uncertainty)
- Via SDPB with measured (Î”_Ïƒ, c_T): Î´Î³â‚ƒ â‰ˆ Â±0.3 (**bypasses** the mapping)

---

## 7. The Optimal Strategy: QMC + Bootstrap

The key insight: **c_T alone** uses the GNâ†’O(N) mapping (20% uncertainty). But **c_T plus Î”_Ïƒ** feeds directly into the conformal bootstrap without any mapping:

$$\text{SDPB}: \quad \Delta_3 \geq f(\Delta_\sigma, c_T) \quad \text{(rigorous bound)}$$

### Workflow

1. **QMC** measures c_T/c_T^free and Î”_Ïƒ [1 week]
2. **SDPB** computes the lower bound on Î”â‚ƒ at the measured (Î”_Ïƒ, c_T) point [2 weeks, parallel]
3. If the bound gives Î”â‚ƒ > 6: **Î³â‚ƒ > 2 is proved rigorously** [done]

The bootstrap inputs are fully specified by the Kâ‚„ algebra:
- d = 3, Zâ‚‚ global symmetry (from active Vâ‚„ channels)
- No parity (from Zâ‚ƒ flux, C = âˆ’2)
- Stress tensor at Î” = 3 (Ward identity)
- Conserved Zâ‚‚ current at Î” = 2 (Vâ‚„ symmetry)
- External scalar Î”_Ïƒ = [measured]
- c_T = [measured]

**Neither QMC nor bootstrap alone closes the gap. Together they should.** QMC gives what it does best (simple ground-state correlators). Bootstrap gives what it does best (rigorous operator dimension bounds).

---

## 8. Decision Tree

```
Phase 0: Sign check
    â”‚
    â”œâ”€â”€ âŸ¨signâŸ© > 0.1 â”€â”€â†’ Comp 1: Find U_c
    â”‚                         â”‚
    â”‚                    â”œâ”€â”€ Continuous â”€â”€â†’ c_T + Î”_Ïƒ measurement
    â”‚                    â”‚                      â”‚
    â”‚                    â”‚                 â”œâ”€â”€ N_eff < 1.3
    â”‚                    â”‚                 â”‚     Î³â‚ƒ > 2 CONFIRMED
    â”‚                    â”‚                 â”‚     â†’ STEP 4 CLOSED
    â”‚                    â”‚                 â”‚
    â”‚                    â”‚                 â”œâ”€â”€ 1.3 < N_eff < 1.8
    â”‚                    â”‚                 â”‚     MARGINAL
    â”‚                    â”‚                 â”‚     â†’ Feed into SDPB
    â”‚                    â”‚                 â”‚          â”‚
    â”‚                    â”‚                 â”‚     â”œâ”€â”€ Î”â‚ƒ > 6 proved â†’ CLOSED
    â”‚                    â”‚                 â”‚     â””â”€â”€ Î”â‚ƒ < 6 allowed â†’ Full Comp 2
    â”‚                    â”‚                 â”‚
    â”‚                    â”‚                 â””â”€â”€ N_eff > 1.8
    â”‚                    â”‚                       Î³â‚ƒ < 2 LIKELY
    â”‚                    â”‚                       â†’ Reconsider program
    â”‚                    â”‚
    â”‚                    â””â”€â”€ First-order â”€â”€â†’ No CFT. Program pause.
    â”‚
    â””â”€â”€ âŸ¨signâŸ© < 0.1 â”€â”€â†’ Switch to DMRG/iDMRG
```

---

## 9. Expected Outcomes

Three independent lines of evidence on the expected N_eff:

**From entanglement (proved):** L_AdS/G_N = 4/3 = O(1). In holographic theories, this ratio scales as c_T. O(1) means N_eff ~ 1â€“2. Points to **Scenario A** (strongly interacting).

**From analogy (SU(2) Hubbard on honeycomb):** Î”_Ïƒ â‰ˆ 0.88. c_T/c_T^free â‰ˆ 0.4â€“0.6 for similar strongly-coupled Dirac systems. N_eff â‰ˆ 1.6â€“2.4.

**From self-duality (proved):** Î» = 1/2 means the theory is at the maximally interacting point for its symmetry class. Pushes c_T/c_T^free toward the lower end. N_eff pushed toward 1.0â€“1.5.

**Combined expectation:** N_eff â‰ˆ 1.0â€“2.0, central value ~1.5. Via the O(N) curve with CS enhancement: Î³â‚ƒ â‰ˆ 2.0â€“3.0. **Likely resolves the question**, especially via the QMC+bootstrap path.

---

## 10. Timeline and Cost

| Week | Activity | GPU-hours |
|------|----------|-----------|
| 0 | Phase 0 sign check | ~50 |
| 1 | Computation 1 (find U_c) | ~200â€“500 |
| 2 | c_T + Î”_Ïƒ measurement at U_c | ~100â€“300 |
| 2â€“3 | SDPB bootstrap with measured inputs (parallel) | ~500â€“2000 CPU |
| 3 | Final analysis and Î³â‚ƒ determination | ~10 |
| **Total** | | **~350â€“850 GPU + ~2000 CPU** |

Compare to the original Computation 2 (direct Î”â‚ƒ): ~2,000â€“4,000 GPU-hours plus delicate operator construction.

---

## 11. Deliverables

**From QMC (Week 2):**

(a) c_T/c_T^free with statistical + systematic errors  
(b) N_eff with error bars  
(c) Î”_Ïƒ from density-density correlator with error bars  
(d) Confirmation: continuous vs. first-order transition  
(e) Consistency check: Î”_T = 3 (Ward identity verification)

**From Analysis (Week 2, same day):**

(f) Î³â‚ƒ estimate from O(N) bootstrap curve Â± propagated errors  
(g) GO/NO-GO decision for full Computation 2 vs. bootstrap path

**From SDPB Bootstrap (Week 3):**

(h) Rigorous lower bound on Î”â‚ƒ at measured (Î”_Ïƒ, c_T)  
(i) If Î”â‚ƒ > 6 proved: **Î³â‚ƒ > 2 ESTABLISHED â†’ Step 4 CLOSED**  
(j) If not: specific Î”â‚ƒ bound value for program assessment

---

## 12. Relationship to Existing Documents

This specification **replaces** the direct Î”â‚ƒ measurement as the first computational step after finding U_c. It does not eliminate Computation 2 â€” it provides a fast triage that either:

- **Resolves the question** (if N_eff < 1.3 or if SDPB proves Î”â‚ƒ > 6), making Computation 2 unnecessary, or
- **Constrains the question** (if N_eff is marginal), making Computation 2 targeted and more efficient (known N_eff constrains the expected Î”â‚ƒ range)

The original computation_requirements_v2 document remains valid for the full Computation 2 if needed. This document specifies the **new fast path** that should be executed first.
