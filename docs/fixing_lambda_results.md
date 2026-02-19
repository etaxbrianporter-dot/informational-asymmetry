---
title: Fixing Λ: What Falls Out
section: K₆ Higgs Sector
status: active
---

# Fixing Λ: What Falls Out

## Brian Porter — February 2026

---

## The Question

One free parameter remains: Λ (the lattice/GUT scale). Can the gauge sector of the spectral action on K₆ × K₄ fix it?

## The Answer

**Λ was never a free parameter of the Higgs mass prediction.** It is a free parameter of the *hierarchy problem*.

The computation reveals that the Higgs mass formula m_H = √(2a₄/a₂) × m_W involves:
- a₄/a₂ = combinatorial invariant of K₆ matchings (Λ-independent)
- m_W = experimental input

Λ cancels in the ratio. The BZ averaging already contains all momentum-shell integration (= RG running). The T² lattice scale IS the EW scale.

---

## What Λ Determines

Λ enters the framework through the gravitational sector:

| Quantity | Formula | Role |
|:---|:---|:---|
| Planck mass | M_Pl² ∝ f₂Λ²a₂ | Hierarchy problem |
| Cosmological constant | Λ_CC ∝ f₄Λ⁴a₀ | CC problem |
| Higgs mass | m_H = √(2a₄/a₂) × m_W | **Λ drops out** |

The remaining free parameter f₂Λ² determines M_Pl/m_W — the gauge hierarchy — but NOT the Higgs mass.

---

## Five Key Results from the Computation

### 1. SM couplings don't unify (well-known)

| Crossing | Scale | log₁₀(μ/GeV) |
|:---|:---|:---|
| g₁ = g₂ | 1.09 × 10¹³ | 13.04 |
| g₁ = g₃ | 1.72 × 10¹⁴ | 14.23 |
| g₂ = g₃ | 2.89 × 10¹⁶ | 16.46 |

Triangle gap spans 3.4 decades. No unification without BSM physics.

### 2. The 7/48 boundary condition lives at the TeV scale

The democratic K₆ ratio λ/g₂² = 7/48 = 0.1458 is reached by SM running at **μ ≈ 2.8 TeV** — not at a GUT scale. The ratio starts at λ/g₂² = 0.304 at m_Z and falls through 7/48 around a TeV before λ goes negative (vacuum instability at ~10⁶ GeV).

This is physically meaningful: 7/48 from the democratic K₆ point is a *threshold* condition near the EW scale, not a high-energy boundary condition.

### 3. The CCM approach (7/48 at GUT scale) gives m_H ≈ 153 GeV

Setting λ/g₂² = 7/48 at the g₂-g₃ crossing (~10¹⁶·⁵ GeV) and running down:

| Λ_GUT | m_H (predicted) |
|:---|:---|
| 10¹³ GeV | 152.6 GeV |
| 10¹⁴ GeV | 152.8 GeV |
| 10¹⁵ GeV | 153.0 GeV |
| 10¹⁶ GeV | 153.1 GeV |
| 10¹⁷ GeV | 153.2 GeV |

The prediction is **insensitive to Λ_GUT** and always gives ~153 GeV — the quasi-fixed point of the top Yukawa. This is the CCM (2006) problem: running from a GUT-scale boundary condition always overshoots.

The BZ formula avoids this entirely because BZ averaging IS the momentum-shell integration.

### 4. The BZ formula reproduces experiment

| Quantity | BZ formula | Experiment |
|:---|:---|:---|
| λ/g₂² | 0.3076 | 0.3038 |
| m_H | 126.10 GeV | 125.09 GeV |
| Residual | +1.01 GeV | — |

Zero free parameters. 0.81% precision. The residual is consistent with O(1%) corrections (2-loop, finite-size, threshold).

### 5. Λ ≈ M_Pl (gravitational sector)

If Λ determines the Planck mass via f₂Λ²a₂ = (6/π³)M_Pl², then:
- Λ ≈ 1.25 × M_Pl for f₂ = 1
- f₂ = π³/(6a₂) ≈ 1.56 if Λ = M_Pl exactly
- t = (m_W/M_Pl)² ≈ 4.3 × 10⁻³⁵ (NOT related to the heat kernel t ≈ 2 from K₄ geodesic)

The heat kernel parameter t from the K₄ moduli space is an **intrinsic** dimensionless quantity (geodesic distance² on S²). It is not 1/Λ² in any physical units.

---

## Updated Parameter Count

| Parameter | Status | What determines it |
|:---|:---|:---|
| a₂ | Fixed | K₆ Gram matrix eigenvalue |
| a₄ | Fixed | K₆ quartic form |
| m_H | Predicted | √(2a₄/a₂) × m_W, zero free parameters |
| N_gen | Predicted | φ(7)/2 = 3 from K₈ vacuum rank |
| Yukawa hierarchy | Predicted | 415:135:1 from K₈ vacuum eigenvector |
| t (heat kernel) | Fixed | d² from K₄ geodesic flow |
| c (normalization) | Fixed | 4/a₂ from BZ formula structure |
| J_C (complex structure) | Fixed | Ψ₂ from K₄ (forced, unique) |
| **Λ** | **Free** | **Determines M_Pl/m_W hierarchy** |
| **f₂** | **Free** | **Degenerate with Λ: only f₂Λ² matters** |

**Effective free parameter count: 1 (the hierarchy ratio f₂Λ²/M_Pl²).**

But this one parameter does NOT affect any spectral prediction. It affects only the gravitational sector. The matching framework's particle physics predictions are **zero-parameter**.

---

## What This Means

The honest statement is:

> **The matching framework has zero free parameters for particle physics and one free parameter for gravity.** That one parameter is the same free parameter every theory of quantum gravity faces: the Planck scale (equivalently, the gauge hierarchy).

The Higgs mass, generation count, Yukawa hierarchy, gauge group, and parity violation are all combinatorial/topological and independent of Λ. The Planck mass, cosmological constant, and EW/Planck hierarchy require knowing Λ, which the matching framework does not determine from internal structure alone.

This is actually the *cleanest possible* outcome. Λ not being internally determined means it's the **interface** between the matching framework (particle physics) and quantum gravity (Planck scale). The matching framework completes the particle physics side. Gravity remains open.

---

## The Revised Architecture

```
MATCHING FRAMEWORK (zero free parameters)
├── K₆ → a₂, a₄ → m_H = 126.1 GeV
├── K₈ → φ(7)/2 = 3 → three generations
├── K₈ vacuum → 415:135:1 → Yukawa hierarchy
├── K₄ → su(2)₊ → parity violation
├── K₆ × K₄ → SU(3) × SU(2) × U(1) → gauge group
└── K₄ moduli → t, c, J_C → spectral geometry complete

    ↕ interface: f₂Λ² (one free parameter)

GRAVITATIONAL SECTOR (one free parameter)
├── f₂Λ²a₂ → M_Pl → Newton's constant
├── f₄Λ⁴a₀ → Λ_CC → cosmological constant
└── M_Pl/m_W ratio → hierarchy problem
```

The interface is clean. The particle physics is complete. The gravity is not.
