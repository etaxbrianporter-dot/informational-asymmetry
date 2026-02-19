---
title: The 2/3 Trace: Verdict
section: Killed Approaches
status: active
---

# The 2/3 Trace: Verdict

## ⟨E⟩ Is NOT 10/11

The BZ average converges to ⟨|d₁(k)|⟩ = 0.90909414... which is NOT
10/11 = 0.90909090... The difference 3.23 × 10⁻⁶ is stable from N=400
through N=1600 (ratio of successive differences → 1.0003). It's a
transcendental number — an elliptic integral over the honeycomb BZ.

Consequently:
```
  Var(E)  = 1 - ⟨E⟩²  = 0.173548    (NOT 21/121 = 0.173554)
  Var(E²) = ⟨E⁴⟩ - 1  = 2/3          (EXACT, combinatorial)
```

## What the 2/3 Does

The 2/3 controls the ±E case (no arrow of time):

```
  f₋(ε) = Var(E²) / (4ε² ln ε) = (2/3) / (4ε⁴ ln ε)
  At Hubble: log₁₀(f₋) = -246.6
```

The arrow of time promotes variance → mean:

```
  f₊(ε) = Var(E) / (4ε² ln ε) = 0.17355 / (4ε² ln ε)
  At Hubble: log₁₀(f₊) = -125.4
```

The SHIFT between them:

```
  Δ(orders) = log₁₀(f₋) - log₁₀(f₊)
            = 2 × log₁₀(ε_H) + log₁₀(Var(E²)/Var(E))
            = 2 × 60.93 + log₁₀(3.841)
            = 121.86 + 0.58
            = 122.4 orders
```

## The Three Numbers

| Quantity | Value | Source | Role |
|----------|-------|--------|------|
| Var(E²) = 2/3 | Exact | K₄ combinatorial | Sets ±E baseline |
| Var(E) = 0.17355 | Transcendental | Elliptic integral on BZ | Sets +E magnitude |
| Var(E²)/Var(E) = 3.841 | Transcendental | Ratio | Controls the shift |

The shift of 122.4 orders comes from 2 × log₁₀(ε_H) + 0.58:

- The 121.86 part is PURE DIMENSIONAL ANALYSIS (2 powers of ε)
- The 0.58 part is log₁₀(3.841) — a combination of 2/3 and ⟨E⟩

So the 2/3 contributes ~0.4 orders out of 122. The overwhelming
majority comes from the dimensional shift ε⁴ → ε².

## Absolute Magnitude: Off by ~2 Orders

| Quantity | log₁₀ |
|----------|-------|
| f₊(ε_H), Gaussian channel | **-125.4** |
| ρ_Λ/M_Pl⁴ (energy density) | -123.2 |
| Λ/M_Pl² = 3H₀² | -121.4 |
| (H₀/M_Pl)² | -121.9 |

The Gaussian channel formula overshoots by 2-4 orders depending
on which convention for Λ_CC we compare against.

The discrepancy is in the PREFACTOR, not the scaling. Possible
sources:
- The factor of 4 in the denominator (2D BZ normalization)
- 6 Dirac points (each cone contributes)
- The mapping f(ε) → Λ_CC involves the spectral action, not
  just information theory
- The precise physical meaning of "resolution ε"

## What IS Precise

The arrow-of-time shift: **122.4 orders**

This is robust to all normalization conventions because it's a RATIO.
Whatever prefactor C multiplies f₊, the same C multiplies f₋.
The ratio f₋/f₊ = (Var(E²)/Var(E)) × ε² is convention-independent.

Evaluated: Δ = 2 × 60.93 + log₁₀(3.841) = 122.4

Observed: the CC is ~122 orders below M_Pl⁴.

This is a zero-parameter result. The 122 comes from:
- The Hubble/Planck ratio (60.93 from cosmology)  
- The combinatorial ratio Var(E²)/Var(E) = 3.841 from K₄
- The factor of 2 from ε⁴ → ε² (the arrow of time)
