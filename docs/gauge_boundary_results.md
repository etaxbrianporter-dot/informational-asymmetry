# Gauge Coupling Self-Consistent Boundary: Results

## Brian Porter — February 20, 2026

---

## The Problem

The spectator mechanism predicts c₃:c₂:c₁ = 4:6:8 at the cutoff Λ. After 2-loop RG running to M_Z, this gives sin²θ_W = 0.2349 (1.6% high) and c₁_eff = 8.20 (2.5% high). Can the next-order Seeley-DeWitt term (f₀) close the gap?

## The Two-Parameter System

Boundary conditions at Λ:

```
1/α_i(Λ) = K₁ · c_i + K₂ · c'_i
```

where K₁ = f₂Λ²/(2π²), K₂ = f₀/(4π²), and:

| Gauge factor | c_i (leading) | c'_i (subleading) | c'_i/c_i |
|:---|:---:|:---:|:---:|
| U(1)_Y | 8 | 17/10 | 0.213 |
| SU(2) | 6 | 3/2 | 0.250 |
| SU(3) | 4 | 2 | 0.500 |

Three observables (α₁, α₂, α₃ at M_Z), three unknowns (K₁, K₂, Λ).

## Solution

| Parameter | Value |
|:---|:---|
| Λ | 3.93 × 10⁸ GeV (log₁₀ = 8.59) |
| K₁ = f₂Λ²/(2π²) | 5.784 |
| K₂ = f₀/(4π²) | 1.319 |
| ρ = f₀/(f₂Λ²) | **0.228** |
| f₂Λ² | 114.2 |
| f₀ | 52.1 |
| y_t²(Λ) | 0.672 |

Boundary conditions at Λ:

| Coupling | 1/α(Λ) | Leading (f₂) | Subleading (f₀) | f₀ fraction |
|:---|:---:|:---:|:---:|:---:|
| U(1) | 48.52 | 46.28 (95.4%) | 2.24 (4.6%) | 4.8% |
| SU(2) | 36.68 | 34.71 (94.6%) | 1.98 (5.4%) | 5.7% |
| SU(3) | 25.78 | 23.14 (89.8%) | 2.64 (10.2%) | 11.4% |

The f₀ correction is 5–11% of the leading term, largest for SU(3).

## Reference: f₀ = 0

| Parameter | f₀ = 0 | f₀ ≠ 0 |
|:---|:---:|:---:|
| Λ | 7.6 × 10⁷ GeV | 3.9 × 10⁸ GeV |
| c₁_eff | 8.29 (+3.7%) | 8.00 (exact) |
| sin²θ_W | 0.2349 (+1.6%) | 0.23122 (exact) |

---

## CRITICAL DISCOVERY: Higgs Mass Consistency Test

The two-parameter solution requires ρ = f₀/(f₂Λ²) = 0.228. This has a sharp consequence for the Higgs sector:

**In standard NCG:**
```
λ/g² = ρ × (a₄/a₂)
m_H² = 2ρ(a₄/a₂)m_W²
```
With ρ = 0.228: **m_H = 60.2 GeV** ← catastrophically wrong

**In the BZ framework:**
```
m_H = √(2 a₄/a₂) × m_W = 126.1 GeV  (ρ-independent)
```

| | ρ required | m_H predicted |
|:---|:---:|:---:|
| Gauge couplings | 0.228 | — |
| Standard NCG Higgs | 0.984 | 60.2 GeV (with gauge ρ) |
| BZ framework Higgs | drops out | 126.1 GeV |

**The tension factor is 4.3×.** Standard NCG needs ρ ≈ 1 for the Higgs mass but the gauge couplings need ρ ≈ 0.23. These are incompatible.

**The BZ framework resolves this completely:** the Higgs mass comes from the K₆ BZ average where a₄/a₂ is a pure combinatorial ratio and f₀/(f₂Λ²) cancels. The gauge couplings see f₀ through the spectator mechanism on the product geometry K₄ × K₆ × K₈.

This is a strong structural argument FOR the BZ interpretation over standard NCG.

---

## Honest Assessment

### What was established

1. The two-parameter system has a **unique physical solution** (K₁ > 0, K₂ > 0, Λ > 0)
2. ρ = 0.228 — a 23% subleading correction (significant but sub-unity)
3. Λ shifts from ~10⁷·⁹ to ~10⁸·⁶ (same ballpark, slightly higher)
4. The Higgs mass consistency **forces** the BZ interpretation

### What was NOT established

1. **No new prediction.** 3 equations, 3 unknowns → always solvable. The system is exactly determined — no overconstrained test.
2. **Higher orders not controlled.** ρ² ≈ 0.052 → the f₋₂/Λ² term could contribute ~5%, comparable to what we're fitting.
3. **Threshold corrections not included.** Step-function decoupling at M_t, M_W, M_H changes the effective β coefficients. Expected O(1%) effect.

### For the ledger

| Claim | Status | Update |
|:---|:---|:---|
| f₀ closes 2.5% gap | **CONSISTENT** | ρ = 0.228 gives exact fit |
| Expansion convergence | **MARGINAL** | ρ = 0.23, ρ² = 0.05 |
| BZ framework required | **STRONGLY SUPPORTED** | Standard NCG gives m_H = 60 GeV with gauge ρ |
| Predictive power | **NONE (this step)** | Need 4th observable |

---

## Next Steps (ranked by value)

1. **Threshold corrections.** Add step-function decoupling at M_t = 173 GeV and M_W = 80.4 GeV. This introduces NO new parameters but shifts ρ and Λ. If ρ decreases → leading-order spectator is more accurate. If ρ increases → expansion diverging.

2. **5/3 GUT normalization.** The factor 5/3 in α₁ comes from the SU(5) embedding of U(1)_Y. Does K₈ modify this? The representation content of K₈ is not necessarily SU(5)-compatible. A different normalization k/3 directly shifts c₁. Even a 3.6% change in the normalization closes the gap without f₀.

3. **Fourth observable.** The top Yukawa y_t²(Λ) = 0.672 is a prediction from the spectral triple on K₈. Comparing this with the K₈ matching algebra prediction would overconstrain the system and provide a genuine test.
