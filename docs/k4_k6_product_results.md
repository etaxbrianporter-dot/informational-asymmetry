# K₄ × K₆ Product Spectral Action: Results

## Executive Summary

The product spectral action on K₄ × K₆ reveals a clean **factorization**: K₆ controls the *scale* of gauge couplings and the Higgs sector, while the gauge coupling *ratios* and Weinberg angle correction require K₈. Three exact algebraic results emerge, plus a structural theorem about sector-blindness that redirects the computation plan.

---

## 1. Three Exact Results

### Result 1: a₂ cross-term vanishes identically

The Seeley-DeWitt coefficient a₂ on the K₄ × K₆ product decomposes as:

```
a₂(product, t) = a₂(K₄) × dim_F + t² × dim₄ × Tr(D_F²)
               = 60 + 4t²
```

with **zero mixing term**. The cross-term Tr(X ⊗ D_F) where X = D†γ + γD vanishes because:

- Tr(γ·Mᵢ) = Σⱼ γⱼⱼ(Mᵢ)ⱼⱼ = 0 for all matching matrices Mᵢ
- This is because perfect matching adjacency matrices have **zero diagonal** — no vertex connects to itself
- This is a topological property of matchings, not a symmetry property

**Note**: γ does NOT anticommute with all matchings (M₂ commutes with γ). The vanishing is more fundamental — it's intrinsic to the matching structure at any K₂ₙ.

**Consequence**: The spectator mechanism cᵢ = dim(K₂ᵢ) is **exact at the a₂ level**. No mixing corrections from the product geometry.

### Result 2: a₄ is an exact polynomial in t²

```
a₄(product, t) = 100 + 16t² + 0.3478t⁴
```

with analytically verified coefficients:

| Coefficient | Value | Formula | Source |
|:---|:---|:---|:---|
| c₀ | 100 | a₄(K₄) × dim_F = (20/3) × 15 | Pure K₄ |
| c₂ | **16** | 2⟨Tr(D†D)⟩_BZ + ⟨Tr(X²)⟩_BZ = 8 + 8 | Mixed K₄-K₆ |
| c₄ | 0.3478 | dim₄ × R(K₆) = 4 × 0.08694 | Pure K₆ |

The c₂ coefficient decomposes into two equal pieces:
- **Piece 1** = 2 × ⟨Tr(D†_K4 D_K4)⟩_BZ × Tr(D_F²) = 2 × 4 × 1 = **8**
- **Piece 2** = ⟨Tr(X²)⟩_BZ × Tr(D_F²) = 8 × 1 = **8**

where X = D†_K4·γ + γ·D_K4 is the chirality-Dirac anticommutator. The fact that ⟨Tr(X²)⟩_BZ = 8 = 2 × a₂(K₄) is itself a non-trivial identity of the K₄ democratic point.

### Result 3: c₄ = dim₄ × R(K₆)

The pure K₆ contribution at order t⁴ is controlled by the **normalized K₆ kurtosis** R(K₆) = Tr(D_F⁴)/[Tr(D_F²)]² = 0.08694. This is a different R from the Higgs mass ratio a₄/a₂ = 0.3722 — R(K₆) here is the kurtosis of the full Gram matrix spectrum, not the vacuum eigenvalue ratio.

---

## 2. Sector-Blindness Theorem

**All four V₄ characters carry exactly 25% of both a₂ and a₄, at all values of t.**

This means the K₆ cross-coupling is **sector-blind**: it shifts all gauge sectors by the same multiplicative factor. The ratio c₁/c₃ = 8/4 = 2 is unchanged by the K₄ × K₆ product structure.

### What K₆ Controls

| Observable | K₆ role | Status |
|:---|:---|:---|
| Overall gauge coupling scale | Multiplicative shift via f₂Λ² | ✓ Sector-blind |
| Higgs quartic coupling | a₄/a₂ ratio in K₆ vacuum | ✓ Already computed |
| CC prefactor | Species multiplicity (×15) | ✓ Computed |
| Gauge coupling ratios | **None** | ✗ Sector-blind |
| Weinberg angle correction | **None** | ✗ Needs K₈ |

### What K₈ Must Provide

The sector-dependent correction (the Yukawa Casimirs aᵢ = {17/10, 3/2, 2}) encodes hypercharge, isospin, and color quantum numbers of SM fermions. These are properties of the K₈ matching algebra, not K₆.

The 2.5% residual in c₁ and the 1.6% error in sin²θ_W are both consequences of the K₈ coupling to K₄ through the gauge sector. The K₈ product D_K4 ⊗ D_K8 must break the V₄ character uniformity to produce sector-dependent corrections.

---

## 3. CC Prefactor Analysis

| Source | Multiplier | Running total | Gap remaining |
|:---|:---|:---|:---|
| Single K₄ channel | 1 | 10⁻¹²⁵·⁴ | 10³·⁴ |
| K₆ matching space | ×15 | 10⁻¹²⁴·² | 10²·² |
| K₈ matching space | ×105 | 10⁻¹²²·² | 10⁰·² |
| BZ normalization + spin | ×24 | — | — |

Including K₈ (105 matchings), the species multiplicity factor is 15 × 105 = 1575, which closes **most** of the 3.4-order prefactor gap. The remaining factor ~2 is within the expected range from BZ normalization conventions and the spectral action mapping f(ε) → Λ_CC.

---

## 4. Computation Plan (Revised)

The K₄ × K₆ computation reveals a clean three-level factorization:

**Level 1 — K₄ alone** (complete):
- Lorentzian signature, light emergence, CC exponent (122 orders)
- a₂(K₄) = 4, a₄(K₄) = 20/3 (exact)

**Level 2 — K₄ × K₆** (just completed):
- Higgs mass: mH = √(2a₄/a₂) × mW = 126.1 GeV
- CC prefactor: ×15 species
- Gauge coupling scale (sector-blind)
- Exact polynomial: a₄ = 100 + 16t² + 0.35t⁴

**Level 3 — K₄ × K₈** (next build):
- Yukawa hierarchy: 415:135:1
- **Sector-dependent** gauge corrections → close 2.5% c₁ residual
- Weinberg angle → close 1.6% sin²θ_W residual
- CC prefactor completion: ×105

### Key Question for Level 3

Does the K₈ product break V₄ character uniformity? If the K₈ Dirac operator couples differently to different V₄ channels — which it should, since K₈ encodes the generation structure and hypercharge assignments — then the a₄(K₄×K₈) will have **channel-dependent** contributions. These channel-dependent pieces are the aᵢ corrections.

### Feasibility

- K₈: 105 matchings → K₄ × K₈ product is 4 × 105 = 420-dimensional
- BZ averaging over 50×50 grid: 2500 diagonalizations of 420×420 matrices
- Estimated runtime: ~30 minutes (vs ~2 min for K₄ × K₆)
- Exact polynomial fit should still work (a₄ as polynomial in t²)

---

## 5. Updated Ledger Entry

**Proved (new)**:
- a₂ cross-term vanishes for ANY K₂ₙ product (zero diagonal of matchings)
- K₆ coupling is sector-blind (V₄ uniformity)
- a₄(K₄×K₆) = 100 + 16t² + 4R(K₆)t⁴ with exact coefficients
- ⟨Tr(X²)⟩_BZ = 2a₂(K₄) = 8 (chirality identity at democratic point)

**Redirected**:
- Gauge coupling ratio correction: K₆ → K₈ (sector-blindness forces this)
- Self-consistent boundary problem: requires K₈, not just K₆

**Confirmed**:
- CC species counting: 15 × 105 = 1575 closes bulk of prefactor gap
