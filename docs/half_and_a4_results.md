# The Factor of ½ and the Seeley-DeWitt a₄ Coefficient

## Part A: The Algebraic Root of ½

### D_F in the SVD Basis

Transforming D_F into the Yukawa SVD basis {u₁,u₂,u₃,u₀,v₁,v₂,v₃,v₀} reveals clean structure:

**Yukawa block (V₊ color → V₋ color):** Exactly diagonal
```
D_pm = diag(0.9377, 0.7656, 0.7438)
```

**Mass blocks:** Non-diagonal, with significant structure
```
||D_cc_pp|| = 0.848  (V₊ color self-coupling)
||D_cc_mm|| = 0.834  (V₋ color self-coupling)
||D_cl_p||  = 0.555  (V₊ color↔lepton mixing)
||D_cl_m||  = 0.688  (V₋ color↔lepton mixing)
||D_ll_pm|| = 0.000  (lepton sector: NO cross-chirality coupling)
```

The lepton null direction has zero Yukawa coupling in both chirality sectors. The lepton is completely decoupled from the Higgs mechanism.

### Block Decomposition of a₃

The total a₃ = Σ_a Tr([D,T_a]²) = 4.170 decomposes as:

| Block | Contribution | Fraction |
|-------|-------------|----------|
| V₊ color↔color | 0.890 | 21.3% |
| V₋ color↔color | 1.138 | 27.3% |
| **Color↔lepton** | **2.008** | **48.1%** |
| Yukawa (cross-chirality) | 0.135 | 3.2% |

**The dominant contribution (48%) comes from color↔lepton mixing**, not from the pure color sector. The SU(3) generators, while acting only on color indices, produce fluctuations that leak into the lepton sector because the SVD "color" directions are not pure Z₇ doublet states (u₃ has 42% singlet content).

### The Yukawa ½: Exact Derivation

For diagonal Yukawa D_pm = diag(σ₁, σ₂, σ₃):

```
[D_pm, T_a]_ij = (σ_i - σ_j) · (T_a)_ij
```

Therefore:
```
Σ_a ||[D_pm, T_a]||² = Σ_{i≠j} (σ_i - σ_j)² × Σ_a |(T_a)_ij|²
```

The off-diagonal Casimir sum for SU(3) fundamental gives:
```
Σ_a |(T_a)_ij|² = ½    (for i ≠ j)
```

So the Yukawa contribution to a₃ is:
```
a₃(Yukawa) = ½ × Σ_{i≠j} (σ_i - σ_j)² = ½ × [6Σσ² - 2(Σσ)²]
```

**The ½ is the off-diagonal SU(3) Casimir coefficient.** It's not a mysterious correction — it's a standard representation-theoretic identity applied to the diagonal Yukawa structure.

### Critical Insight: Gauge Coupling Measures Variance

When all Yukawa eigenvalues are equal (σ₁ = σ₂ = σ₃ = σ):
- [diag(σ,σ,σ), T_a] = 0 for all T_a
- **a₃(Yukawa) = 0 exactly**

The Yukawa contribution to gauge coupling is proportional to the **variance** of the Yukawa spectrum, not its magnitude. More precisely:

```
a₃(Yukawa) = ½ × (2N·Σσ² - 2(Σσ)²) where N = 3
           = ½ × 6σ² × [1 - (Σσ)²/(NΣσ²)]
```

Our 10.6% Yukawa spread gives (Σσ)²/Σσ² = 2.966 (close to N=3 degenerate limit), making a₃(Yukawa) = 0.068 — only 1.6% of total a₃. The gauge coupling is overwhelmingly determined by the mass blocks D₊₊, D₋₋ and the color-lepton mixing, NOT by the Yukawa sector.

### Resolution of the ½ Puzzle

The original puzzle: a₃(SVD) / [C₂ × Tr(P_color·D²)] = 0.502.

This is NOT a simple factor of ½. The full decomposition shows:
- 48.5% of a₃ comes from within-chirality color blocks
- 48.1% from color-lepton mixing  
- 3.2% from Yukawa (cross-chirality)

The Casimir formula C₂ × Tr(P_color·D²) assumes D_F acts entirely within the color sector. But color-lepton mixing redistributes half the spectral weight, giving the apparent ½.

---

## Part B: Seeley-DeWitt a₄ Coefficient

### Raw Numbers

| Quantity | Single copy (ℝ⁸) | Doubled (ℝ¹⁶) |
|----------|------------------|----------------|
| a₀ = Tr(I) | 8 | 16 |
| a₂ = Tr(D²) | 9.247 | 18.495 |
| a₄ = Tr(D⁴) | 16.081 | 32.161 |
| a₄/a₂ | 1.739 | 1.739 |
| a₄/a₂² | 0.188 | 0.094 |

### Z₇ Sector Decomposition

| Sector | a₂ | a₄ | a₄/a₂ |
|--------|-----|-----|--------|
| d₁ | 2.677 | 4.912 | 1.835 |
| d₂ | 2.649 | 4.569 | 1.725 |
| d₃ | 2.600 | 4.239 | 1.630 |
| triv | 0.165 | 0.295 | **1.787** |
| hub | 1.156 | 2.066 | **1.787** |

**Discovery: triv and hub have identical a₄/a₂ = 1.787.** The singlet sector is spectrally uniform despite triv being distributed across all vertices and hub being localized on v=7. This is a non-trivial consistency check — the SU(2) singlet sector has a well-defined spectral ratio.

The three color doublets have a₄/a₂ spanning [1.630, 1.835], giving another measure of the color-sector splitting.

### Chirality Decomposition

a₄(γ=+1) = 8.147 (50.7%)
a₄(γ=-1) = 7.933 (49.3%)

Nearly symmetric — the 1.4% asymmetry is consistent with the mild Yukawa hierarchy.

### Yukawa Sector (Higgs Content)

| Quantity | Value | Fraction of total |
|----------|-------|------------------|
| Σσ² (Higgs mass²) | 2.019 | 21.8% of a₂ |
| Σσ⁴ (Higgs quartic) | 1.423 | 8.8% of a₄ |
| Σσ⁴/Σσ² | 0.705 | — |
| (Σσ²)²/Σσ⁴ (flatness) | 2.864 | — |

The flatness ratio 2.864 (close to the degenerate maximum of 3) reflects the near-degeneracy of the three Yukawa eigenvalues. This controls the Higgs quartic coupling in the spectral action.

### Connection to K₆ Higgs

The K₆ paper found a₄/a₂ = 0.3722 for the Higgs sector. Our K₈ Yukawa gives:
- Σσ⁴/Σσ² = 0.705 (Yukawa-only ratio)
- Full a₄/a₂ = 1.739 (includes all sectors)

The K₈ Yukawa ratio is ~2× the K₆ value. This factor of 2 may arise from the K₈ having two chirality sectors contributing independently, versus K₆'s single Higgs doublet.

---

## Synthesis: Duality Between Gauge Coupling and Yukawa Hierarchy

Both the gauge coupling traces and the Higgs potential are controlled by the same spectral parameter: the **splitting pattern of Yukawa singular values**.

| Limit | σ₁=σ₂=σ₃ (degenerate) | σ₁ ≫ σ₂,σ₃ (hierarchical) |
|-------|----------------------|--------------------------|
| SU(3) symmetry | Exact | Broken |
| a₃(Yukawa) | 0 | Maximal |
| Flatness (Σσ²)²/Σσ⁴ | 3 (flat) | 1 (peaked) |
| Higgs quartic | Minimal | Maximal |
| Physical regime | Our K₈ (r=2.97) | SM top quark (r≈1) |

Our K₈ sits very close to the degenerate limit (r = 2.97 vs maximum 3), which is why:
- SU(3) is nearly exact (10.6% spread)
- The Higgs quartic is nearly minimal
- The Yukawa contribution to gauge coupling is only 3.2%

In the SM, the top quark dominance pushes r toward 1, breaking SU(3) flavor symmetry maximally. The K₈ vertex-space geometry represents the UV (near-degenerate) starting point; the matching-space hierarchy 415:135:1 represents the IR (hierarchical) endpoint. The RG flow from one to the other is the physical Yukawa running.

---

## Ledger Update

**CONFIRMED:**
- ½ traced to off-diagonal SU(3) Casimir: Σ_a |T^a_ij|² = ½ for i≠j
- a₃(Yukawa) ∝ variance of {σ_i}, vanishes for degenerate spectrum
- Dominant a₃ contribution (48%) from color-lepton mixing, not pure color
- Singlet sector (triv, hub) has identical a₄/a₂ = 1.787
- Chirality decomposition of a₄ nearly symmetric (50.7%/49.3%)

**KEY NUMBERS:**
- a₂ = 9.247, a₄ = 16.081, a₄/a₂ = 1.739, a₄/a₂² = 0.188
- Yukawa: Σσ² = 2.019 (21.8% of a₂), Σσ⁴ = 1.423 (8.8% of a₄)
- Flatness ratio: 2.864 (near-degenerate)
- D_ll_pm = 0.000 (lepton completely decoupled from Higgs)

**OPEN:**
- Factor of 2 between K₈ Yukawa ratio and K₆ Higgs ratio
- How does the color-lepton mixing (48% of a₃) relate to quark-lepton universality?
- RG flow connecting vertex-space flatness (r=2.97) to matching-space hierarchy (r≈1)
- Proper extraction of sin²θ_W accounting for the mixed chirality structure of SU(2)
