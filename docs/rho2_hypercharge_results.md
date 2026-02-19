---
title: ρ₂ Hypercharge: Results
section: 
status: active
---

# ρ₂ Hypercharge: Results

## The ℤ₇ Decomposition is Perfect

The 6-fold vacuum decomposes as ρ₁ ⊕ ρ₂ ⊕ ρ₃ with exactly equal weight:
```
ρ₁ (k=1,6): weight = 1/3
ρ₂ (k=2,5): weight = 1/3  
ρ₃ (k=3,4): weight = 1/3
```
All six ℤ₇ eigenvalues are unit-magnitude roots of unity. Democratic.

## KILL: ρ₂ Is NOT an SU(3)×SU(2) Singlet

Normalized gauge kinetic coefficients per doublet:

| Sector | c_Y | c_su3 | c_su2 | c_Y/c_su3 |
|--------|-----|-------|-------|-----------|
| ρ₁ | 0.359 | 0.412 | 0.108 | 0.87 |
| ρ₂ | 0.257 | 0.262 | 0.124 | 0.98 |
| ρ₃ | 0.382 | 0.359 | 0.087 | 1.06 |

**All three doublets carry comparable SU(3), SU(2), and U(1) charges.**
The ratio c_Y/c_su3 ≈ 1 for all sectors — there is no SU(3)×SU(2)-singlet
sector with enhanced hypercharge.

This kills the hypothesis that ρ₂ alone fills the c₁ gap while keeping c₃ = 4.

## WHY: Y Breaks ℤ₇

The hypercharge operator Y = diag(+1/2, -1/2, +1/3, -1/3, +1/3, -1/3, +1/3, -1/3)
is NOT ℤ₇-symmetric. Hub vertices carry Y = ±1/2, K₆ vertices carry Y = ±1/3.

The ℤ₇-averaged Y is proportional to identity on vertices 0-6:
  Y_avg(0,...,6) = 1/21 ≈ 0.0476

The ℤ₇-symmetric part of C_Y contributes EQUALLY to all three doublets:
```
C_Y(sym): ρ₁ = ρ₂ = ρ₃ = 0.151  (exact equality)
```

The ℤ₇-BREAKING part distributes unequally but does NOT single out ρ₂:
```
C_Y(break): ρ₁ = 1.254, ρ₂ = 0.856, ρ₃ = 1.346
```

## Key Structural Insight

The double commutator measures gauge-charge "disruption" of matchings.
Since Y acts on VERTICES while ℤ₇ permutes vertices cyclically,
and since the hub vertices (0,1) are NOT at a fixed point of ℤ₇,
the Y charge cannot align with the ℤ₇ decomposition.

In other words: the physical gauge quantum numbers (Y, color, isospin)
and the generation quantum numbers (ℤ₇ irreps) live in ORTHOGONAL
structures. No ρ_k can be a pure gauge singlet because gauge charges
aren't ℤ₇-symmetric.

## What Survives

The double commutator traces on the full vacuum:
```
Tr(C_Y|vac)   = 3.909  (= 6 × 0.652)
Tr(C_su3|vac) = 4.047  (= 6 × 0.674)  
Tr(C_su2|vac) = 1.251  (= 6 × 0.209)
```

Ratios: C_Y : C_su3/8 : C_su2/3 = 0.652 : 0.084 : 0.070
Per-generator: U(1) contributes 8-9× more than non-abelian per generator.

This enhanced U(1) response is real but doesn't solve the problem because
the non-abelian traces are also nonzero and would shift c₃ and c₂.

## Status of c₁ = 8

| Hypothesis | Status |
|-----------|--------|
| A. K₈ ℂ⁸ added to ℋ_F | KILLED (breaks c₃) |
| B. Double commutator on vacuum | KILLED (gives zero vac→vac) |
| B'. Nonzero C_Y projection exists | TRUE but also shifts c₃, c₂ |
| C. ρ₂ is SU(3)×SU(2)-singlet | KILLED (all ρ_k carry all charges) |
| D. Vertex-space Tr(Y²D²) correlation | CHECKED: ~zero correlation |

## The Remaining Path

The computation shows that within matching space, gauge charges are
democratic across ℤ₇ sectors. The c₁ = 8 prediction cannot come from
K₈ matching space alone.

This leaves two possibilities:
1. c₁ = 8 comes from a mechanism OUTSIDE the matching algebra
   (e.g., anomaly cancellation, index theorem, or topological constraint)
2. c₁ = 26/5 is the tree-level value, and c₁ = 8 is an effective value
   including loop corrections (threshold effects at Λ ~ 10⁸ GeV)
3. The spectral triple construction couples K₈ to K₆×K₄ through a
   mechanism we haven't identified (the off-diagonal D_F block)
