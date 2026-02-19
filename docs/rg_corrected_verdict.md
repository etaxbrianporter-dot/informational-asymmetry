---
title: RG Running: Corrected Verdict
section: Killed Approaches
status: active
---

# RG Running: Corrected Verdict

## Bug Fix

Previous run had d(g²)/dt = bg⁴/(16π²) — should be bg⁴/(8π²).
The factor-of-2 error made couplings run at half speed, shifting
the crossing scale from 10⁸ to 10¹⁴ and making c₁ appear closer to 8.
Now verified: numerical 1-loop matches analytic to 4 significant figures.

## Results

### 1-loop (analytic, verified numerically)
```
At c₂/c₃ = 3/2: c₁ = 8.069,  Λ = 1.73 × 10⁸ GeV  (0.87% from 8)
At c₁ = 8:       c₂ = 5.972,  Λ = 2.01 × 10⁸ GeV  (0.47% from 6)
```

### 2-loop (with top Yukawa)
```
At c₂/c₃ = 3/2: c₁ = 8.199,  Λ = 1.03 × 10⁸ GeV  (2.5% from 8)
At c₁ = 8:       c₂ = 5.920,  Λ = 1.57 × 10⁸ GeV  (1.3% from 6)
```

### Forward prediction (c₁:c₂:c₃ = 8:6:4 exact at Λ = 1.03 × 10⁸)
```
α₂(M_Z) = 0.033801  (obs: 0.033801,  Δ = -0.0004%)  ← exact by construction
α₃(M_Z) = 0.117999  (obs: 0.118000,  Δ = -0.001%)   ← near-exact
α₁(M_Z) = 0.017298  (obs: 0.016943,  Δ = +2.1%)     ← 2.1% high
sin²θ_W = 0.23492   (obs: 0.23122,   Δ = +1.6%)      ← 1.6% high
```

## Honest Assessment

**The 2-loop corrections move c₁ AWAY from 8, not toward it.**

1-loop: c₁ = 8.07 (0.87% off)
2-loop: c₁ = 8.20 (2.5% off)

The dominant 2-loop effect is the top Yukawa, which modifies the
U(1) running more than SU(2) or SU(3). This pushes 1/α₁ up relative
to 1/α₃, increasing the effective c₁ at the matching scale.

## What's Still Missing

- **Threshold corrections:** Step-function decoupling at M_t, M_W, M_H
  changes the effective beta coefficients between thresholds. These are
  O(α) corrections, potentially ~1-2%.

- **GUT normalization factor:** We assumed the standard 5/3 for U(1).
  The framework might modify this if the K₈ representation content
  differs from minimal SU(5) embedding. This would directly shift c₁.

- **Spectral action higher moments:** The relation 1/α_i = K c_i assumes
  only the leading f₂ Λ² term. The next term (f₀, scale-independent)
  adds c_i-dependent corrections.

## What IS Established

1. **The pattern c_i = dim(K_{2i}) works to 2.5%.** This is a
   zero-parameter prediction from pure combinatorics. The framework
   says c₃:c₂:c₁ = 4:6:8 and nature says ≈ 4:6:8.2.

2. **The cutoff scale is ~10⁸ GeV.** NOT the GUT scale 10¹⁶.
   This is qualitatively different from standard NCG.

3. **Non-unification is predicted.** The three couplings never meet.
   They are distinct at every scale, with ratios approaching 4:6:8
   only near 10⁸ GeV.

4. **sin²θ_W is predicted to ~1.6%.** Zero parameters. The observed
   value 0.23122 vs predicted 0.23492 is a meaningful test.

5. **The 2.5% residual is real.** It's not numerical noise, not a
   normalization ambiguity, not an input uncertainty. It requires
   either threshold corrections, modified GUT normalization, or
   spectral action higher-order terms to close.

## For the Ledger

| Claim | Status |
|-------|--------|
| c₁ = 8 = dim(K₈) | **2.5% (2-loop), direction: too high** |
| sin²θ_W = 0.23122 | **1.6% (zero parameters)** |
| Λ ~ 10⁸ GeV | **PREDICTION (not GUT scale)** |
| Non-unification | **PREDICTED** |
| 2-loop closes gap | **KILLED (moves wrong way)** |
