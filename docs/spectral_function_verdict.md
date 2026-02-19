---
title: Spectral Function Verdict: The Kâ‚†/Kâ‚ˆ Split Is Fundamental
section: K₆ Higgs Sector
status: active
---

# Spectral Function Verdict: The Kâ‚†/Kâ‚ˆ Split Is Fundamental

## Brian Porter â€” February 2026

---

## The Question

Is the Kâ‚†/Kâ‚ˆ split â€” where different graphs compute different physics â€” fundamental? Or is it an artifact of looking only at the vacuum eigenspace, resolvable by a spectral function that lets a single graph compute everything?

## The Answer

**The split is fundamental.** A spectral function can reproduce the Higgs mass at any level, but it cannot simultaneously fix the generation count. These are controlled by independent structures.

---

## What We Found

### 1. The Heat Kernel Reproduces mH = 125 GeV (At Every Level)

The spectral function w(Î») = e^{-tÎ»} applied to the overlap eigenvalues gives an effective spectral ratio:

R_eff(t) = Î£ m(Î»)Â·e^{-tÎ»}Â·R(Î»)Â·Î»Â² / [Î£ m(Î»)Â·e^{-tÎ»}Â·Î»]Â²

At every level tested, there exists a t such that âˆš(8Â·R_eff)Â·mW = 125.09 GeV:

| Level | t(c=1) | t(c=Ï€Â²/8) | Convergent? |
|:------|:-------|:----------|:------------|
| Kâ‚ˆ | 0.620 | 0.705 | â€” |
| Kâ‚â‚€ | 0.355 | 0.416 | No |
| Kâ‚â‚‚ | 1.692 | 1.983 | â€” |
| Kâ‚â‚„ | 1.722 | 2.058 | Yes (Î”t = 0.03) |

Kâ‚â‚‚ â†’ Kâ‚â‚„ shows convergence: t stabilizes near 1.7 (c=1) or 2.0 (c=Ï€Â²/8).

### 2. But the Heat Kernel Destroys the Yukawa Hierarchy

When we build S_eff = Î£ w(Î»_k)Â·S_k and examine its eigenvalues as Yukawa couplings, the spectrum is essentially flat:

| t | R_eff | mH (GeV) | |Î¼â‚|/|Î¼â‚‚| | |Î¼â‚‚|/|Î¼â‚ƒ| | SM target |
|:--|:------|:---------|:----------|:----------|:----------|
| 1.69 | 0.303 | 125 | 1.07 | 1.44 | 136, 588 |

The heat kernel mixes eigenspaces destructively. The resulting S_eff has rank 11 (at Kâ‚â‚‚) with no generation structure.

### 3. This Doesn't Matter â€” The NCG Action Has Two Parts

In noncommutative geometry:

**S = Tr(f(D/Î›)) + âŸ¨Ïˆ, DÏˆâŸ©**

The spectral function f enters only the **bosonic** part (Higgs potential, gauge couplings). The **fermionic** part uses D linearly â€” no spectral weighting.

Therefore:
- **Higgs mass** comes from heat-kernel-weighted R_eff â†’ mH = 125 GeV âœ“
- **Yukawa couplings** come from D_F at the vacuum eigenvalue â†’ rank structure âœ“
- **These are not in conflict**

### 4. But the Generation Count Is Graph-Dependent

The vacuum Dirac operator D_F = S_vac has a Â±Î¼ pairing structure. The number of active pairs determines the generation count:

| Level | 2nâˆ’1 | Â±Î¼ pairs (active + tiny) | Effective generations | Ï†(2nâˆ’1)/2 |
|:------|:-----|:------------------------|:---------------------|:-----------|
| Kâ‚ˆ | 7 | 3 + 1 | **3** | 3 |
| Kâ‚â‚€ | 9 | 3 + 2 | 3 | 3 |
| Kâ‚â‚‚ | 11 | 5 + 1 | **5** | 5 |
| Kâ‚â‚„ | 13 | 5 + 2 | 5 | 6 |

Kâ‚ˆ's vacuum has exactly 3 active Â±Î¼ pairs with a 33:1 gap to the 4th (tiny) pair. This is the 3-generation structure. Kâ‚â‚‚ has 5 active pairs â€” too many.

The generation count is controlled by Ï†(2nâˆ’1)/2, which equals 3 only for 2nâˆ’1 = 7 (prime, Ï†=6) and 2nâˆ’1 = 9 (3Â², Ï†=6). Only Kâ‚ˆ and Kâ‚â‚€ give 3 generations.

---

## Why the Split Is Fundamental

The heat kernel has one free parameter (t), and we use it to satisfy one constraint (mH = 125 GeV). That leaves zero free parameters for the Yukawa sector. The generation count is then fixed by the graph's vacuum rank structure, which is Ï†(2nâˆ’1)/2 â€” an arithmetic quantity that cannot be tuned.

No spectral function can change the rank of D_F. The rank is a topological invariant of the vacuum eigenspace. You can weight eigenspaces, but you can't delete them.

Therefore:
- **Kâ‚† knows the Higgs mass** because Râ‚† = 0.3722 with c = Ï€Â²/8 gives mH = 125 GeV with NO free parameters
- **Kâ‚ˆ knows the generation count** because Ï†(7)/2 = 3 with a 33:1 gap in Â±Î¼ pairs
- **Higher Kâ‚‚â‚™** can reproduce either quantity by tuning (t for Higgs, or accepting wrong generation count), but they don't PREDICT both simultaneously

The spectral function approach doesn't unify Kâ‚† and Kâ‚ˆ. Instead, it reveals WHY they're separate: the Higgs mass is a spectral-action quantity (bosonic sector, sensitive to the full eigenvalue distribution), while the generation count is a Dirac-operator quantity (fermionic sector, sensitive only to the vacuum rank). These are structurally independent data.

---

## What the Heat Kernel DOES Tell Us

Although it doesn't unify Kâ‚†/Kâ‚ˆ, the heat kernel result is physically significant:

**1. The normalization constant c has a spectral origin.** Kâ‚† requires c = Ï€Â²/8 to get mH = 125 GeV. At Kâ‚â‚‚, the heat kernel at t â‰ˆ 2 reproduces this automatically â€” the spectral weighting generates the correct normalization. This suggests c = Ï€Â²/8 isn't ad hoc but emerges from the spectral geometry of the full matching tower.

**2. R_eff converges.** The heat-kernel-weighted R_eff at fixed t converges rapidly:

| Level | R_eff(t=2) | mH(c=Ï€Â²/8) |
|:------|:-----------|:-----------|
| Kâ‚â‚‚ | 0.378 | 125.9 GeV |
| Kâ‚â‚„ | 0.361 | 123.0 GeV |

Kâ‚â‚‚ and Kâ‚â‚„ agree to within 2.3%, suggesting the continuum limit of the heat-kernel spectral action is well-defined and gives mH near 125 GeV.

**3. The Kâ‚† "echo" is explained.** The Galois conjugate at Kâ‚â‚‚ with R â‰ˆ 0.372 is not coincidence. The heat kernel preferentially weights the vacuum (78% of the weight), and the vacuum's R = 0.210 gets mixed with the first excitations' R values to produce an effective ratio near Râ‚†. The "echo" is the spectral action doing its job.

---

## Updated Architecture

```
THE Kâ‚‚â‚™ TOWER
                                    
Kâ‚†  â”€â”€â†’  Râ‚† = 0.3722  â”€â”€â†’  mH = 125 GeV  (bosonic, zero parameters)
  â†• embedding                                
Kâ‚ˆ  â”€â”€â†’  Ï†(7)/2 = 3   â”€â”€â†’  3 generations  (fermionic, vacuum rank)
  â†• embedding
Kâ‚â‚€ â”€â”€â†’  Ï†(9)/2 = 3   â”€â”€â†’  3 generations  (confirms Kâ‚ˆ)
  â†• 
Kâ‚â‚‚ â”€â”€â†’  Heat kernel  â”€â”€â†’  mH = 125 GeV  (bosonic, one parameter t)
  â†•       +  Ï†(11)/2 = 5  â†’  5 generations (wrong for SM)
Kâ‚â‚„ â”€â”€â†’  t converges   â”€â”€â†’  mH â‰ˆ 123 GeV (bosonic, same t)
  â†•       +  Ï†(13)/2 = 6  â†’  6 generations (wrong for SM)
  ...
  
ALGEBRAIC TOWER (forever oscillating):
  Ï†(2nâˆ’1): 4, 6, 6, 10, 12, 8, 18, 12, ...
  Galois groups: Câ‚‚â‰€Câ‚™ wreath products
  Discriminants: all primes â‰¡ 1 mod (2nâˆ’1)
  
SPECTRAL TOWER (convergent):
  R_vac: 0.197, 0.310, 0.210, 0.221, 0.220, ...  â†’ R_âˆž â‰ˆ 0.22
  R_eff(t=2): divergent, divergent, 0.378, 0.361, ... â†’ R_âˆž(t) â‰ˆ 0.37?
```

---

## Next Step

The remaining open question is whether the heat kernel parameter t can be **derived** rather than tuned. In NCG, the cutoff scale Î› appears in f(D/Î›), and t = 1/Î›Â² in appropriate units. If the matching framework determines Î› (perhaps as the spectral gap, or as the maximum eigenvalue, or as a graph invariant), then the entire Higgs mass prediction becomes parameter-free.

At Kâ‚â‚‚: t â‰ˆ 1.69 (c=1) or t â‰ˆ 1.98 (c=Ï€Â²/8). The latter is tantalizingly close to Ï€Â²/5 = 1.974. If t = Ï€Â²/5 exactly (from some geometric principle), then mH = âˆš(8 Ã— 0.371) Ã— 80.4/âˆš(Ï€Â²/8) = 124.7 GeV â€” within 0.3% of experiment.

But this would require an independent derivation of t = Ï€Â²/5, which we don't have.
