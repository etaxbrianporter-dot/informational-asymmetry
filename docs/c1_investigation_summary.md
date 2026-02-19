---
title: câ‚ = 8 Investigation: Results and Status
section: Program Overview
status: active
---

# câ‚ = 8 Investigation: Results and Status

## Three Major Discoveries

### 1. The Kâ‚ˆ vacuum IS a gauge singlet
The single commutator test:
```
||A_Y @ v_k||   = 0.000000  (all 6 vacuum directions)
||A_su3 @ v_k|| = 0.000000  (all 6 vacuum directions)
||A_su2 @ v_k|| = 0.000000  (all 6 vacuum directions)
```
The gauge transform [T, D_vac] is **orthogonal to matching space** â€” it doesn't vanish in vertex space (||[Y, D_eff]|| â‰ˆ 0.82), but its projection onto the matching algebra is zero.

**Physical meaning:** The vacuum superposition Î£ cáµ¢ Máµ¢ is gauge-invariant *within the matching algebra*, even though individual matchings are not.

### 2. Double commutator gives nonzero vacuum projection
Despite gauge singlet property, the *double* commutator projects back:
```
Tr(C_Y|vac)   = 3.909   (normalized: 0.332 per direction)
Tr(C_su3|vac) = 4.047   (normalized: 0.043 per SU(3) gen per direction)
Tr(C_su2|vac) = 1.251   (normalized: 0.035 per SU(2) gen per direction)
```
U(1) is 8-10Ã— stronger than non-abelian traces per generator on vacuum.

### 3. The coupling factor is exactly dim(Kâ‚„) = 4
For câ‚_GUT = 8 from the vertex trace Tr(YÂ²) = 7/6 on â„‚â¸:
```
Î”câ‚_SM needed = 14/3
Tr(YÂ²|â„‚â¸)    = 7/6
Factor needed = (14/3) / (7/6) = 4 = dim(Kâ‚„)
```
This is the Kâ‚„ spectator dimension â€” the same pattern as câ‚ƒ and câ‚‚.

## The Constraint Triangle

To simultaneously achieve câ‚ƒ=4, câ‚‚=6, câ‚_GUT=8:

| Requirement | Value | Constraint on Kâ‚ˆ sector |
|-------------|-------|------------------------|
| Î”câ‚ƒ = 0 | Must not change | SU(3)-singlet |
| Î”câ‚‚ = 0 | Must not change | SU(2)-singlet |
| Î”câ‚_SM = 14/3 | Must add exactly this | Tr(YÂ²) = 14/3 |

**The Kâ‚ˆ fermions must be SU(3)Ã—SU(2) singlets with nonzero hypercharge.**

## What Each Sub-Hypothesis Gives

### Full â„‚â¸ âŠ— â„‚â´ (all Kâ‚ˆ vertices Ã— Kâ‚„ spectator)
- Î”câ‚_SM = (7/6) Ã— 4 = 14/3 âœ“
- Î”câ‚ƒ = 1 Ã— 4 = 4 âœ— (doubles câ‚ƒ)
- Î”câ‚‚ = (1/4) Ã— 4 = 1 âœ— (shifts câ‚‚)
- **KILLED by SU(3) constraint**

### Hub-only â„‚Â² âŠ— â„‚â´ Ã— 3 generations
- Î”câ‚_SM = (1/2) Ã— 4 Ã— 3 = 6
- Î”câ‚ƒ = 0 âœ“ (hub is SU(3)-singlet)
- Î”câ‚‚ = ? (hub IS the SU(2) doublet â€” depends on identification)
- câ‚_GUT = (3/5)(26/3 + 6) = 8.8 â€” **overshoots by 10%**

### Ïâ‚‚ sector (SU(3)Ã—SU(2)-singlet, proved inaccessible from Higgs)
- Need: Tr(YÂ²|Ïâ‚‚) Ã— dim(Kâ‚„) Ã— 3 gen = 14/3
- â†’ YÂ² per Ïâ‚‚ state = 7/60 â‰ˆ 0.117
- â†’ Y = Â±0.342 (between 1/3 and 1/2)
- Î”câ‚ƒ = 0 âœ“, Î”câ‚‚ = 0 âœ“ (by construction)
- **Viable if Y assignment is correct**

## The Structural Picture

```
câ‚ƒ = Tr(TÂ²_su3|â„‚â¶) Ã— dim(â„‚â´) = 1 Ã— 4 = 4    [Kâ‚† gauge, Kâ‚„ spectator]
câ‚‚ = Tr(Ï„Â²_su2|â„‚â´) Ã— dim(â„‚â¶) = 1 Ã— 6 = 6    [Kâ‚„ gauge, Kâ‚† spectator]
câ‚ = Tr(YÂ²|â„‚â¶âŠ—â„‚â´) + Tr(YÂ²|singlet sector Ã— â„‚â´)
   = 26/3         + 14/3
   = 40/3                                         [Kâ‚†Ã—Kâ‚„ + Kâ‚ˆ singlet]
```

The U(1) is special because Y = Tâ‚ƒ_R + Q_BL/2 spans BOTH Kâ‚„ and Kâ‚†.
Neither serves as pure spectator. Kâ‚ˆ provides the additional states.

## Status

| Item | Status |
|------|--------|
| câ‚ƒ = 4 from â„‚â¶âŠ—â„‚â´ | **PROVED** |
| câ‚‚ = 6 from â„‚â¶âŠ—â„‚â´ | **PROVED** |
| câ‚_SM = 26/3 from â„‚â¶âŠ—â„‚â´ | **PROVED** |
| câ‚_GUT = 8 fits experiment (0.9%) | **OBSERVED** |
| Kâ‚ˆ vacuum is gauge singlet | **PROVED** |
| Double commutator nonzero on vacuum | **COMPUTED** |
| Coupling factor = 4 = dim(Kâ‚„) | **OBSERVED** |
| Î”câ‚ = 14/3 from Kâ‚ˆ singlet sector | **REQUIRED** (not derived) |
| Identification of which Kâ‚ˆ states contribute | **OPEN** |

## Next Move
The Ïâ‚‚ sector is the prime candidate. It is:
- Proved SU(3)Ã—SU(2)-singlet (from inaccessibility theorem)
- 2 states per generation Ã— 3 generations = 6 states
- These 6 states ARE the 6-fold vacuum degeneracy

Need to compute: what hypercharge does Ïâ‚‚ carry? 
If YÂ²(Ïâ‚‚) = 7/60 per state â†’ câ‚_GUT = 8 exactly.
