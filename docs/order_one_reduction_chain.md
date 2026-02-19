---
title: Order-One Condition and the Algebra Reduction Chain at Kâ‚ˆ
section: Gravitational Sector
status: active
---

# Order-One Condition and the Algebra Reduction Chain at Kâ‚ˆ

## Brian Porter â€” February 2026

---

## Summary

The order-one condition from noncommutative geometry, applied to the Kâ‚ˆ matching construction with the natural charge conjugation operator J (the Â±Î¼ eigenvalue swap), produces a **dimensional reduction chain that matches the NCG Standard Model at every level**:

| Level | Matching construction | NCG Standard Model | dim |
|:------|:---------------------|:-------------------|:----|
| 0 | R âŠ• Mâ‚‡(R) | â€” | 50 |
| 1 | Mâ‚„(R) âŠ• Mâ‚‚(R) âŠ• Mâ‚‚(R) | C âŠ• H âŠ• Mâ‚ƒ(C) | **24** |
| 2 | Mâ‚‚(R)Â² âŠ• RÂ² âŠ• RÂ² | Lie(A_F) | **12** |

The dimensions match exactly at levels 1 and 2. The algebra types differ (real vs complex/quaternionic), with the Zâ‚‡ Fourier complex structure providing the missing ingredient to bridge the gap.

---

## I. Setup: The Spectral Triple Data

**Hilbert space:** H = Râ¸ (vertex space of Kâ‚ˆ)

**Dirac operator:** D = S_vac (the vacuum matching matrix), symmetric with eigenvalues:
Î¼ = âˆ’0.693, âˆ’0.599, âˆ’0.372, âˆ’0.012, 0.000, +0.403, +0.567, +0.706

**Real structure:** J defined by the approximate Â±Î¼ pairing in the eigenbasis of D. The swap Ïƒ: i â†” 7âˆ’i pairs eigenvalues by magnitude:

| Pair | Î¼â» | Î¼âº | |Î¼âº/Î¼â»| |
|:-----|:----|:----|:---------|
| 0 | âˆ’0.693 | +0.706 | 1.018 |
| 1 | âˆ’0.599 | +0.567 | 0.947 |
| 2 | âˆ’0.372 | +0.403 | 1.085 |
| 3 | âˆ’0.012 | +0.000 | â‰ˆ0 (Majorana) |

The corresponding J satisfies JÂ² = I and {D, J} has norm only 0.068 (relative to ||D|| = 1.40), so **J nearly anticommutes with D** â€” the Â±Î¼ pairing is 95% exact. The residual Î´D = D âˆ’ Dâ‚€ (where Dâ‚€ has exact Â±Î¼ pairing) has ||Î´D||/||Dâ‚€|| = 2.4%, and is identified with the Majorana mass term.

## II. The Order-One Condition

The NCG order-one condition is [[D, a], JbJâ»Â¹] = 0 for all a, b âˆˆ A.

In the D-eigenbasis (where D = diag(Î¼â‚,...,Î¼â‚ˆ) and J acts as the permutation Ïƒ), this becomes an **index-set closure condition**:

If (i,j) and (k,l) are both in A's index set, and j = Ïƒ(k), then (i, Ïƒ(l)) must also be in A's index set (and symmetrically for the other term).

**Result:** The maximal subalgebras satisfying this closure are precisely the **Ïƒ-invariant block-diagonal subalgebras** â€” those whose blocks in the eigenbasis are unions of complete Ïƒ-pairs {i, Ïƒ(i)}.

The four Ïƒ-pairs give 15 possible block structures (the Bell number Bâ‚„ = 15). All 15 satisfy the order-one condition. Their dimensions range from 16 (four separate Mâ‚‚(R) blocks) to 64 (the full Mâ‚ˆ(R)).

## III. The Dimension-24 Subalgebras

Six of the 15 block structures have dimension 24, all of the form Mâ‚„(R) âŠ• Mâ‚‚(R) âŠ• Mâ‚‚(R). They differ in which two Ïƒ-pairs are merged into the Mâ‚„(R) block:

| # | Merged pairs | Mâ‚„ eigenvalues | Physical motivation |
|:--|:-------------|:---------------|:-------------------|
| 0 | (0,1) | Â±0.70, Â±0.58 | Largest eigenvalues |
| 1 | (0,2) | Â±0.70, Â±0.39 | Extreme + middle |
| 2 | (0,3) | Â±0.70, Â±0.01 | Vâ‚ƒ-dominant (Fourier) |
| 3 | (1,2) | Â±0.58, Â±0.39 | Middle pairs |
| 4 | (1,3) | Â±0.58, Â±0.01 | Vâ‚ + near-zero |
| 5 | (2,3) | Â±0.39, Â±0.01 | Vâ‚‚ + Majorana |

**The dimension 24 = dim_R(C âŠ• H âŠ• Mâ‚ƒ(C))** is an exact match with the NCG Standard Model finite algebra A_F.

Choice #2 (merging pairs 0 and 3) is distinguished by the Zâ‚‡ representation theory: both Ïƒ-pairs are Vâ‚ƒ-dominant in the Fourier basis, so merging them respects the Zâ‚‡ structure.

## IV. Level 2: The J-Commutant (dim 12)

Within each block, the commutant of Ïƒ (matrices commuting with charge conjugation) further reduces the algebra:

**Mâ‚„(R) block:** Ïƒ acts as the antidiagonal swap on Râ´ (no fixed points). The commutant consists of 4Ã—4 persymmetric matrices â€” those satisfying A[i,j] = A[3âˆ’i, 3âˆ’j]. This is an 8-dimensional subalgebra isomorphic to **Mâ‚‚(R) âŠ• Mâ‚‚(R)** (block-diagonal in the Ïƒ-eigenbasis, where Ïƒ = diag(+1,+1,âˆ’1,âˆ’1)).

**Mâ‚‚(R) blocks:** Ïƒ acts as the swap on RÂ². The commutant consists of persymmetric 2Ã—2 matrices [[a,b],[b,a]], which form a 2-dimensional algebra isomorphic to **R âŠ• R** (via the change of variables a Â± b).

**Total:** 8 + 2 + 2 = **12 = dim(su(3) âŠ• su(2) âŠ• u(1))**.

This matches the dimension of the SM gauge algebra.

## V. The Algebra Type Gap

While dimensions match at every level, the algebra types differ:

| Level | Matching | SM | Same dim? | Same type? |
|:------|:---------|:---|:----------|:-----------|
| 1 | Mâ‚„(R) âŠ• Mâ‚‚(R) âŠ• Mâ‚‚(R) | C âŠ• H âŠ• Mâ‚ƒ(C) | Yes (24) | No |
| 2 | Mâ‚‚(R)Â² âŠ• (RâŠ•R)Â² | Lie(A_F) | Yes (12) | No |

The fundamental difference: the matching construction is **purely real** (all algebras are M_n(R) or R), while A_F involves complex (C, Mâ‚ƒ(C)) and quaternionic (H) algebras.

## VI. The Zâ‚‡ Complex Structure as Bridge

The Zâ‚‡ action on Râ· provides a **natural complex structure** J_C on each 2-dimensional irrep V_k (k = 1, 2, 3):

J_k acts as rotation by 2Ï€k/7 on V_k, satisfying J_kÂ² = âˆ’I.

This complex structure is additional spectral triple data beyond D and J. Its effect on the algebras:

- End_R(V_k) = Mâ‚‚(R) [dim 4] â†’ End_C(V_k) = C [dim 2] under J_C
- Vâ‚€ and handle have no complex structure (they're 1-dimensional real)

When applied to the order-one subalgebra:

**Mâ‚„(R) block** â†’ Ïƒ-commutant Mâ‚‚(R) âŠ• Mâ‚‚(R):
- One Mâ‚‚(R) acts on the Ïƒ = +1 eigenspace (which may carry complex structure from Zâ‚‡) â†’ becomes C or stays Mâ‚‚(R) depending on Zâ‚‡ content
- The other Mâ‚‚(R) acts on the Ïƒ = âˆ’1 eigenspace â†’ similar reduction

**Mâ‚‚(R) blocks** â†’ Ïƒ-commutant R âŠ• R:
- With Zâ‚‡ complex structure: one R becomes part of C, the other stays R

The **asymmetric** application of complex structure (present on Vâ‚, Vâ‚‚, Vâ‚ƒ but not Vâ‚€ or handle) is precisely what distinguishes C âŠ• H âŠ• Mâ‚ƒ(C) from a symmetric product of matrix algebras.

## VII. The Majorana Mass Term

The eigenvalue pair (Î¼â‚ƒ, Î¼â‚„) = (âˆ’0.012, 0.000) has maximal charge conjugation asymmetry â€” one is exactly zero (the all-ones kernel) while the other is small but nonzero. This pair:

- Dominates the Vâ‚ƒ irrep (73% weight) and Vâ‚‚ (25%)
- Has |Î¼_near-zero|/|Î¼_max| = 0.017

In the SM Dirac operator, this corresponds to the Majorana mass term for the right-handed neutrino, which produces the see-saw mechanism for neutrino masses. The matching construction naturally generates this structure: the zero mode is the Zâ‚‡-invariant "spectator" (Vâ‚€ + handle), while the near-zero mode lives in Vâ‚ƒ â€” the same irrep that mediates the coupling chain between Vâ‚‚ (the hub) and Vâ‚.

## VIII. Summary of the Reduction Chain

```
R âŠ• Mâ‚‡(R)                              [dim 50]  Matching algebra
     â”‚
     â”‚ order-one (Ïƒ-invariant blocks)
     â–¼
Mâ‚„(R) âŠ• Mâ‚‚(R) âŠ• Mâ‚‚(R)               [dim 24]  = dim(A_F)
     â”‚
     â”‚ Ïƒ-commutant (charge conjugation)
     â–¼
(Mâ‚‚(R) âŠ• Mâ‚‚(R)) âŠ• (RâŠ•R) âŠ• (RâŠ•R)    [dim 12]  = dim(gauge algebra)
     â”‚
     â”‚ Zâ‚‡ complex structure (conjectured)
     â–¼
C âŠ• H âŠ• Mâ‚‚(C)  ???                    [dim ??]  â†’ A_F ?
```

**What works:**
- Dimensions 50 â†’ 24 â†’ 12 match the NCG framework exactly
- The order-one condition selects a 3-summand algebra from a 2-summand one
- Charge conjugation J arises naturally from Â±Î¼ pairing (95% exact)
- The Majorana mass term appears automatically as the asymmetric near-zero pair

**What's missing:**
- Algebra types are real (M_n(R)) instead of complex/quaternionic
- The Zâ‚‡ complex structure has not been formally incorporated into the order-one framework
- Mâ‚‚(C) at the bottom of the chain needs to become Mâ‚ƒ(C) â€” this likely requires larger Kâ‚‚â‚™
- The 6-fold degeneracy of dim-24 subalgebras needs a physical selection principle

## IX. Open Directions

1. **Incorporate J_C into the order-one framework.** The Zâ‚‡ complex structure J_C should enter as a grading or additional real structure in an extended spectral triple. The combined action of (D, J, J_C) may select a unique dim-24 subalgebra and fix its type.

2. **Higher Kâ‚‚â‚™.** At Kâ‚â‚‚ (p = 11), the Zâ‚â‚ representation has 5 non-trivial 2-dim irreps. The Mâ‚„(R) block could grow to Mâ‚†(R), whose Ïƒ-commutant Mâ‚ƒ(R) âŠ• Mâ‚ƒ(R) might, with complex structure, become Mâ‚ƒ(C) â€” exactly the color algebra. This is the most promising direction for recovering A_F exactly.

3. **The 6-fold degeneracy.** The six dim-24 subalgebras correspond to six choices of which Ïƒ-pairs to merge. In the SM, this might relate to the six types of quarks (up, down Ã— 3 colors) or to the broken flavor symmetry. The vacuum selection (Vâ‚‚ irrep) may break this degeneracy.

4. **Continuum limit.** As Kâ‚‚â‚™ â†’ âˆž, the matching algebra R âŠ• M_p(R) grows, the Ïƒ-pair structure becomes richer, and the order-one condition may converge to a unique algebra. Whether this limit is A_F = C âŠ• H âŠ• Mâ‚ƒ(C) is the central open question.
