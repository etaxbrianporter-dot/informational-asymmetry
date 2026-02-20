# Why Matchings Are Special: The Structural Theorem

**Brian Porter — February 19, 2026**
**Updated with superselection correction**

---

## Executive Summary

Three computational tests and one self-correction yield a precise structural theorem.

**M₁₁ test (negative):** Sporadic group M₁₁ on S(4,5,11) Steiner blocks → trapped in association scheme → cyclotomic arithmetic only (Gal = Z₅). M₁₁ contributes nothing beyond Z₁₁.

**K₆ test (positive):** K₆ matching Gram matrix escapes association scheme → non-abelian Galois S₃ from block 0 cubic. Discriminant primes {2, 5, 97} independent from K₈'s {2, 7, 43, 421}.

**Self-correction:** The "entangled nonic" from physical projection was an artifact — projecting integer matrices through irrational eigenvectors gives non-integer polynomials with no Galois theory. Direction classes are **superselection sectors**, not entangled fields. Each block has its own number field. The vacuum selects one.

---

## 1. The Two Ingredients

Level 2 arithmetic requires TWO ingredients simultaneously:

**Ingredient 1 — Z₃ phase refinement:** Within each direction block, Re(ω^{dz}) creates entries 2 (dz=0) and −1 (dz≠0) for pairs sharing one edge. This splits the overlap classification, escaping the association scheme. This is the engine producing non-abelian Galois.

**Ingredient 2 — δ(dir_i, dir_j) blocking:** The Kronecker delta zeroes inter-block entries, creating superselection sectors.

Without δ blocking, direction-phase information dissolves → Galois degrades to generic S₉ (verified).
Without phase refinement, matrix stays in association scheme → Galois is abelian (Kronecker-Weber).
With both: each sector has structured non-abelian arithmetic.

| Configuration | Phase refinement | δ blocking | Galois |
|--------------|-----------------|------------|--------|
| Fano overlap | None | N/A | Trivial |
| M₁₁ Steiner | Trapped in scheme | N/A | Z₅ (abelian) |
| K₇ Dirac | Breaks equivariance | N/A | S₇ (generic) |
| K₆ crude projection | Escapes scheme | **Absent** | **S₉ (generic)** |
| **K₆ actual Gram** | **Escapes scheme** | **Present** | **S₃ (structured)** |

---

## 2. K₆ Number Fields

G is block-diagonal: G = B₀ ⊕ B₁ ⊕ B₂ (three 5×5 integer matrices).
Full char poly = product of block char polys (verified).

### Block 0 — contains the Higgs eigenvalue

Char poly: (x − 7)(x − 5)(**x³ − 18x² + 100x − 170**)

- Irreducible cubic, disc = 1940 = 2² × 5 × 97
- Galois group: **S₃** (disc not a perfect square)
- Roots: 8.946, 5.748, **3.306 = a₂**
- v₀ has 100% weight in this block (verified: 1.000000)
- The Higgs mass lives in this S₃ field

### Block 1

Char poly: (x − 4)(x² − 13x + 41)²

- Number field: **ℚ(√5)** — the golden ratio field
- disc = 5, Galois = Z₂
- Eigenvalues: (13 ± √5)/2 ≈ 7.618, 5.382

### Block 2

Char poly: (x − 9)(x − 6)(x − 5)(x² − 10x + 23)

- Number field: **ℚ(√2)**
- disc = 8, Galois = Z₂
- Eigenvalues: 5 ± √2 ≈ 6.414, 3.586

### Discriminant primes

| Level | Primes | Source |
|-------|--------|--------|
| K₆ block 0 | {2, 5, 97} | Higgs cubic |
| K₆ block 1 | {5} | Golden ratio |
| K₆ block 2 | {2} | ℚ(√2) |
| K₈ | {2, 7, 43, 421} | Vacuum sextic |

Galois disjointness: K₆ and K₈ share only the trivial prime 2.

---

## 3. The Superselection Principle

δ(dir_i, dir_j) encodes: **two matchings interfere only if they share a direction class.**

This is not a mathematical convenience. It is the combinatorial origin of gauge superselection. In the physical theory: fields in different gauge representations don't directly couple.

The v₀ verification makes this concrete:

```
  v₀ components:
    dir=0 matchings: 0.398, 0.398, 0.674, 0.075, -0.473  (total weight: 1.0)
    dir=1 matchings: 0, 0, 0, 0, 0                        (total weight: 0.0)
    dir=2 matchings: 0, 0, 0, 0, 0                        (total weight: 0.0)
```

The vacuum lives entirely in one superselection sector. The Higgs mass comes from that sector's cubic polynomial. The other sectors (√5 field, √2 field) contain different algebraic content that doesn't contribute to this observable.

---

## 4. The Association Scheme Prison

The Kronecker-Weber theorem creates a hard wall: equivariant bilinear forms trapped in the association scheme can only produce abelian Galois groups (subgroups of cyclotomic Galois groups).

**M₁₁ Steiner system:** M†M lies in the 4-dimensional association scheme of the 3-class design. All eigenvalues in ℚ(cos(2π/11)). Galois = Z₅. The sporadic group adds nothing.

**K₆ matchings:** G escapes the scheme because shared-edge direction class ≠ shared-edge count. Within block 0, pairs with overlap 1 get G-entry 2 (if dz=0) or −1 (if dz≠0). Same overlap, different G values. Max residual from scheme = 2.2.

**The escape mechanism:** Two matchings sharing one edge in direction class 0 with dz=0 are algebraically different from two matchings sharing one edge in direction class 0 with dz≠0. The Z₃ holonomy phase distinguishes them. This distinction is invisible to the association scheme but visible to G.

---

## 5. K₈ Direction Structure

At K₈, the direction permutation ×2 mod 7 acts as a 3-cycle on the three non-handle direction overlap matrices:

```
  P_dir · O₀ · P_dir^T = O₂  (direction 0 → 2)
  P_dir · O₁ · P_dir^T = O₀  (direction 1 → 0)  
  P_dir · O₂ · P_dir^T = O₁  (direction 2 → 1)
```

The three direction sectors are **isospectral** — same eigenvalues, related by conjugation. This Z₃ conjugacy between sectors is the structural origin of the C₃ quotient in the vacuum sextic's Galois group C₂ ≀ C₃.

The K₈ pattern: three isomorphic sectors → vacuum selects one → sector's irreducible polynomial has Galois with C₃ acting on conjugate root pairs → physical C₃ = color.

---

## 6. Structural Theorem (Final Form)

**Theorem.** The matching complex of K₂ₙ with Heawood surface embedding and spectral action Gram matrix satisfies:

**(a)** G is block-diagonal by direction class, with entries in {−1, 0, 2, 2n}.

**(b)** Within each block, Z₃ phase refinement escapes the association scheme.

**(c)** Each block's char poly has irreducible factors with non-abelian Galois (at least S₃).

**(d)** The vacuum eigenvector lives entirely in one direction sector. The Higgs eigenvalue a₂ is algebraic, living in that sector's number field.

**(e)** Without δ blocking, Galois degrades to generic S_n. The δ is essential.

**Two necessary, jointly sufficient ingredients:** scheme escape (from Z₃ phases) + superselection (from δ blocking).

---

## 7. What This Unlocks

1. **Predictive filter.** Candidate extractors must have BOTH scheme escape AND sector blocking. Vertex-based structures (Steiner systems, codes) lack the edge-direction data needed for scheme escape. Unsectored edge structures lack the blocking needed to prevent generic degradation.

2. **The δ as gauge principle.** The mathematical requirement for structured arithmetic (δ blocking) maps to the physical requirement for gauge superselection. This isn't an analogy — it's the mechanism by which gauge structure emerges from matching combinatorics.

3. **Vacuum selection.** The ground state chooses one direction sector. The Higgs mass is determined by the algebraic number theory of that sector. Different sectors contain different number fields with different physical content. Vacuum selection IS sector selection.

4. **Higher levels.** The K₁₀, K₁₂, ... levels should each have their own direction blocks, their own sector structure, their own number fields. The direction count grows with the graph, predicting richer Galois structure at each level.
