---
title: K₅: The Shadow Graph Between Spacetime and Matter
section: K₄ Spacetime
status: active
---

# K₅: The Shadow Graph Between Spacetime and Matter

## Computation Results — February 17, 2026

---

## The Discovery

K₄ (spacetime) and K₆ (matter) are the two physical graphs in the K₄D framework. Between them sits K₅ — a graph that **cannot** appear as a physical factor (no Pfaffian, no perfect matchings, odd vertex count) but whose shadow is everywhere in the relations between the two physical graphs.

---

## 1. Parity Twinning (PROVED)

**Theorem.** For N even, |b₀/a₀|(K_N) = |b₀/a₀|(K_{N+1}) = 1/(N+1).

*Proof.* K_N even: b₀ = −N, a₀ = N(N+1), ratio = 1/(N+1).
K_{N+1} odd: b₀ = N+2, a₀ = (N+1)(N+2), ratio = 1/(N+1). ∎

**Consequence:** K₄ and K₅ are **parity twins** — they share |b₀/a₀| = 1/5 exactly.

| Graph | a₀ | b₀ | |b₀/a₀| | Role |
|-------|-----|-----|---------|------|
| K₄ | 20 | −4 | 1/5 | Spacetime (Pfaffian, matchings, physical) |
| K₅ | 30 | +6 | 1/5 | Shadow (no Pfaffian, no matchings, structural) |

This means the birefringence coupling β/α_T = b₀/a₀ = −1/5 is not an isolated property of K₄. It's a **pair invariant** shared with K₅. The twin inherits the asymmetry without the geometry.

Similarly, K₆ and K₇ are twins with |b₀/a₀| = 1/7.

---

## 2. ℤ₃ Flux Twinning (PROVED)

The number of triangles C(N,3) mod 3 determines the net ℤ₃ flux:

| Graph | Triangles | mod 3 | Net ℤ₃ flux |
|-------|-----------|-------|-------------|
| K₃ | 1 | 1 | 2π/3 |
| K₄ | 4 | 1 | 2π/3 |
| K₅ | 10 | 1 | 2π/3 |
| K₆ | 20 | 2 | 4π/3 ≡ −2π/3 |
| K₇ | 35 | 2 | −2π/3 |

K₃, K₄, and K₅ all carry the **same** net ℤ₃ flux. K₆ carries the **opposite**. The parity twins are also flux twins.

---

## 3. The Higgs Coset (STRUCTURAL)

The Lie algebra chain:

```
so(4) ⊂ so(5) ⊂ so(6)
  K₄      K₅      K₆
space   bridge   matter
```

where dim(so(N)) = |E(K_N)| always.

**The key:** K₅ = K₄ + star(vertex 5). Adding the 5th vertex creates 4 new edges.

|E(K₅)| − |E(K₄)| = 10 − 6 = **4** = dim_ℝ(Higgs doublet)

These 4 edges correspond to generators J₁₅, J₂₅, J₃₅, J₄₅ of so(5)/so(4).

The coset so(5)/so(4) = S⁴ is exactly the **Higgs vacuum manifold** in SO(5) gauge-Higgs unification (Manton). The Higgs doublet has 4 real components = 1 complex doublet (H⁺, H⁰).

**The Higgs IS the edge-difference between K₅ and K₄.**

In NCG language: gauge fields come from inner automorphisms (the so(4) part), while the Higgs comes from inner fluctuations (the so(5)/so(4) coset). K₅ provides exactly this decomposition.

---

## 4. The Matching Bridge (PROVED)

matchings(K₆) / matchings(K₄) = 5!! / 3!! = 15/3 = **5** = |V(K₅)|

General: (2n+1)!! / (2n−1)!! = 2n+1 = |V(K_{2n+1})|

The odd graph K_{2n+1} is the **multiplicative bridge** between the matching counts of K_{2n} and K_{2n+2}. K₅ is the factor by which K₆'s Pfaffian has more terms than K₄'s.

---

## 5. The Petersen Graph (COMBINATORIAL)

The Petersen graph = Kneser graph K(5,2):
- **Vertices** = 2-element subsets of V(K₅) = edges of K₅ → **10** = |E(K₅)|
- **Edges** = pairs of disjoint K₅ edges → **15** = |E(K₆)|

Properties:
- 3-regular (degree 3 = N_c = QCD colors)
- Chromatic number 3 (needs exactly N_c colors)
- Girth 5 (shortest cycle = pentagon = |V(K₅)|)

The Petersen graph encodes the **matching structure** of K₅. Its vertex count equals K₅'s edge count, and its edge count equals K₆'s edge count. It literally bridges the combinatorics of the two physical graphs.

**Stars of K₅:** Each vertex v ∈ K₅ has 4 incident edges, forming an independent set of size 4 in the Petersen graph. There are 5 such stars, and 5 × 4 = 20 = a₀(K₄).

---

## 6. Genus Boundary (TOPOLOGICAL)

| Graph | Genus | Surface | Physical role |
|-------|-------|---------|--------------|
| K₄ | 0 | Sphere S² | Spacetime — last planar |
| K₅ | 1 | Torus T² | Shadow — first toroidal |
| K₆ | 1 | Torus T² | Matter — still toroidal |
| K₇ | 1 | Torus T² | K₆'s twin — still toroidal |
| K₈ | 2 | Double torus | Beyond physical chain |

K₅ marks the **genus transition** from S² to T². The K₄D boundary IS a torus. K₅ is the first complete graph that requires the same topology as the boundary.

---

## 7. Line Graph L(K₅) and the Boundary Lattice

L(K₅) has:
- 10 vertices (= |E(K₅)|)
- 30 edges (= a₀(K₅))
- **6-regular** = same coordination number as the triangular lattice

The boundary theory (K₄ on triangular lattice) lives on a 6-regular graph. L(K₅) is the minimal 6-regular graph in the line-graph family. The boundary lattice is, in a precise sense, the thermodynamic limit of L(K₅)-like graphs.

---

## 8. The Birefringence Formula Decomposition

$$\alpha_0 = \frac{(2-\sqrt{3})^2}{10} \cdot F(m_A)$$

Every factor is a graph invariant:

| Factor | Value | Origin |
|--------|-------|--------|
| (2−√3)² | 0.0718 | K₄ D² eigenvalue ratio = 1/R |
| 1/5 | 0.2000 | |b₀/a₀| = K₄−K₅ twin pair invariant |
| 1/2 | 0.5000 | Initial displacement convention (= dim(doublet)/dim(coset)?) |
| 10 | denom | 2 × |V(K₅)| = |E(K₅)| = dim(so(5)) |
| F(m_A) | ~0.73 | Cosmological evolution (not a graph invariant) |

---

## 9. What K₅ Cannot Do

- **No Pfaffian:** odd dimension → no Pfaffian → no signature selection → no spacetime
- **No perfect matchings:** odd vertex count → 0 perfect matchings → no gauge algebra
- **No physical factor:** cannot appear in the spectral triple K₄ × K₆ × K₄

---

## 10. What K₅ Does

- **Bridges** K₄ and K₆ through every algebraic and topological invariant
- **Provides** the Higgs coset: so(5)/so(4) = S⁴ = 4 generators = Higgs doublet
- **Shares** K₄'s parity ratio (1/5) and ℤ₃ flux (2π/3) exactly
- **Encodes** the matching ratio 15/3 = 5 between the two physical graphs
- **Marks** the genus transition from S² to T² (sphere → boundary torus)
- **Generates** the Petersen graph as its matching structure, with |V| = |E(K₅)|, |E| = |E(K₆)|

**K₅ is the invisible graph that connects spacetime to matter. It cannot exist as physics, but it IS the relation between the two things that can.**

---

## Status and Implications

**PROVED:** Parity twinning, ℤ₃ flux twinning, matching bridge, Petersen connection, genus boundary, line graph coordination

**STRUCTURAL (not proved but mathematically precise):** Higgs coset so(5)/so(4) = 4 generators = 4 real Higgs components = |E(K₅)| − |E(K₄)|

**OPEN:** Whether the so(5) coset structure can be derived from the spectral action formalism (rather than identified post hoc), and whether K₅ plays a role in the normalization constant c that converts K₆ spectral ratios to the physical Higgs mass.

**KILLED:** The 7/8 coincidence (b₀²/a₀(K₇) = 8/7, not 7/8; the K₆ spectral ratio 7/8 has a different origin from the physical D†D eigenvalues).
