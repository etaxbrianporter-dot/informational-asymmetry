---
title: Binary Coloring of K₆: What the Vacuum Sees
section: K₈ Fermion Sector
status: active
---

# Binary Coloring of K₆: What the Vacuum Sees

## The Discovery

Color each of K₆'s 15 edges RED (appears in ≥1 active matching) or BLUE (silent). The result:

```
11 RED edges, 4 BLUE edges

BLUE edges = {(0,3), (0,5), (1,3), (1,5)} = K₂,₂ between {0,1} and {3,5}
```

The 4 silent edges form a perfect bipartite graph connecting the **doublet** {0,1} to a specific **generation pair** {3,5}.

## The Two Hub Edges

The vacuum has TWO edges each appearing in 3 of the 5 active matchings:

| Hub edge | In matchings | Role |
|----------|-------------|------|
| (0,1) | M₀, M₁, M₂ | Doublet edge |
| (3,5) | M₁, M₄, M₉ | Generation edge |

Each hub is in exactly 3 = matchings(K₄) active matchings. The counting:
- 5 matchings × 3 edges = 15 edge-slots
- 2 hubs × (3−1) excess = 4 wasted slots  
- 15 − 4 = 11 distinct edges covered, 4 silent

## The Three-Pair Hierarchy

The dominant matching M₁ = (01)(24)(35) partitions K₆ into three pairs with distinct roles:

| Pair | Red degree | Blue degree | Role | Y†Y eigenvalue |
|------|-----------|------------|------|---------------|
| {0,1} | 3 each | 2 each | Doublet | — |
| {2,4} | 5 each | 0 each | Heavy generation | ~1.172 |
| {3,5} | 3 each | 2 each | Light generation | ~0.145 |

Vertices {2,4} are **universal connectors** — every edge touching them is active. They participate in all interactions. Vertices {0,1} and {3,5} have the same degree pattern (3 red, 2 blue) but different structural positions: {0,1} is the doublet hub, {3,5} is the generation hub.

## The Ramsey Structure

The blue subgraph K₂,₂ is a 4-cycle — triangle-free. By R(3,3) = 6, the red subgraph MUST contain triangles. It contains 8: all 4 triangles within the generation K₄, plus 4 mixed triangles involving one doublet vertex.

The generation space K₄ on {2,3,4,5} is **entirely red** — a complete active subgraph.

## The Identity That Kills Distinguishability

Two hypotheses for the "4" in λ/g² = a₄/(4a₂):

**(A)** dim(so(5)/so(4)) = 4 — the Higgs coset dimension  
**(B)** |V(K₆)| − 2 = 4 — vertices minus the doublet

These are **algebraically identical for all K₂ₙ**:

```
dim(so(2n−1)/so(2n−2)) = (2n−2)(2n−1)/2 − (2n−2)(2n−3)/2
                        = (2n−2) × 1 = 2n − 2 = |V(K₂ₙ)| − 2
```

For K₄: both give 2. For K₆: both give 4. For K₈: both give 6.

The number cannot distinguish the mechanisms. The question is which **derivation path** is correct — top-down from the Lie algebra, or bottom-up from the vacuum's edge selection.

## What This Means

The binary coloring reveals that the vacuum doesn't just pick a direction or a matching — it picks a **graph-theoretic structure**: two hub edges creating a K₂,₂ of silence in a sea of activity. The 4 silent edges are the suppressed Yukawa channels between the Higgs doublet and the lighter generation.

The factor of 4 appears as:
- 4 silent edges (= K₂,₂ size)
- 4 generation vertices (= |V| − 2)
- 4 Higgs components (= dim(coset))
- All algebraically the same number, all arising from removing a single edge from K₆
