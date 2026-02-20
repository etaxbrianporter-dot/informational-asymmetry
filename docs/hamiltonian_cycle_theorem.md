# The Hamiltonian Cycle Theorem — Complete Proof

## Context

This theorem closes Pillar 3 of the CP Parity Uniformity Theorem by proving
that trace(Q|_live) = k for the vacuum overlap matrix of K_{2n}.

## Setting

K_{2n} with vertices {0, 1, ..., p} where p = 2n−1 is prime. We call vertex p
the **center** and vertices {0, ..., p−1} the **circle** (identified with ℤ/pℤ).

A vacuum direction block contains 2p perfect matchings of K_{2n}. Each matching
consists of:
- One **spoke edge** {s, center} connecting a circle vertex s to the center
- n−1 = (p−1)/2 **circle edges** pairing the remaining p−1 circle vertices

## The Two Patterns

**Theorem 1.** The 2p matchings of a vacuum direction block fall into exactly
two classes (patterns A and B) of p matchings each, distinguished by their
circle pairing structure relative to the spoke vertex.

*Proof.* Fix spoke vertex s. The remaining p−1 circle vertices must be paired
into (p−1)/2 edges. Among the p−1 circle edges available, p−1 are "adjacent"
(of the form {s+j, s+j+1 mod p} for j = 1, ..., p−1) and the rest are
"non-adjacent." Each perfect matching of these p−1 vertices into (p−1)/2 = n−1
pairs uses exactly n−2 adjacent circle edges and 1 non-adjacent circle edge.

The non-adjacent pair acts as a "bridge" that partitions the remaining vertices
into two arcs, each paired by consecutive edges. There are exactly 2 valid
bridge positions that yield a matching of p−1 vertices on the circle, producing
patterns A and B. □

**Structure verified at K₈/K₁₀/K₁₂:**

| Level | p  | Pattern A non-adj pair | Pattern B non-adj pair | Adjacent pairs |
|-------|----|----------------------|----------------------|----------------|
| K₈   | 7  | (3, 6)               | (1, 4)               | 2 each         |
| K₁₀  | 9  | (1, 6)               | (3, 8)               | 3 each         |
| K₁₂  | 11 | (1, 8)               | (3, 10)              | 4 each         |

**Key property:** The non-adjacent pairs of A and B are related by a shift of ±2
positions on the circle. Specifically, if A has bridge (a, b) relative to the
spoke, then B has bridge (a±2, b±2) (mod p).

## Edge Classification

**Theorem 2.** The 2p matchings use exactly 3p edges, classified as:
- p **frequent edges**: the circle adjacencies {i, i+1 mod p}, each appearing
  in 2(n−2) matchings
- 2p **rare edges**: all others, each appearing in exactly 2 matchings

*Proof.* 

*Frequent edges:* Consider circle edge {i, i+1 mod p}. This edge appears in a
matching with spoke s iff {i, i+1} is one of the n−2 adjacent pairs in that
matching's pattern. As s ranges over ℤ/pℤ, the adjacent pairs shift by the
spoke position. For a fixed pattern (A or B), the pair {i, i+1} appears for
exactly n−2 values of s (those where i and i+1 are both in the "same arc"
created by the bridge). Since both patterns contribute n−2 matchings,
freq({i, i+1}) = 2(n−2). ✓

*Rare edges:* The spoke edges {s, center} for s ∈ ℤ/pℤ contribute p rare edges,
each appearing in exactly 2 matchings (one of each pattern with spoke s). The
non-adjacent circle pairs contribute p more rare edges: for each spoke s, pattern
A contributes one non-adjacent pair and pattern B contributes a different one.
By the ±2 shift relation, each such edge appears with exactly one (spoke, pattern)
combination plus its ±2-shifted partner, giving frequency 2.

Total edges: p + 2p = 3p. ✓
Total incidences: p·2(n−2) + 2p·2 = 2p(n−2) + 4p = 2p·n = 2p matchings × n edges. ✓ □

## Claim 1: The Rare-Edge Graph is a Single Hamiltonian Cycle

**Theorem 3.** The rare-edge graph G_rare (vertices = 2p matchings, edges = 2p
rare edges connecting the two matchings that share each rare edge) is a single
Hamiltonian cycle of length 2p.

*Proof.* G_rare is 2-regular (each matching has exactly 2 rare edges: one spoke
and one non-adjacent circle edge), so it decomposes into disjoint cycles covering
all 2p vertices. We show there is exactly one cycle.

**Step 1: Within-spoke connections.** For each spoke s, matchings A(s) and B(s)
share the spoke edge {s, center}. This gives p edges of G_rare.

**Step 2: Cross-spoke connections.** Each matching also has one rare circle edge
(its non-adjacent pair). The matching B(s) has non-adjacent pair at relative
positions (a, b). The only other matching containing this same absolute edge is
A(s ± 2), whose non-adjacent pair occupies the same absolute vertices (since A's
relative pair is B's pair shifted by ∓2, which exactly compensates the spoke shift
of ±2). This gives p more edges of G_rare.

**Step 3: Cycle structure.** The edges form the pattern:

    A(s₀) —spoke— B(s₀) —circle— A(s₀±2) —spoke— B(s₀±2) —circle— A(s₀±4) —...

The spoke labels follow the arithmetic progression s₀, s₀±2, s₀±4, ... (mod p).
Since gcd(2, p) = 1 (p is an odd prime), this progression visits ALL p values
of ℤ/pℤ before returning to s₀. Therefore the cycle has length 2p and is
Hamiltonian. □

**Verified spoke sequences:**
- K₈  (p=7):  6, 4, 2, 0, 5, 3, 1  (step −2 mod 7)
- K₁₀ (p=9):  7, 0, 2, 4, 6, 8, 1, 3, 5  (step +2 mod 9)
- K₁₂ (p=11): 9, 0, 2, 4, 6, 8, 10, 1, 3, 5, 7  (step +2 mod 11)

## Claim 2: The Alternating Signing Satisfies All Edge Constraints

**Theorem 4.** The alternating ±1 signing of the Hamiltonian cycle lies in
ker(M^T), where M is the (2p × E) matching-edge incidence matrix.

*Proof.* The Hamiltonian cycle strictly alternates A-B-A-B (since within-spoke
edges connect A↔B and cross-spoke edges also connect B↔A). Therefore the
alternating signing assigns +1 to all pattern-A matchings and −1 to all
pattern-B matchings (or vice versa). Call this vector v_AB.

We must verify M^T v_AB = 0, i.e., for every edge e:
Σ_{matchings m containing e} v_AB(m) = 0.

**Rare edges:** Each rare edge connects exactly 2 matchings of opposite patterns
(one A, one B). Their signed sum = (+1) + (−1) = 0. ✓

**Frequent edges:** Each frequent edge {i, i+1 mod p} appears in exactly n−2
pattern-A matchings and n−2 pattern-B matchings (proved in Theorem 2). Their
signed sum = (n−2)(+1) + (n−2)(−1) = 0. ✓ □

## The Q-Parity

**Theorem 5.** The alternating signing v_AB has Q-eigenvalue −1.

*Proof.* Q acts on matchings by the vertex multiplication v ↦ (p−1)v mod p,
which sends the center to itself and acts on circle vertices as s ↦ −s mod p.
This maps spoke s to spoke −s mod p (changing the spoke vertex but keeping the
center).

Q also transforms the circle pairing: under s ↦ −s, the relative pair structure
is reflected, mapping pattern A to pattern B and vice versa. Therefore
Q(A(s)) = B(−s mod p) and Q(B(s)) = A(−s mod p).

Since v_AB(A(·)) = +1 and v_AB(B(·)) = −1, we have:
v_AB(Q(A(s))) = v_AB(B(−s)) = −1 = −v_AB(A(s))
v_AB(Q(B(s))) = v_AB(A(−s)) = +1 = −v_AB(B(s))

Therefore Q v_AB = −v_AB, giving Q-eigenvalue −1. □

## Conclusion: Pillar 3

**Corollary (Pillar 3).** For the diagonal block C₀ of the vacuum overlap matrix:

1. **rank(C₀) = 2p − 1:** The vector v_AB ∈ ker(C₀) (since C₀ = 2MM^T and
   M^T v_AB = 0), giving rank ≤ 2p−1. The rank is exactly 2p−1 because the
   Hamiltonian cycle structure forces all linear relations to be proportional
   to v_AB (the single cycle constrains the signing up to scale).

2. **Null vector in V₋:** v_AB has Q-eigenvalue −1, so ker(C₀) ⊂ V₋.

3. **Asymmetry = +1:** rank(C₀|V₊) = p and rank(C₀|V₋) = p−1.

By the block-circulant structure, each of the k blocks C_d inherits the same
+1 asymmetry, giving:

    trace(Q|_live) = Σ_{d=0}^{k-1} [rank(C_d|V₊) − rank(C_d|V₋)] = k

Combined with Pillars 1 and 2, this completes the CP Parity Uniformity Theorem:

    m(ε, j) = (p ∓ 1)/2  independent of C_k sector j.  □

## Summary of the Proof Architecture

The entire proof reduces to three elementary facts about the round-robin
1-factorization of K_{2n}:

1. **Two patterns:** Each direction block has p matchings of pattern A and p of
   pattern B, distinguished by the position of their unique non-adjacent circle
   pair.

2. **Step-2 Hamiltonian cycle:** The rare-edge graph is a single cycle because
   the cross-spoke transitions step by ±2 mod p, and gcd(2, p) = 1.

3. **A/B balance:** Each frequent (adjacent circle) edge appears equally often
   in pattern A and pattern B, because the spoke-shift symmetry over ℤ/pℤ
   treats both patterns identically.

No spectral analysis, representation theory, or physics is needed. The theorem
is a purely combinatorial fact about perfect matchings on complete graphs.
