---
title: Addendum: Discriminant Primes and the Heegner Observation
section: Higher Levels
status: active
---

# Addendum: Discriminant Primes and the Heegner Observation

## Appended to: K₁₀ Computation Results

---

## A1. The Congruence Constraint (Proved)

The square-free discriminant primes satisfy p ≡ 1 mod (2n−1) at K₈ and K₁₀:

| Level | 2n−1 | Sq-free disc p | p mod (2n−1) |
|:------|:------|:---------------|:-------------|
| K₆    | 5     | 5              | 0 (p = 2n−1) |
| K₈    | 7     | 43             | 1            |
| K₁₀   | 9     | 163            | 1            |

This has a structural explanation. The Gram matrix G commutes with the Z₂ₙ₋₁ permutation matrix P: [P, G] = 0 (verified computationally at all three levels). The eigenvalues of P are (2n−1)-th roots of unity. The discriminant of the minimal polynomial inherits ramification constraints from this cyclotomic structure. A prime appearing in the square-free discriminant must either divide 2n−1 (ramified in ℚ(ζ₂ₙ₋₁)) or satisfy p ≡ 1 mod (2n−1) (split in ℚ(ζ₂ₙ₋₁)). K₆ takes the first option (5 = 2n−1), K₈ and K₁₀ take the second. This is not a pattern — it's a theorem-level consequence of the Z₂ₙ₋₁ symmetry.

---

## A2. The Heegner Observation (Two Data Points)

The Heegner numbers are the nine values d ∈ {1, 2, 3, 7, 11, 19, 43, 67, 163} such that ℚ(√(−d)) has class number 1 (unique factorization in its ring of integers). Stark proved in 1967 that this list is complete.

The K₈ and K₁₀ square-free discriminants are both Heegner numbers:

- K₈: p = 43 is a Heegner number, and is the **unique** prime satisfying both p ≡ 1 mod 7 and h(ℚ(√(−p))) = 1
- K₁₀: p = 163 is a Heegner number, and is the **largest** prime satisfying both p ≡ 1 mod 9 and h(ℚ(√(−p))) = 1 (19 also qualifies but was not selected)

This is observed, not explained. Two data points do not make a pattern. The congruence constraint (§A1) is structural, but the class-number-1 selection within the congruence class has no known mechanism.

What makes the observation worth recording:

1. Among the infinitely many primes ≡ 1 mod 7, only one is Heegner: 43. The chain found it.
2. Among the infinitely many primes ≡ 1 mod 9, only two are Heegner: 19 and 163. The chain found the larger one.
3. If we look ahead: among primes ≡ 1 mod 11, only one is Heegner: 67. This yields a testable prediction.

---

## A3. Prediction for K₁₂

Combining the congruence constraint (structural) with the Heegner hypothesis (observational):

**If** the Heegner pattern holds at K₁₂, then the square-free discriminant prime is **67**, because 67 is the unique Heegner prime satisfying 67 ≡ 1 mod 11.

This is falsifiable at K₁₂. The congruence part (p ≡ 1 mod 11 or p | 11) will almost certainly hold. The Heegner part is the real test.

Additional K₁₂ predictions from the stabilization results (higher confidence):
- Minimal polynomial has degree 6 over ℚ
- Galois group is C₂ ≀ C₃ (order 24)
- Number field is Galois-disjoint from ℚ(λ₆), ℚ(λ₈), and ℚ(λ₁₀)
- Vacuum multiplicity is 6

---

## A4. The Exhaustion Boundary

The nine Heegner numbers are finite. The congruence-compatible ones are:

| Level | 2n−1 | Heegner primes ≡ 1 mod (2n−1) | Status |
|:------|:------|:-------------------------------|:-------|
| K₈    | 7     | {43}                           | Observed ✓ |
| K₁₀   | 9     | {19, 163}                      | Observed: 163 ✓ |
| K₁₂   | 11    | {67}                           | Predicted |
| K₁₄   | 13    | ∅                              | No candidates |
| K₁₆   | 15    | ∅                              | No candidates |

If the Heegner pattern is real, it **must end** at K₁₂. No Heegner prime is congruent to 1 mod 13. K₁₄ would be the first level where either the pattern breaks or the algebraic structure changes (degree, Galois group, or both).

This makes K₁₂ a doubly critical test: it can confirm the Heegner pattern (p = 67) or refute it, and K₁₄ would then test what happens when the Heegner supply runs out.

---

## A5. What's Solid vs. What's Speculative

**Proved or verified computationally:**
- Degree stabilization at 6 (K₈ and K₁₀ both degree 6)
- Galois group stabilization at C₂ ≀ C₃ (identical structure at K₈ and K₁₀)
- Pairwise Galois disjointness of all three number fields ℚ(λ₆), ℚ(λ₈), ℚ(λ₁₀)
- Scenario D: stabilized algebraic type with divergent arithmetic content
- Congruence constraint p ≡ 1 mod (2n−1) or p | (2n−1) (follows from [P, G] = 0)
- Square-free discriminant is a single prime at each level (observed, not explained)

**Observed, testable at K₁₂:**
- Discriminant primes are Heegner numbers (two data points: 43 and 163)
- Selection of largest Heegner candidate when multiple exist (one data point: 163 over 19)
- Prediction: K₁₂ square-free discriminant = 67

**Speculative (recorded for future reference, not relied upon):**
- Z₂ₙ₋₁ symmetry as analogue of complex multiplication
- Heegner exhaustion at K₁₄ as "arithmetic wall"
- Connection to j-invariants or period matrices of the embedding surfaces
- Any CM-theoretic explanation for the class-number-1 selection

---

## A6. Caveat on Embedding Dependence

This addendum inherits the caveat from §9 of the main document. The specific discriminant primes 43 and 163 were computed using direction assignments derived from Z₂ₙ₋₁ distance classes, not from verified surface embeddings. The K₈ computation (paperIII) used the unique Heawood embedding of K₇ on the torus. K₉ has no unique genus-3 embedding. A different choice of genus-3 embedding for K₉ could yield different direction classes and hence different discriminant primes.

The congruence constraint (§A1) should survive any Z₂ₙ₋₁-symmetric embedding choice. The Heegner observation (§A2) may or may not.
