---
title: Complete Discriminant Atlas: K₆ through K₁₂
section: Arithmetic Structure
status: active
---

# Complete Discriminant Atlas: K₆ through K₁₂

## Full Discriminant Factorizations

| Level | m = 2n−1 | deg | disc(f) | Sq-free Δ |
|:------|:---------|:----|:--------|:----------|
| K₆ | 5 | 2 | 2² × 5 | 5 |
| K₈ | 7 | 6 | 2³⁶ × 7⁴ × 43 × 421² | 43 |
| K₁₀ | 9 | 6 | 2⁷⁸ × 3¹⁰ × 17² × 163 | 163 |
| K₁₂ | 11 | 10 | 2¹⁹⁰ × 11⁸ × 23 × 353 × 439² × 64217² × 478838249² | 8119 = 23 × 353 |

## New Results (this session)

1. **439² ∥ disc(K₁₂)** — confirmed as genuine spectral prime, double root at x ≡ 24 mod 439
2. **64217² ∥ disc(K₁₂)** — newly discovered, double root at x ≡ 21290 mod 64217
3. **478838249² ∥ disc(K₁₂)** — newly discovered, both are prime
4. **disc(K₈) = 2³⁶ × 7⁴ × 43 × 421²** — complete factorization via sympy
5. Congruence law validated (one-directional); biconditional form falsified by 421 at K₈

## The Congruence Law — Correct Statement

**Theorem (verified K₆–K₁₂).** Let f₂ₙ(x) be the vacuum characteristic polynomial at level K₂ₙ, and let m = 2n−1. If p is a prime with v_p(disc(f₂ₙ)) odd and p ∤ 2m, then p ≡ 1 mod m.

**Equivalently:** The square-free part of disc(f₂ₙ) is supported on {2, m} ∪ {primes p ≡ 1 mod m}.

**Evidence:**
- K₆: sq-free = 5 (level prime)
- K₈: sq-free = 43, and 43 ≡ 1 mod 7 ✓
- K₁₀: sq-free = 163, and 163 ≡ 1 mod 9 ✓
- K₁₂: sq-free = 23 × 353, and 23 ≡ 353 ≡ 1 mod 11 ✓

**The biconditional is false:** p = 421 at K₈ satisfies 421 ≡ 1 mod 7 (would predict odd exponent) but v₄₂₁(disc) = 2 (even). So p ≡ 1 mod m does NOT force odd multiplicity.

## Splitting Type Analysis (K₈)

The distinction between odd and even multiplicity relates to how completely the polynomial splits mod p:

**p = 43 (v = 1, odd):** The sextic splits into 6 linear factors mod 43:
(x−24)²(x−7)(x−35)(x−41)(x−42). Complete splitting; the double root contributes v_p = 1.

**p = 421 (v = 2, even):** The sextic splits as (x−407)²(x−55)(x−73) × (irred. quadratic mod 421). Partial splitting: 3 roots visible, 1 double, but an irreducible quadratic persists. The residual non-split factor pushes v_p from 1 to 2.

**Mechanism:** The squared discriminant factor arises when a spectral prime simultaneously ramifies (double root) AND fails to split completely (residual irreducible factor). At p = 43, there's ramification but complete splitting, giving v = 1. At p = 421, there's ramification with incomplete splitting, giving v = 2.

## Discriminant Prime Classification

For each non-structural prime p dividing disc(f₂ₙ):

| Category | Congruence | Exponent | Splitting | Examples |
|:---------|:-----------|:---------|:----------|:---------|
| Structural | p = 2 | high even | — | 2 at all levels |
| Level prime | p = m | high | — | 5, 7, 3, 11 |
| Sq-free | p ≡ 1 mod m | odd | complete + 1 double | 43(K₈), 163(K₁₀), 23,353(K₁₂) |
| Squared-split | p ≡ 1 mod m | even | partial + 1 double | 421(K₈) |
| Squared-inert | p ≡ −1 mod m | even | partial + 1 double | 17(K₁₀), 439,64217,478838249(K₁₂) |

The "squared-split" category (421 at K₈) is the one that breaks the biconditional. It shows that even ≡1 primes can appear squared, when the relative discriminant of the quadratic extension over the cubic subfield ramifies there.

## Connection to Field Tower

For C₂ ≀ C_k groups, the number field has a tower structure:

ℚ → ℚ(θ) → ℚ(λ)  (degree k → degree 2k)

The discriminant of the sextic decomposes as:

disc(f) = disc(cubic)² × N_{ℚ(θ)/ℚ}(relative disc)

- disc(cubic) contributes the structural primes (2, m)
- The relative discriminant contributes the spectral primes
- Primes that split in the cubic AND ramify in the quadratic appear with even exponent
- Primes that DON'T fully split in the cubic but still create a double root appear with odd exponent

This explains the one-directional nature: v_p odd requires the prime to not fully split in the base field, which (for these cyclotomic-adjacent extensions) forces p ≡ 1 mod m. But p ≡ 1 allows both complete and partial splitting, so either parity is possible.

## Open Questions

1. Does K₁₄ (m=13) follow the same pattern? Sq-free = 417041 = 79 × 5279, both ≡ 1 mod 13. ✓ But what are the squared primes?

2. Is there a formula predicting WHICH p ≡ 1 primes appear squared vs simple? (Likely related to splitting behavior in the cubic subfield.)

3. The 2-adic valuation grows rapidly: 2, 36, 78, 190. Is there a formula for v₂(disc(f₂ₙ))?

4. The level-prime exponents: 5¹, 7⁴, 3¹⁰, 11⁸. Pattern?
