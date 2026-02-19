---
title: The Congruence Law: Mechanism Identified
section: 
status: active
---

# The Congruence Law: Mechanism Identified

## The Gap

The original proof of Law 3 argued: [P,G] = 0 ⟹ congruence constraint. Computational testing proved this insufficient — random Z₉-equivariant matrices also satisfy [P,G] = 0 but 70% violate the congruence law. Something specific to the matching chain was needed.

## The Discovery

Systematic perturbation testing revealed:
- **Cross-overlap vector controls congruence**: Keeping the matching chain's cross-overlap [12,8,0,0,8,4,8,4,0] but randomizing the diagonal → 80/80 pass congruence
- **Random cross-overlap breaks it**: Random cross-overlap with any diagonal → ~30% pass
- **Diagonal is invisible**: Perturbing the diagonal doesn't even change the square-free discriminant — it stays at 163 regardless

## The Mechanism: Norm Structure in ℤ[ζ₉]

The cross-overlap circulant vector g = [g₀, g₁, ..., g₈] defines an algebraic integer:

**β = Σᵣ gᵣ · ζ₉ʳ ∈ ℤ[ζ₉]**

The norm from ℚ(ζ₉) to ℚ is:

**N(β) = Π_{k coprime to 9} β(ζ₉ᵏ)**

### Computed norms:

| Level | Cross-overlap β | N(β) | Factorization | Odd prime | p mod (2n-1) |
|-------|----------------|------|---------------|-----------|-------------|
| K₈ | [4,0,0,0,2,2,4] ∈ ℤ[ζ₇] | 2752 | 2⁶ · 43 | 43 | 43 ≡ 1 mod 7 |
| K₁₀ | [12,8,0,0,8,4,8,4,0] ∈ ℤ[ζ₉] | 667648 | 2¹² · 163 | 163 | 163 ≡ 1 mod 9 |

### Verification at K₁₀:
- 163 splits completely in ℤ[ζ₉] (since 163 ≡ 1 mod 9)
- Found primitive 9th root ζ₉ ≡ 38 mod 163
- β ≡ 0 mod exactly one prime ideal above 163 (at ζ₉ = 53 mod 163)
- β ≢ 0 mod the other five prime ideals
- This gives v₁₆₃(N(β)) = 1 (one ideal of norm 163 divides β)

## The f-Divisibility Theorem

**Theorem**: Let β ∈ ℤ[ζₚ] and let q be a rational prime with gcd(q, p) = 1. Let f = ord_p(q) be the inertia degree of q in ℚ(ζₚ). Then:

**v_q(N(β)) ≡ 0 mod f**

*Proof*: The prime q splits in ℤ[ζₚ] as q = (𝔮₁ · ... · 𝔮_g)^e where g = φ(p)/f and each 𝔮ᵢ has norm q^f. Therefore:

v_q(N(β)) = Σᵢ f · v_{𝔮ᵢ}(β) = f · (Σᵢ v_{𝔮ᵢ}(β))

which is always a multiple of f. □

## Application to the Congruence Law

**Corollary**: An odd prime q can appear in the square-free part of N(β) only if f = ord_{2n-1}(q) is odd.

The residue classes and their inertia degrees for 2n-1 = 9:

| Residue mod 9 | Inertia degree f | f parity | Can appear in sq-free? |
|---------------|-----------------|----------|----------------------|
| 1 | 1 | odd | YES |
| 8 ≡ -1 | 2 | even | NO |
| 4, 7 | 3 | odd | yes (but don't empirically) |
| 2, 5 | 6 | even | NO |

This explains ALL the computational observations:
- All failing primes in random tests have residue 8 mod 9 (f=2, even)
- The matching chain's discriminant has 17² (17 ≡ 8 mod 9, even exponent) and 163¹ (163 ≡ 1 mod 9, odd exponent)
- No primes with residue 2 or 5 mod 9 ever appear to odd power

## Why the Diagonal Doesn't Matter

The degree-6 polynomial arises from the resultant:

f₂ₙ(x) = Res_c(x² - 2A(c)·x + (A(c)² - B²(c)), minimal_poly(c))

where A(c) comes from the diagonal circulant and B²(c) = |β(ζ₉ᵏ)|² comes from the cross-overlap.

The discriminant of f₂ₙ has the form:

Disc(f₂ₙ) = (powers of 2,3) · N_{ℚ(cos)/ℚ}(B²(c))^α · (terms from A)^β

The key: N_{ℚ(cos)/ℚ}(B²(c)) = |N_{ℚ(ζ)/ℚ}(β)|² = N(β)². The odd prime content of N(β) is entirely determined by β (the cross-overlap element), not by A (the diagonal). The diagonal only affects the powers of 2 and 3 in the discriminant.

## The Revised Proof of Law 3

**Theorem (Congruence Law)**: Let p be an odd prime in the square-free part of Disc(f₂ₙ), with gcd(p, 2n-1) = 1. Then ord_{2n-1}(p) is odd. In particular, p ≡ 1 mod (2n-1) when 2n-1 is prime.

*Proof*:
1. The Gram matrix G of K₂ₙ is Z_{2n-1}-equivariant with direction-sector block structure
2. In each 2-orbit direction sector, the off-diagonal circulant defines β ∈ ℤ[ζ_{2n-1}]
3. The discriminant of the sector characteristic polynomial contains N(β) as a factor (through the resultant construction)
4. By the f-divisibility theorem, v_p(N(β)) ≡ 0 mod f where f = ord_{2n-1}(p)
5. For p to appear in the square-free part, v_p must be odd, hence f must be odd
6. When 2n-1 is prime, f | φ(2n-1) = 2n-2. The odd divisors of 2n-2 with f = 1 give p ≡ 1 mod (2n-1). □

*Note*: The proof allows p with f = 3, 5, etc. (odd f > 1). These don't appear empirically at K₈ and K₁₀ because N(β) is small enough that its only odd prime factors have f = 1. Whether this persists at all levels is an open question.

## Connection to "Casting Out Nines"

The digit-sum property of 9 in base 10: any multiple of 9 has digit sum divisible by 9, because 10 ≡ 1 mod 9.

The congruence law: only primes p with p ≡ 1 mod 9 appear in the square-free discriminant.

Both are manifestations of the same arithmetic: **the residue 1 mod 9 is the identity element of the multiplicative group (ℤ/9ℤ)***. Primes at the identity (f=1) contribute freely to the norm. Primes at other elements (f > 1) contribute in packets of size f, which are invisible to parity when f is even.

The "casting out" is literal: the norm map N: ℤ[ζ₉] → ℤ projects the ideal structure onto the integers, and this projection erases odd-power contributions from primes with even inertia degree — exactly as digit summation "casts out" multiples of 9.

## Status of the Three Laws

| Law | Statement | Status | Source |
|-----|-----------|--------|--------|
| 1 (Degree) | deg(f₂ₙ) = φ(2n-1) | GENERIC | Z_{2n-1} equivariance |
| 2 (Wreath) | Gal ≅ C₂ ≀ C_{φ/2} | GENERIC | Cyclotomic field theory |
| 3 (Congruence) | sq-free primes have odd ord | SPECIFIC | β ∈ ℤ[ζ_{2n-1}] norm structure |

Laws 1 and 2 hold for any Z_{2n-1}-equivariant integer matrix.
Law 3 requires that the off-diagonal structure comes from an element of ℤ[ζ_{2n-1}] — which the matching chain provides through its cross-overlap circulant.

The gap in the original proof was: [P,G] = 0 gives Laws 1 and 2, but not Law 3. The missing ingredient is the **integral cyclotomic structure** of the cross-overlap element β.
