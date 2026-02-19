---
title: Six-Level Motivic Chain: Complete Arithmetic Structure
section: Program Overview
status: active
---

# Six-Level Motivic Chain: Complete Arithmetic Structure

## Brian Porter â€” February 2026

---

## Summary of New Results

The spectral extractor has been computed through six levels (Kâ‚†, Kâ‚ˆ, Kâ‚â‚€, Kâ‚â‚‚, Kâ‚â‚„, Kâ‚â‚†). This document collects all arithmetic invariants in one place, reports three new findings, and corrects a claim in the continuum limit analysis.

**New findings:**

1. **Law 4 (Constant Term Law):** The constant term of every vacuum polynomial factors as 2^{deg} Ã— (odd part), with the odd part governed by the same congruence arithmetic as the discriminant.

2. **Complete discriminant factorization:** All six discriminants now fully factored through their square-free parts. The congruence law (square-free primes p â‰¡ 1 mod (2nâˆ’1)) holds at every level without exception.

3. **Correction to 2/Ï€ conjecture:** The raw vacuum eigenvalues Î»_min converge to ZERO as ~n^{âˆ’2.22}, not to 2/Ï€. The Kâ‚â‚‚ match was a numerical coincidence on three data points. The values 0.672 (Kâ‚ˆ) and 0.659 (Kâ‚â‚€) reported in the continuum limit analysis use an undocumented normalization inconsistent with the raw sector eigenvalues.

---

## 1. The Vacuum Polynomials

All six minimal polynomials of the vacuum eigenvalue, with integer coefficients throughout:

**Kâ‚†** (n=3, q=5, deg=2):
fâ‚†(x) = xÂ² âˆ’ 10x + 20

**Kâ‚ˆ** (n=4, q=7, deg=6):
fâ‚ˆ(x) = xâ¶ âˆ’ 44xâµ + 720xâ´ âˆ’ 5648xÂ³ + 22512xÂ² âˆ’ 43456x + 31808

**Kâ‚â‚€** (n=5, q=9, deg=6):
fâ‚â‚€(x) = xâ¶ âˆ’ 48xâµ + 828xâ´ âˆ’ 6560xÂ³ + 24624xÂ² âˆ’ 40320x + 21312

**Kâ‚â‚‚** (n=6, q=11, deg=10):
fâ‚â‚‚(x) = xÂ¹â° âˆ’ 96xâ¹ + 3616xâ¸ âˆ’ 70976xâ· + 803584xâ¶ âˆ’ 5466944xâµ + 22570944xâ´ âˆ’ 55624448xÂ³ + 77834240xÂ² âˆ’ 55385088x + 14879744

**Kâ‚â‚„** (n=7, q=13, deg=12):
fâ‚â‚„(x) = xÂ¹Â² âˆ’ 128xÂ¹Â¹ + 6272xÂ¹â° âˆ’ 162480xâ¹ + 2518144xâ¸ âˆ’ 24714048xâ· + 157488384xâ¶ âˆ’ 655492864xâµ + 1764928256xâ´ âˆ’ 2987878400xÂ³ + 2995546112xÂ² âˆ’ 1568845824x + 311554048

**Kâ‚â‚†** (n=8, q=15, deg=8):
fâ‚â‚†(x) = xâ¸ âˆ’ 120xâ· + 3984xâ¶ âˆ’ 56640xâµ + 382496xâ´ âˆ’ 1230720xÂ³ + 1808064xÂ² âˆ’ 1113600x + 215296

All polynomials are **monic, integer-coefficient, and irreducible over â„š**.

---

## 2. The Four Laws

### Law 1: Cyclotomic Degree

deg(f_{2n}) = Ï†(2nâˆ’1), with the degenerate case Kâ‚† (where Ï†(5) = 4 but deg = 2).

| Level | 2nâˆ’1 | Ï†(2nâˆ’1) | deg(f) | Note |
|-------|-------|---------|--------|------|
| Kâ‚†   | 5     | 4       | 2      | Degenerate (real subfield) |
| Kâ‚ˆ   | 7     | 6       | 6      | âœ“ |
| Kâ‚â‚€  | 9     | 6       | 6      | âœ“ |
| Kâ‚â‚‚  | 11    | 10      | 10     | âœ“ |
| Kâ‚â‚„  | 13    | 12      | 12     | âœ“ |
| Kâ‚â‚†  | 15    | 8       | 8      | âœ“ (first degree DROP: Ï†(15) < Ï†(13)) |

### Law 2: Wreath Product

Gal(K_n/â„š) â‰… Câ‚‚ â‰€ C_{d/2} where d = deg(f_{2n}).

| Level | Galois group | Order |
|-------|-------------|-------|
| Kâ‚†   | Câ‚‚          | 2     |
| Kâ‚ˆ   | Câ‚‚ â‰€ Câ‚ƒ     | 24    |
| Kâ‚â‚€  | Câ‚‚ â‰€ Câ‚ƒ     | 24    |
| Kâ‚â‚‚  | Câ‚‚ â‰€ Câ‚…     | 160   |
| Kâ‚â‚„  | Câ‚‚ â‰€ Câ‚†     | 384   |
| Kâ‚â‚†  | Câ‚‚ â‰€ Câ‚„     | 64    |

### Law 3: Congruence (Discriminant)

Every odd prime p appearing with ODD exponent in disc(f_{2n}) â€” that is, every prime in the square-free part of the discriminant â€” satisfies p â‰¡ 1 mod (2nâˆ’1).

Complete factorization of all six discriminants:

| Level | q  | disc(f) factored | Square-free part | Primes, residues mod q |
|-------|-----|-----------------|-----------------|----------------------|
| Kâ‚†   | 5  | 2Â² Â· 5           | 5               | 5 = q [level prime] |
| Kâ‚ˆ   | 7  | 2Â³â¶ Â· 7â´ Â· 43 Â· 421Â² | 43         | 43 â‰¡ 1 mod 7 âœ“ |
| Kâ‚â‚€  | 9  | 2â´â¸ Â· 3Â¹â° Â· 17Â² Â· 163 | 163        | 163 â‰¡ 1 mod 9 âœ“ |
| Kâ‚â‚‚  | 11 | 2Â¹â°â° Â· 11â¸ Â· 23 Â· 353 Â· 439Â² Â· 64217Â² Â· 478838249Â² | 23 Ã— 353 | 23 â‰¡ 1, 353 â‰¡ 1 mod 11 âœ“ |
| Kâ‚â‚„  | 13 | 2Â¹âµâ¶ Â· 5Â² Â· 13Â¹â° Â· 31Â² Â· 79 Â· 5279 Â· 32941Â² Â· 414259Â² Â· 1256726743Â² | 79 Ã— 5279 | 79 â‰¡ 1, 5279 â‰¡ 1 mod 13 âœ“ |
| Kâ‚â‚†  | 15 | 2â¶â´ Â· 3â¸ Â· 5â¶ Â· 11Â² Â· 29Â² Â· 31 Â· 2371 Â· 2069873161Â² | 31 Ã— 2371 | 31 â‰¡ 1, 2371 â‰¡ 1 mod 15 âœ“ |

**Heegner connection:** Kâ‚ˆ and Kâ‚â‚€ have prime discriminant square-free parts (43, 163) which are both Heegner numbers. From Kâ‚â‚‚ onward, the square-free part has two prime factors.

### Law 4: Constant Term (NEW)

The constant term of f_{2n}(x) factors as:

**f_{2n}(0) = 2^{deg(f)} Ã— (odd part)**

where the 2-adic valuation equals the degree EXACTLY. The odd part carries the same congruence structure as the discriminant:

| Level | q  | deg | Constant | vâ‚‚ | Odd part | Odd factored |
|-------|-----|-----|----------|-----|----------|-------------|
| Kâ‚†   | 5  | 2   | 20       | 2   | 5        | 5            |
| Kâ‚ˆ   | 7  | 6   | 31808    | 6   | 497      | 7 Ã— 71      |
| Kâ‚â‚€  | 9  | 6   | 21312    | 6   | 333      | 3Â² Ã— 37     |
| Kâ‚â‚‚  | 11 | 10  | 14879744 | 10  | 14531    | 11 Ã— 1321   |
| Kâ‚â‚„  | 13 | 12  | 311554048| 12  | 76063    | 13 Ã— 5851   |
| Kâ‚â‚†  | 15 | 8   | 215296   | 8   | 841      | 29Â²         |

Structure of the odd part:

- Kâ‚† through Kâ‚â‚„ (q prime or q = 9): odd part = q Ã— pâ‚‚ where pâ‚‚ â‰¡ 1 mod q.
  Specifically: pâ‚‚ = 1 (Kâ‚†), 71 (Kâ‚ˆ), 37 (Kâ‚â‚€), 1321 (Kâ‚â‚‚), 5851 (Kâ‚â‚„).

- Kâ‚â‚† (q = 15, first composite): odd part = 29Â² with ordâ‚â‚…(29) = 2, e = 2.
  The f-divisibility condition e â‰¡ 0 mod ord_q(p) is satisfied: 2/2 = 1.

**Interpretation:** The constant term equals the norm N_{â„š(Î»)/â„š}(Î»_vac), i.e., the product of all Galois conjugates of the vacuum eigenvalue. The exact power vâ‚‚ = deg shows that 2 ramifies in the number field K_n with a specific, predictable pattern.

---

## 3. Galois Disjointness

All 15 pairs of number fields are Galois disjoint, verified by non-vanishing resultants:

Res(f_i, f_j) â‰  0 for all i â‰  j.

The compositum through Kâ‚â‚† has degree:
2 Ã— 6 Ã— 6 Ã— 10 Ã— 12 Ã— 8 = **69,120**

This is the product of all six degrees â€” maximal, with no collapse.

---

## 4. Spectral Data

### Vacuum eigenvalues (smallest root of each f_{2n})

| Level | n | Î»_min | Î»_max | Î»_min/Î»_max |
|-------|---|-------|-------|-------------|
| Kâ‚†   | 3 | 2.76393 | 7.23607 | 0.38197 |
| Kâ‚ˆ   | 4 | 1.95951 | 17.8001 | 0.11008 |
| Kâ‚â‚€  | 5 | 0.97929 | 20.9487 | 0.04675 |
| Kâ‚â‚‚  | 6 | 0.63663 | 36.0272 | 0.01767 |
| Kâ‚â‚„  | 7 | 0.45886 | 54.7809 | 0.00838 |
| Kâ‚â‚†  | 8 | 0.34877 | 77.0171 | 0.00453 |

### Asymptotic behavior

Best-fit power law: **Î»_min(n) â‰ˆ 35.5 Ã— n^{âˆ’2.22}**

The vacuum eigenvalue converges to ZERO, not to 2/Ï€ â‰ˆ 0.63662. The Kâ‚â‚‚ numerical coincidence (Î»_min = 0.63663 â‰ˆ 2/Ï€ to 4 digits) does not persist: Kâ‚â‚„ and Kâ‚â‚† continue decreasing.

The ratio Î»_min/Î»_max also converges to zero, showing the spectral gap widens relative to the spectral radius.

### Power sums (all exact integers)

| Level | pâ‚ (trace) | pâ‚‚ | pâ‚ƒ | pâ‚„ |
|-------|-----------|-----|-----|-----|
| Kâ‚†   | 10  | 60      | 400      | 2800      |
| Kâ‚ˆ   | 44  | 496     | 7088     | 113216    |
| Kâ‚â‚€  | 48  | 648     | 11040    | 209760    |
| Kâ‚â‚‚  | 96  | 1984    | 56256    | 1825792   |
| Kâ‚â‚„  | 128 | 3840    | 176144   | 9186816   |
| Kâ‚â‚†  | 120 | 6432    | 463680   | 35283328  |

### R-ratio test

The ratio R(n) = eâ‚„/eâ‚‚Â² does NOT converge:
Kâ‚ˆ: 0.0434, Kâ‚â‚€: 0.0359, Kâ‚â‚‚: 0.0615, Kâ‚â‚„: 0.0640, Kâ‚â‚†: 0.0241.

No stable pattern. The ratio pâ‚„/pâ‚‚Â² also oscillates: 0.778, 0.460, 0.500, 0.464, 0.623, 0.853. The R-ratio test is INCONCLUSIVE â€” neither convergent nor divergent in an obvious way.

---

## 5. Correction to Continuum Limit Analysis

The earlier document "The Continuum Limit: Analytic Convergence vs Algebraic Divergence" reported:

- Kâ‚ˆ: Î»_vac â‰ˆ 0.672
- Kâ‚â‚€: Î»_vac â‰ˆ 0.659
- Kâ‚â‚‚: Î»_vac â‰ˆ 0.6366

and conjectured convergence to 2/Ï€.

**Status:** The Kâ‚â‚‚ value (0.6366) matches the raw sector eigenvalue. The Kâ‚ˆ value (0.672) and Kâ‚â‚€ value (0.659) do NOT match their raw eigenvalues (1.9595, 0.9793 respectively). The normalization that produced those two values is undocumented.

With six levels now computed, the raw eigenvalue sequence (2.764, 1.960, 0.979, 0.637, 0.459, 0.349) is unambiguously decreasing toward zero. The conjecture "Î»_vac â†’ 2/Ï€" does not hold for the raw eigenvalues.

**Open question:** Is there a NATURAL normalization of the vacuum eigenvalue that converges to a nonzero limit? The quantity Î»_min Ã— n^{2.22} â‰ˆ 35 is approximately constant, but the exponent 2.22 is empirical. A theoretical derivation of the decay rate (or identification of a convergent normalized quantity) remains open.

---

## 6. Johnson Scheme Budget Identity

Confirmed at all six levels:

Î»_mid Ã— d_phys = 2(nâˆ’1) Ã— Î»_max

where d_phys = C(2nâˆ’1, nâˆ’1)/(2nâˆ’1) is the physical sector dimension and Î»_max = C(2nâˆ’2, nâˆ’1) is the largest Johnson eigenvalue.

| Level | d_phys | Î»_max_J | Budget |
|-------|--------|---------|--------|
| Kâ‚†   | 2      | 6       | 24     |
| Kâ‚ˆ   | 5      | 20      | 120    |
| Kâ‚â‚€  | 14     | 70      | 560    |
| Kâ‚â‚‚  | 42     | 252     | 2520   |
| Kâ‚â‚„  | 132    | 924     | 11088  |
| Kâ‚â‚†  | 429    | 3432    | 48048  |

---

## 7. Open Problems

### Resolved
1. âœ… Galois disjointness through Kâ‚â‚† (all 15 pairs)
2. âœ… Discriminant congruence at all 6 levels
3. âœ… Constant term structure (Law 4)
4. âœ… Degree drop at Kâ‚â‚† (Ï†(15) = 8 < Ï†(13) = 12)

### Still Open
1. **Period structure:** Can Î»_vac be expressed as a period (integral of algebraic form over algebraic cycle)?
2. **Decay rate:** Why Î»_min ~ n^{âˆ’2.22}? Is the exponent rational?
3. **Kâ‚â‚ˆ computation:** 2nâˆ’1 = 17 (prime), predicted deg = 16. Would test Law 1 at a new prime level and extend the chain to 7 levels.
4. **Normalization identification:** What quantity (if any) from the spectral extractor converges to a nonzero transcendental limit?
5. **Cross-structure test:** Apply the extractor to Eâ‚ˆ root system, Young tableaux, or other representation-theoretic sources to test universality.
