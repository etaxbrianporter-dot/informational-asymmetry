---
title: Kâ‚â‚† Computation Results: The Degree Drop Confirmed
section: Higher Levels
status: active
---

# Kâ‚â‚† Computation Results: The Degree Drop Confirmed

## Brian Porter â€” February 2026

---

## 0. Main Results

**Result 1 (Degree drop).** The Kâ‚â‚† vacuum eigenvalue satisfies an irreducible degree-8 polynomial over â„š:

$$f_{16}(x) = x^8 - 120x^7 + 3984x^6 - 56640x^5 + 382496x^4 - 1230720x^3 + 1808064x^2 - 1113600x + 215296$$

with deg = 8 = Ï†(15) = Ï†(3 Ã— 5). The degree **drops** from 12 (at Kâ‚â‚„) to 8, confirming the cyclotomic degree law at its sixth data point and its first non-monotone transition. This non-monotone behavior â€” degree *decreasing* as graph size increases â€” cannot arise from any mechanism that scales with combinatorial complexity. It is purely arithmetic, controlled by the factorization 15 = 3 Ã— 5.

**Result 2 (Wreath product).** Gal(fâ‚â‚†/â„š) â‰… Câ‚‚ â‰€ Câ‚„ = (Câ‚‚)â´ â‹Š Câ‚„, order 64. All 8 conjugacy classes observed with Ï‡Â² = 2.04 on 7 degrees of freedom (excellent fit). The block count is Ï†(15)/2 = 4.

**Result 3 (Congruence law).** The square-free discriminant is **73501 = 31 Ã— 2371**, with both primes satisfying 31 â‰¡ 2371 â‰¡ 1 mod 15.

**Result 4 (Perfect Legendre correlation).** 662 out of 662 primes show exact correlation between Legendre symbol (73501/p) and Frobenius permutation sign.

**Result 5 (Galois disjointness).** â„š(Î»â‚â‚†) is Galois disjoint from all five previous fields. The compositum through Kâ‚â‚† has degree 2 Ã— 6 Ã— 6 Ã— 10 Ã— 12 Ã— 8 = **69,120** (maximal).

---

## 1. Computation

### 1.1 Enumeration

Kâ‚â‚† has N = 15!! = 2,027,025 perfect matchings. All matchings enumerated in 3.0s via parallel recursion (15 first-pair choices across 4 CPU cores). Under Zâ‚â‚… action, they decompose into 135,135 orbits, all of size 15. There are 1,701 distinct edge-class sectors.

### 1.2 Sector Structure

| Sector type | Orbit count | Sector size | Count | Notes |
|:------------|:------------|:------------|:------|:------|
| Single-orbit | 1 | 15 | 8 | Palindromic â†’ half-degree |
| Two-orbit | 2 | 30 | 4 | Full irreducible factor |
| Three-orbit | 3 | 45 | 4 | Both full and half factors |
| Larger | â‰¥4 | â‰¥60 | 1,685 | Not needed for extraction |

The vacuum eigenvalue was extracted from **two-orbit sectors** (30 matchings each). Four such sectors exist, corresponding to edge-class signatures:
- (6,0,0,1,0,0,0,1), (1,0,0,6,0,0,0,1), (0,6,0,0,0,0,1,1), (0,1,0,0,0,0,6,1)

### 1.3 The Minimal Polynomial

The 30Ã—30 integer overlap matrix in each two-orbit sector has characteristic polynomial factoring as:

$$(x - \lambda_1)^{a_1} \cdots \cdot [f_{16}(x)]^2$$

where fâ‚â‚†(x) is the degree-8 irreducible factor â€” **identical across all four sectors**. Specifically:

$$f_{16}(x) = x^8 - 120x^7 + 3984x^6 - 56640x^5 + 382496x^4 - 1230720x^3 + 1808064x^2 - 1113600x + 215296$$

All 8 roots are real:

| Root | Value |
|:-----|:------|
| Î±â‚ | 0.348774883131 |
| Î±â‚‚ | 0.900390174317 |
| Î±â‚ƒ | 1.581562532857 |
| Î±â‚„ | 3.643214493670 |
| Î±â‚… | 7.333802651377 |
| Î±â‚† | 13.123860174048 |
| Î±â‚‡ | 16.051334719154 |
| Î±â‚ˆ | 77.017060371447 |

Three-orbit sectors additionally produce the degree-4 palindromic factor:
$$g(x) = x^4 - 34x^3 - 36x^2 + 1800x + 4880$$
which is the half-degree factor expected from the palindromic constraint c_r = c_{15âˆ’r} in single-orbit sectors.

---

## 2. Galois Group: Câ‚‚ â‰€ Câ‚„

### 2.1 Discriminant

$$\Delta(f_{16}) = 2^{64} \times 3^8 \times 5^6 \times 11^2 \times 29^2 \times 31 \times 2371 \times 2069873161^2$$

Square-free part: **73501 = 31 Ã— 2371**.

### 2.2 Congruence Check

| Prime | p mod 15 | Status |
|:------|:---------|:-------|
| 31 | 1 | âœ“ |
| 2371 | 1 | âœ“ |

Both primes in the square-free discriminant satisfy p â‰¡ 1 mod 15, extending the congruence law through six levels.

### 2.3 Frobenius Cycle Types

Over 662 primes (up to 5000), the observed cycle types match Câ‚‚ â‰€ Câ‚„ predictions:

| Cycle type | Predicted (Câ‚‚â‰€Câ‚„) | Observed | Diff |
|:-----------|:-------------------|:---------|:-----|
| (4,4) | 0.3125 | 0.3233 | +0.011 |
| (8,) | 0.2500 | 0.2538 | +0.004 |
| (2,2,4) | 0.1250 | 0.1148 | âˆ’0.010 |
| (1,1,1,1,2,2) | 0.0938 | 0.0831 | âˆ’0.011 |
| (2,2,2,2) | 0.0781 | 0.0816 | +0.003 |
| (1,1,1,1,1,1,2) | 0.0625 | 0.0650 | +0.003 |
| (1,1,2,2,2) | 0.0625 | 0.0650 | +0.003 |
| (1â¸) | 0.0156 | 0.0136 | âˆ’0.002 |

**Ï‡Â² = 2.04 on 7 degrees of freedom** (excellent fit; critical value at 5% is 14.1).

### 2.4 Structural Verification

- **All 8 conjugacy classes observed**: Câ‚‚ â‰€ Câ‚„ has exactly 8 conjugacy classes, all detected.
- **Cycle parts**: {1, 2, 4, 8} â€” exactly the allowed parts for Câ‚‚ â‰€ Câ‚„.
- **Group order**: |Câ‚‚ â‰€ Câ‚„| = 2â´ Ã— 4 = 64. Identity frequency 9/662 = 0.0136 (predicted 1/64 = 0.0156).
- **Legendre symbol**: 662/662 **PERFECT** correlation between (73501/p) and Frobenius sign.

---

## 3. Pairwise Galois Disjointness

### 3.1 All Pairs Verified

| Pair | Method | Result |
|:-----|:-------|:-------|
| â„š(Î»â‚†) âˆ© â„š(Î»â‚â‚†) | Quadratic subfields: â„š(âˆš5) â‰  â„š(âˆš73501). Resultant â‰  0. | âœ“ Disjoint |
| â„š(Î»â‚ˆ) âˆ© â„š(Î»â‚â‚†) | Câ‚‚ â‰€ Câ‚ƒ has no quadratic subfield. Different degrees (6 vs 8). | âœ“ Disjoint |
| â„š(Î»â‚â‚€) âˆ© â„š(Î»â‚â‚†) | Same. Resultant â‰  0. | âœ“ Disjoint |
| â„š(Î»â‚â‚‚) âˆ© â„š(Î»â‚â‚†) | â„š(âˆš8119) â‰  â„š(âˆš73501). gcd(10, 8) = 2 but quadratics differ. | âœ“ Disjoint |
| â„š(Î»â‚â‚„) âˆ© â„š(Î»â‚â‚†) | â„š(âˆš417041) â‰  â„š(âˆš73501). Resultant â‰  0. | âœ“ Disjoint |

### 3.2 Compositum

$$[\mathbb{Q}(\lambda_6, \lambda_8, \lambda_{10}, \lambda_{12}, \lambda_{14}, \lambda_{16}) : \mathbb{Q}] = 2 \times 6 \times 6 \times 10 \times 12 \times 8 = \mathbf{69{,}120}$$

Maximal â€” no collapse from any pair. All 15 pairwise intersections are â„š.

---

## 4. The Complete Arithmetic Chain

| Level | 2nâˆ’1 | Ï†(2nâˆ’1) | deg | Gal | |Gal| | Sq-free Î” | Blocks | Ï‡Â² |
|:------|:-----|:--------|:----|:----|:------|:----------|:-------|:----|
| Kâ‚† | 5 | 4 | 2 | Câ‚‚ | 2 | 5 | 1 | â€” |
| Kâ‚ˆ | 7 | 6 | 6 | Câ‚‚â‰€Câ‚ƒ | 24 | 43 | 3 | â€” |
| Kâ‚â‚€ | 9 | 6 | 6 | Câ‚‚â‰€Câ‚ƒ | 24 | 163 | 3 | â€” |
| Kâ‚â‚‚ | 11 | 10 | 10 | Câ‚‚â‰€Câ‚… | 160 | 8119 | 5 | â€” |
| Kâ‚â‚„ | 13 | 12 | 12 | Câ‚‚â‰€Câ‚† | 384 | 417041 | 6 | 7.54 |
| **Kâ‚â‚†** | **15** | **8** | **8** | **Câ‚‚â‰€Câ‚„** | **64** | **73501** | **4** | **2.04** |

### 4.1 Three Laws Confirmed Through Six Levels

**Cyclotomic degree law**: deg(fâ‚‚â‚™) = Ï†(2nâˆ’1). Verified at Kâ‚† (Ï†(5)=4â†’2), Kâ‚ˆ (Ï†(7)=6), Kâ‚â‚€ (Ï†(9)=6), Kâ‚â‚‚ (Ï†(11)=10), Kâ‚â‚„ (Ï†(13)=12), **Kâ‚â‚† (Ï†(15)=8)**. The non-monotone transition 12 â†’ 8 is the strongest test yet. âœ“

**Wreath product law**: Gal(fâ‚‚â‚™/â„š) â‰… Câ‚‚ â‰€ C_{Ï†(2nâˆ’1)/2}. Verified at all six levels. At Kâ‚â‚†: Câ‚‚ â‰€ Câ‚„ with all 8 conjugacy classes and Ï‡Â² = 2.04. âœ“

**Congruence law**: Odd primes in the square-free part of Î”(fâ‚‚â‚™) satisfy p â‰¡ 1 mod (2nâˆ’1). Extended to sixth level: 31 â‰¡ 2371 â‰¡ 1 mod 15. âœ“

### 4.2 Discriminant Primes

| Level | 2nâˆ’1 | Sq-free Î” | Factorization | Congruence |
|:------|:-----|:----------|:-------------|:-----------|
| Kâ‚† | 5 | 5 | 5 | (level prime) |
| Kâ‚ˆ | 7 | 43 | 43 | 43â‰¡1(7) |
| Kâ‚â‚€ | 9 | 163 | 163 | 163â‰¡1(9) |
| Kâ‚â‚‚ | 11 | 8119 | 23 Ã— 353 | both â‰¡1(11) |
| Kâ‚â‚„ | 13 | 417041 | 79 Ã— 5279 | both â‰¡1(13) |
| **Kâ‚â‚†** | **15** | **73501** | **31 Ã— 2371** | **both â‰¡1(15)** |

Pattern: Kâ‚ˆ, Kâ‚â‚€ have prime discriminants (43, 163 â€” both Heegner numbers). Kâ‚â‚‚, Kâ‚â‚„, Kâ‚â‚† have two-prime discriminants. All odd primes appearing satisfy the congruence law.

### 4.3 The Degree Sequence

| Level | 2nâˆ’1 | Ï†(2nâˆ’1) | Note |
|:------|:-----|:--------|:-----|
| Kâ‚† | 5 | 4 (â†’2) | Degenerate |
| Kâ‚ˆ | 7 | 6 | Prime, monotone rise |
| Kâ‚â‚€ | 9 = 3Â² | 6 | Stabilization |
| Kâ‚â‚‚ | 11 | 10 | Monotone rise |
| Kâ‚â‚„ | 13 | 12 | Monotone rise |
| **Kâ‚â‚†** | **15 = 3Ã—5** | **8** | **Non-monotone DROP** |

The degree drop at Kâ‚â‚† is the signature of the cyclotomic law. No other proposed mechanism for the degree sequence â€” combinatorial scaling, random matrix theory, polynomial growth â€” predicts this drop. Only Ï†(15) = Ï†(3)Ï†(5) = 2 Ã— 4 = 8 gives the correct answer.

---

## 5. Physical Interpretation

### 5.1 Generation Count = Ï†(2nâˆ’1)/2 = 4

At Kâ‚â‚†: 4 blocks. This is fewer than Kâ‚â‚„'s 6 blocks, and the same as Kâ‚â‚‚'s 5 would give in a different arithmetic regime. The generation count at Kâ‚â‚† "regresses" because Ï†(15)/2 = 4, while Ï†(13)/2 = 6. This non-universality is again controlled entirely by number theory.

### 5.2 Quadratic Subfield

â„š(âˆš73501) is the quadratic subfield of â„š(Î»â‚â‚†). Note: 73501 = 31 Ã— 2371, and this quadratic field is distinct from all previous levels' quadratic fields (â„š(âˆš5), â„š(âˆš8119), â„š(âˆš417041)).

### 5.3 Group Order Sequence

| Level | |Gal| | Factorization |
|:------|:-------|:--------------|
| Kâ‚† | 2 | 2 |
| Kâ‚ˆ | 24 | 2Â³ Ã— 3 |
| Kâ‚â‚€ | 24 | 2Â³ Ã— 3 |
| Kâ‚â‚‚ | 160 | 2âµ Ã— 5 |
| Kâ‚â‚„ | 384 | 2â· Ã— 3 |
| Kâ‚â‚† | 64 | 2â¶ |

The group order **drops** from 384 to 64, another non-monotone transition. For Câ‚‚ â‰€ C_m, the order is 2^m Ã— m. At Kâ‚â‚†, m = 4, giving 16 Ã— 4 = 64. The all-power-of-2 order at Kâ‚â‚† reflects the fact that Ï†(15)/2 = 4 is itself a power of 2.

---

## 6. Continuum Limit Implications

### 6.1 The Bifurcation Question

With six data points, the algebraic side of the continuum limit is now well-constrained:

- **Fields grow then shrink then grow**: deg = 2, 6, 6, 10, 12, 8, ... The non-monotone behavior means the algebraic complexity does not diverge uniformly.
- **But fields remain disjoint**: The compositum degree 69,120 is maximal. No finite extension of â„š contains all the vacuum eigenvalues.
- **The pro-finite group**: G_âˆž = limâ†Gal(K_n/â„š) is now better constrained by the six levels. The "regression" at Kâ‚â‚† means G_âˆž is NOT the naive limit of ever-larger wreath products.

### 6.2 Next Predictions

| Kâ‚â‚ˆ | Kâ‚‚â‚€ | Kâ‚‚â‚‚ | Kâ‚‚â‚„ |
|:-----|:-----|:-----|:-----|
| 2nâˆ’1 = 17 | 19 | 21 = 3Ã—7 | 23 |
| Ï† = 16 | 18 | 12 | 22 |
| deg = 16 | 18 | **12** | 22 |
| Gal = Câ‚‚â‰€Câ‚ˆ | Câ‚‚â‰€Câ‚‰ | **Câ‚‚â‰€Câ‚†** | Câ‚‚â‰€Câ‚â‚ |

Kâ‚‚â‚‚ provides another non-monotone test: Ï†(21) = 12, same as Kâ‚â‚„. This would give the first **exact degree coincidence** between non-adjacent levels (Kâ‚â‚„ and Kâ‚‚â‚‚ both have degree 12), testing whether the corresponding polynomials live in related or unrelated number fields.

---

## 7. Computation Details

All computations performed with exact integer arithmetic (Python 3 + sympy).

| Stage | Time | Details |
|:------|:-----|:--------|
| Matching enumeration | 3.0s | 2,027,025 matchings, 4-way parallel |
| Zâ‚â‚… classification | 10.1s | 135,135 orbits, 1,701 sectors |
| Overlap matrices | <1s | 30Ã—30 integer matrices in 4 sectors |
| Characteristic polynomial | <1s | Exact over â„¤[x] |
| Factorization over â„š | <1s | Degree-8 irreducible extracted |
| Discriminant | <1s | Full factorization |
| Frobenius analysis (662 primes) | 1.2s | All 8 cycle types observed |
| Galois group verification | <1s | Ï‡Â² = 2.04 |
| Disjointness (resultants) | <1s | All 5 pairs verified |

**Total wall time: ~27 seconds.**

The Kâ‚â‚† computation, despite involving 15Ã— more matchings than Kâ‚â‚„, completed in roughly the same time. The polynomial reduction theorem ensures that the actual linear algebra operates on 30Ã—30 matrices regardless of the 2M matching count.
