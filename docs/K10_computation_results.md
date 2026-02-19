---
title: K₁₀ Computation: Degree Stabilization and the Fourth Scenario
section: Higher Levels
status: active
---

# K₁₀ Computation: Degree Stabilization and the Fourth Scenario

## The Arithmetic Chain K₄ → K₆ → K₈ → K₁₀

Brian Porter — February 2026

---

## 0. Main Result

The K₁₀ vacuum eigenvalue satisfies an **irreducible degree-6 polynomial** over ℚ, with Galois group **C₂ ≀ C₃** — identical in degree and Galois structure to K₈. The degree sequence 1 → 2 → 6 → 6 has **stabilized**. The Galois group sequence trivial → C₂ → C₂≀C₃ → C₂≀C₃ has stabilized. Yet the number fields at each level remain **pairwise Galois disjoint**: they share no algebraic structure beyond ℚ. 

This combination — **stabilized algebraic structure with divergent arithmetic content** — is a fourth scenario not anticipated in the prior analysis. It reveals that the polynomial/series wall is about **field disjointness**, not algebraic complexity.

---

## 1. Computation

### 1.1 Setup

K₁₀ has 945 = 9!! perfect matchings. The Z₉-symmetric direction/phase structure was constructed following the paperIII prescription:

- **Edge classes**: K₉ edges partitioned into 4 classes of 9 by Z₉ distance (plus 9 handle edges to vertex 9)
- **Phase assignment**: ω^(class mod 3) where ω = e^(2πi/3)
- **Direction vectors**: d₀ = (1,0,0), d₁ = (0,1,0), d₂ = (0,0,1), d₃ = (-1,-1,-1)
- **Gram matrix**: G_{ij} = Tr(M_i M_j) · Re(ζ_j* ζ_i) · δ(D_i, D_j)

All Johnson scheme predictions verified exactly: λ_max = 1050, λ_mid = 240, d_phys = 35, budget identity 240 × 35 = 8400 = 2 × 4 × 1050.

### 1.2 Vacuum Eigenvalue

The 945×945 Gram matrix G has vacuum eigenvalue:

$$\lambda_{\text{vac}}(K_{10}) \approx 0.97929166601$$

with multiplicity 6, living in Z₉ sectors {1,2,4,5,7,8} (three conjugate pairs).

### 1.3 Direction Sector Reduction

The vacuum eigenvectors span three direction sectors: (2,−1,−1), (1,3,0), (−3,−2,−3), each containing 18 matchings. Within each 18-matching sector, the Gram matrix is an **18×18 integer matrix** (after multiplying G by 2).

### 1.4 Minimal Polynomial

The characteristic polynomial of the 18×18 integer block factors completely:

$$\det(xI - 2G|_{\text{sector}}) = x(x - 88)(x - 36)^2(x - 4)^2 \cdot [x^6 - 96x^5 + 3312x^4 - 52480x^3 + 393984x^2 - 1290240x + 1363968]^2$$

The vacuum eigenvalue (of 2G) satisfies the **irreducible sextic**:

$$f_{10}(x) = x^6 - 96x^5 + 3312x^4 - 52480x^3 + 393984x^2 - 1290240x + 1363968$$

Equivalently, the vacuum of G satisfies:

$$g_{10}(y) = y^6 - 48y^5 + 828y^4 - 6560y^3 + 24624y^2 - 40320y + 21312$$

---

## 2. Galois Group

The Galois group of f₁₀(x) over ℚ is:

$$\text{Gal}(f_{10}/\mathbb{Q}) \cong C_2 \wr C_3 = (C_2)^3 \rtimes C_3$$

of order 24. This is **identical** to the K₈ Galois group.

Structural properties:
- Order: 24
- Abelian: No
- Solvable: Yes
- Primitive: No (imprimitive on 6 points)
- Derived subgroup: V₄ = C₂ × C₂ (order 4)
- Abelianization: order 6
- Not contained in A₆ (discriminant not a perfect square)

The six roots form three blocks of two, cyclically permuted by the C₃ factor:

$$\{1.959, 17.785\}, \quad \{5.138, 41.897\}, \quad \{8.993, 20.229\}$$

Block sums: 19.743, 47.035, 29.222 — roots of the cubic resolvent t³ − 96t² + 2880t − 27136 = 0.

---

## 3. Discriminant

$$\Delta(f_{10}) = 2^{78} \times 3^{10} \times 17^2 \times 163$$

**Square-free part: 163** (the largest Heegner number).

| Level | Square-free disc | New primes | Structural prime (2n−1) |
|:---|:---|:---|:---|
| K₆ | 5 | {5} | 5 (prime) |
| K₈ | 43 | {7, 43, 421} | 7 (prime) |
| K₁₀ | 163 | {3, 17, 163} | 9 = 3² (composite) |

The square-free discriminants 5, 43, 163 are **distinct primes** at each level. Each level introduces genuinely new arithmetic content.

---

## 4. Galois Disjointness

### 4.1 Results

All three pairs of number fields are **Galois disjoint over ℚ**:

- **ℚ(λ₆) ∩ ℚ(λ₈) = ℚ**: proved previously (ℚ(λ₈) has no quadratic subfield)
- **ℚ(λ₈) ∩ ℚ(λ₁₀) = ℚ**: cubic subfields are distinct (GCD = 1, resultant = 2¹⁸ × 23689 ≠ 0)
- **ℚ(λ₆) ∩ ℚ(λ₁₀) = ℚ**: ℚ(λ₁₀) has no quadratic subfield (C₂ ≀ C₃ lattice)

The compositum has maximal degree:

$$[\mathbb{Q}(\lambda_6, \lambda_8, \lambda_{10}) : \mathbb{Q}] = 2 \times 6 \times 6 = 72$$

### 4.2 Cubic Subfield Comparison

| Level | Cubic resolvent | Cubic disc |
|:---|:---|:---|
| K₈ | t³ − 44t² + 608t − 2624 | 2¹² × 7² |
| K₁₀ | t³ − 96t² + 2880t − 27136 | 2¹⁸ × 3⁴ |

Different discriminants, coprime over ℚ → the cubic subfields are **algebraically independent**.

---

## 5. The Fourth Scenario

### 5.1 What Was Expected

Three scenarios were considered for the limit n → ∞:

- **A**: Degree stabilizes, fields eventually embed → limit algebraic
- **B**: Degree grows without bound → limit transcendental
- **C**: No convergence

### 5.2 What We Found

**Scenario D**: Degree stabilizes at 6, Galois group stabilizes at C₂ ≀ C₃, but fields remain **pairwise disjoint**. Each level generates a structurally isomorphic but arithmetically independent extension.

This means:
1. The **type** of algebraic number is fixed (degree 6, C₂ ≀ C₃ structure, three-block partition)
2. The **specific** algebraic number changes at each level (different primes, different cubic subfields)
3. The limit, if it exists, is **not** in any of these degree-6 fields → **transcendental**
4. But it's "uniformly structured" — each approximation has identical Galois type

### 5.3 Implications for the Wall

The polynomial/series wall is not about **algebraic complexity** (degree growth). It's about **arithmetic disjointness** (independent field extensions). You can have polynomials of FIXED degree whose roots generate INDEPENDENT fields — and the sequence of roots converges to a transcendental limit.

This is the matching chain's precise version of the RH wall: the Frobenius eigenvalues at each prime p are algebraic numbers of bounded type (they lie on |α| = q^{1/2}), but the fields they generate are independent. The limit (the zeros of L(s,χ)) is transcendental — not because the Frobenius eigenvalues become more complex, but because they become more **numerous and independent**.

---

## 6. The Full Comparison Table

| Property | K₄ | K₆ | K₈ | K₁₀ |
|:---|:---|:---|:---|:---|
| λ_vac | 1 | 5−√5 ≈ 2.764 | ≈ 1.960 | ≈ 0.979 |
| Degree over ℚ | 1 | 2 | **6** | **6** |
| Galois group | {1} | C₂ | **C₂≀C₃** | **C₂≀C₃** |
| |Gal| | 1 | 2 | 24 | 24 |
| Sq-free disc | — | 5 | 43 | 163 |
| Disc primes | — | {5} | {2,7,43,421} | {2,3,17,163} |
| Cubic subfield | — | — | t³−44t²+608t−2624 | t³−96t²+2880t−27136 |
| Quadratic sub? | — | ℚ(√5) | No | No |
| Z₉/Z₇ struct | — | Z₅ | Z₇ | Z₉ |
| Multiplicity | 1 | 2 | 6 | 6 |
| d_phys | 2 | 9 | 20 | 35 |

---

## 7. Discriminant Primes: The Arithmetic Fingerprint

The square-free discriminant primes 5, 43, 163 are all **distinct primes** and follow no obvious pattern. However:

- **5** = structural prime of K₆ (Z₅ symmetry)
- **43** = spectral prime of K₈ (no combinatorial origin)
- **163** = spectral prime of K₁₀ (and the largest Heegner number)

The appearance of 163 — the largest integer d such that ℚ(√(−d)) has class number 1 — is striking but may be coincidental. It would be extremely interesting to see whether K₁₂'s square-free discriminant continues the pattern.

The **structural primes** (from the cyclic symmetry 2n−1):
- K₆: 2n−1 = 5 (prime)
- K₈: 2n−1 = 7 (prime)  
- K₁₀: 2n−1 = 9 = 3² (composite — first non-prime)

K₁₀ is the first level where 2n−1 is composite. The discriminant contains 3 (from 9 = 3²) as a new prime, reflecting this.

---

## 8. What This Means for the Program

### 8.1 The Wall Is Real But Benign

The degree/Galois stabilization means the physics at each level has **identical algebraic structure**: a degree-6 extension with C₂ ≀ C₃ Galois group and a three-block partition. The three blocks are the three generations. The C₂ factors are within-generation conjugation. The C₃ factor is generation permutation.

The wall blocks only the LIMITING value, which is transcendental. But:
- The algebraic structure (degree 6, C₂≀C₃, three generations) is settled by K₈ and confirmed at K₁₀
- The numerical convergence is guaranteed by self-focusing (1/n² corrections)
- By K₁₀, η = 8/9 ≈ 88.9% of spectral weight is captured in the physical sector

### 8.2 The Three-Generation Structure Is Galois-Theoretic

The Galois group C₂ ≀ C₃ has three blocks — this IS the three-generation structure. It appeared at K₈ (where the Yukawa matrix was found) and **persists at K₁₀**. The stabilization means three generations is not an accident of K₈ but a persistent feature of the matching chain's arithmetic.

### 8.3 What to Compute Next

1. **K₁₂**: Does the degree remain 6? Does the Galois group remain C₂ ≀ C₃? What is the new square-free discriminant?
2. **Cubic subfield pattern**: The cubic resolvents at K₈ and K₁₀ have different discriminants. Does the cubic discriminant follow a pattern?
3. **R ratio**: Compute R(K₁₀) = a₄/a₂² restricted to V_phys. Does it converge toward R(K₈)?
4. **Period test**: With the degree-6 minimal polynomial known, compute λ_vac(K₁₀) to 100+ digits and run PSLQ against period basis.

---

## 9. Caveats

The K₁₀ computation uses the **natural Z₉-symmetric direction assignment** (edge classes by distance mod 9), not a direction assignment derived from a specific genus-3 surface embedding of K₉. The genus-3 embedding may give a different direction/phase structure and hence a different Gram matrix. The qualitative results (degree stabilization, Galois group stabilization) should be robust to the choice of embedding, but the specific numerical values and discriminant primes may change.

The K₈ computation in paperIII used the Heawood embedding of K₇ on the torus, which is unique. K₉'s genus-3 embedding is not unique, and the "canonical" choice requires further analysis. This is the primary remaining task for making the K₁₀ computation fully rigorous.
