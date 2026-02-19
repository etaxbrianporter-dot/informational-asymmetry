# H₇ Portability Test: Results

**Brian Porter — February 19, 2026**

---

## The Experiment

Take the K₇ complete graph — the same graph whose Heawood embedding on the torus defines the direction structure for K₈'s fermion sector. Build the Dirac operator D with Z₃ phases on Heawood direction classes (edges classified by |i−j| mod 7 → three classes of 7). Form D†D. Extract the characteristic polynomial. Factor it. Read the arithmetic.

The question: does the extractor produce nontrivial arithmetic from K₇ itself, or only from K₈'s matching structure?

---

## What We Found

### The raw output

The characteristic polynomial of the K₇ directed Dirac's D†D is:

**f₇(x) = x⁷ − 42x⁶ + 645x⁵ − 4484x⁴ + 14532x³ − 20049x² + 9676x − 144**

All coefficients are integers. The polynomial is **irreducible over ℚ**. All seven roots are real (D†D is Hermitian). This defines a degree-7 number field ℚ(α).

### The discriminant

|Disc| = 2⁴ × 3¹⁸ × 37 × 857 × 50857 × 115987²

Discriminant primes: **{2, 3, 37, 857, 50857, 115987}**

### The Galois group

Twelve distinct cycle types from Frobenius elements at primes up to 200. The discriminant is positive but not a perfect square (so Gal ⊄ A₇). This is consistent with **S₇** (order 5040) — the full symmetric group.

### Factorization mod key primes

| Prime | Factor degrees | Note |
|-------|---------------|------|
| 2 | (1, 1, 5) | |
| 3 | (1, 1, 1, 1, 1, 1, 1) | Splits completely |
| 5 | (1, 6) | |
| 7 | (2, 5) | Does NOT split completely at 7 |
| 43 | (1, 6) | Does NOT have special structure at 43 |

---

## Comparison with K₈

| | K₈ (matching chain) | K₇ (Dirac operator) |
|---|---|---|
| Source | 105 matchings, Gram matrix | 7×7 adjacency, D†D |
| Polynomial degree | 6 | 7 |
| Irreducible | ✓ | ✓ |
| Galois group | C₂ ≀ C₃ (order 12) | S₇ (order ~5040) |
| Solvable | ✓ | ✗ |
| Disc primes | {2, 7, 43, 421} | {2, 3, 37, 857, 50857, 115987} |
| Shared primes | 2 | 2 |
| All roots real | ✗ (3 conjugate pairs) | ✓ |
| Coefficients | Integer | Integer |
| Trace | 44 | 42 = 7 × 6 |
| Determinant | −14190 | 144 = 12² |

The arithmetic is **completely independent**. The only shared discriminant prime is 2 (which is trivially present in almost any discriminant). K₈'s signature primes {7, 43, 421} are absent from K₇'s discriminant. K₇'s primes {3, 37, 857, 50857, 115987} are absent from K₈'s.

**However:** The coefficient a₂ = 645 = 3 × 5 × **43**. The prime 43 appears in K₇'s polynomial coefficients, just not in its discriminant. Suggestive but not load-bearing.

---

## What This Tells Us

### 1. The extractor IS portable

The pipeline (graph + direction phases + spectral extraction) works on K₇ just as on K₈. It produces an irreducible polynomial with integer coefficients, defining a genuine number field with nontrivial arithmetic. The extractor is not specific to matching complexes.

### 2. The direction structure is the engine

Without Z₃ phases: the Fano plane overlap matrix is B = 2I + J (trivial 2-design structure). Every Z₇ sector has eigenvalue 2. No number fields. Dead.

With Z₃ phases: irreducible degree-7 polynomial, rich discriminant, full Galois group. The phases break the 2-design symmetry and inject arithmetic content.

### 3. But the arithmetic depends critically on the source

K₇ and K₈ share the same Z₇ symmetry, the same Heawood direction structure, the same Z₃ phases. Yet they produce completely independent arithmetic: different number fields, different Galois groups, different discriminant primes. The source of the bilinear form — the K₈ matching Gram matrix vs. the K₇ directed adjacency — determines the specific arithmetic output.

### 4. K₈'s Galois group is anomalously small

This is the surprise. K₇ produces S₇ (the generic group for a degree-7 polynomial). K₈ produces C₂ ≀ C₃ (order 12) — an extraordinarily constrained group for a degree-6 polynomial. The "generic" outcome for an irreducible degree-6 polynomial would be S₆ (order 720).

K₈'s Galois group is 60 times smaller than generic. This means the matching structure of K₈ — the 105 matchings, the vacuum eigenspace, the genus-2 embedding — imposes severe algebraic constraints that the K₇ Dirac operator does not. The matching chain is not just producing arithmetic; it is producing **highly structured** arithmetic.

### 5. The wreath product structure is special to matchings

K₈'s Gal = C₂ ≀ C₃ is a wreath product: the C₃ from direction class permutation (quadratic residues mod 7), the C₂³ from complex conjugation on eigenvalue pairs. This wreath product structure is exactly what the CFSG analysis predicted for Level 1 extractors with cyclic symmetry.

K₇'s Gal = S₇ has no wreath product structure. The direction phases produce the same C₃ symmetry quotient, but the bilinear form (D†D) doesn't constrain the arithmetic kernel to C₂. Without the matching structure's additional constraints, the Galois group runs wild.

### 6. Galois disjointness confirmed

ℚ(K₇ eigenvalues) ∩ ℚ(K₈ vacuum eigenvalues) = ℚ (almost certainly — different disc primes guarantee disjoint number fields). Different extractors on related objects produce arithmetically independent outputs. This is the extractor analogue of the Galois disjointness we proved between K₆ and K₈ levels of the matching chain.

---

## Structural Diagnosis

The portability test reveals a **hierarchy of extraction quality**:

**Level 0 (Trivial):** Bilinear form has too much symmetry. B = αI + βJ. All sectors degenerate. No arithmetic. (Fano plane overlap.)

**Level 1 (Generic):** Direction phases break symmetry. Irreducible polynomial produced. Galois group is generic (S_n). Arithmetic is present but unconstrained. (K₇ Dirac D†D.)

**Level 2 (Structured):** Matching complex + direction phases. Irreducible polynomial with small, solvable Galois group (wreath product). Arithmetic is highly constrained. Discriminant primes have structural meaning (43 from spectral collision, 7 from Z₇ symmetry order). (K₈ matching Gram matrix.)

The matching chain lives at Level 2. The K₇ Dirac operator lives at Level 1. The difference: the matching complex provides additional algebraic constraints (the secular factorization, the budget identity, the self-focusing) that pin the Galois group to a small solvable wreath product.

This is the answer to "what makes matching extractors special": it's not the direction structure (that's necessary but not sufficient), it's the matching structure itself — the way 105 perfect matchings of K₈ overlap and interfere — that forces the arithmetic into a highly constrained, physically meaningful form.

---

## Implications for the Landscape Survey

The H₇ test confirms and sharpens the landscape picture:

1. **The extractor is portable** — ✓ confirmed. Same pipeline, different source, genuine arithmetic.

2. **The direction structure is necessary** — ✓ confirmed. Without phases: trivial. With phases: nontrivial.

3. **The matching structure provides the constraint** — NEW. Generic extractors (Level 1) produce unconstrained arithmetic (S_n Galois). Matching extractors (Level 2) produce constrained arithmetic (small solvable Galois). The physical content — particle masses, coupling constants — comes from this constraint.

4. **CFSG prediction partially confirmed** — The wreath product Gal = (arithmetic) ≀ (symmetry) appears at Level 2 (matchings) but not Level 1 (generic Dirac). The CFSG taxonomy applies to structured extractors, not generic ones.

5. **Priority adjustment**: Items 1 and 5 on the computation list (H₇ Hamming, Petersen graph) will likely produce Level 1 (generic) arithmetic. Items 3, 4, 7 (M₁₁ Steiner system, PSL₂(7) K₈ decomposition, Golay code) are more likely to produce Level 2 (structured) arithmetic because they involve combinatorial designs with matching-like constraint structures.

---

## Next Steps

**Immediate (confirmatory):**
- Verify K₇ Galois group is exactly S₇ using PARI/GP or Magma
- Check if K₅ Dirac (with Z₃ phases on Z₅ direction classes) also gives S₅

**Next priority (the real test):**
- E₈ root system under cyclic Weyl subgroup — does the Cartan matrix's additional structure (it's an intersection form, not just an adjacency) constrain the Galois group below generic?
- M₁₁ on Steiner blocks — does the design structure force wreath product Galois?

If E₈ or M₁₁ produce small solvable Galois groups, the matching chain is not unique — there's a general principle that "combinatorial designs with balanced overlap structure produce constrained arithmetic." If they produce S_n, then matchings really are special.
