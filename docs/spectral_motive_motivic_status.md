---
title: Motivic Status of the Spectral Motive: Terminal Analysis
section: Killed Approaches
status: active
---

# Motivic Status of the Spectral Motive: Terminal Analysis

## Summary

The spectral motive M(Kâ‚ˆ) does not fit any existing motivic framework. Four candidate identifications were pursued to termination, each blocked by a specific computable obstruction. The spectral motive is a genuinely new arithmetic object â€” more structured than an Artin motive, but not a geometric motive in Grothendieck's sense.

---

## 1. The Setup

The spectral motive M(Kâ‚ˆ) is built from:
- **E**: CM elliptic curve with j = 0, CM by â„¤[Ï‰] (Ï‰ = e^{2Ï€i/3})
- **Câ‚‡ âŠ‚ E[7]**: cyclic 7-torsion subgroup, kernel of the isogeny (3+Ï‰): E â†’ E'
- **Heawood map**: triangulation of Kâ‚‡ on the torus â„‚/((3+Ï‰)â„¤[Ï‰])
- **Gram matrix**: 105Ã—105 integer matrix on perfect matchings of Kâ‚ˆ, weighted by direction classes from Fâ‚‡*/Â±1 â‰… â„¤â‚ƒ
- **Vacuum sextic**: f(x) = xâ¶ âˆ’ 44xâµ + 720xâ´ âˆ’ 5648xÂ³ + 22512xÂ² âˆ’ 43456x + 31808, irreducible over â„š
- **Splitting field**: degree 6 over â„š, Gal = Câ‚‚ â‰€ Câ‚ƒ, contains â„š(âˆš43)
- **Discriminant**: Î”(f) = 2Â³â¶ Ã— 7â´ Ã— 43 Ã— 421Â²

The question: is there an algebraic variety X/â„š whose cohomology H^i(X) realizes Vâ‚â‚†,â‚‚â‚Ž (the 20-dimensional physical eigenspace) as a Galois module, with the Gram matrix as an endomorphism?

---

## 2. Approach 1: hÂ¹(Eâ¿) â€” CM Cohomology

**Hypothesis**: The Gram matrix is an endomorphism of EÂ³, acting on HÂ¹(EÂ³) â‰… â„š(Ï‰)Â³.

**Test**: If this were true, the sextic would factor as g(x)Â·á¸¡(x) over â„š(âˆšâˆ’3), where g âˆˆ â„š(âˆšâˆ’3)[x] has degree 3 and á¸¡ is its Galois conjugate.

**Obstruction**: The constant term satisfies g(0)Â·á¸¡(0) = |g(0)|Â² = eÂ² + 3fÂ² = 31808. But 31808 = 2â¶ Ã— 7 Ã— 71, and 71 â‰¡ 2 mod 3 appears to odd power. The form xÂ² + 3yÂ² represents n only if every prime p â‰¡ 2 mod 3 divides n to even power. Therefore 31808 is not a norm from â„š(âˆšâˆ’3).

**Stronger obstruction**: The splitting field contains â„š(âˆš43). The CM field is â„š(âˆšâˆ’3). These are disjoint: â„š(âˆš43) âˆ© â„š(âˆšâˆ’3) = â„š. No piece of h*(Eâ¿) can produce the number field â„š(âˆš43).

**Status**: TERMINATED. The prime 71 in f(0) and the disjointness of â„š(âˆš43) from â„š(âˆšâˆ’3) are sharp arithmetic obstructions.

---

## 3. Approach 2: Configuration Space Cohomology on E[7]

**Hypothesis**: Vâ‚â‚†,â‚‚â‚Ž = H^k(X) for some configuration space X built from E and E[7].

**Candidate**: HÂ³(SymÂ³(E)) appeared promising because bâ‚ƒ of a 3-dimensional abelian variety is C(6,3) = 20 = dim Vâ‚â‚†,â‚‚â‚Ž.

**Computation**: HÂ³(SymÂ³(E)) = HÂ³(EÂ³)^{Sâ‚ƒ}. Decomposing:
- (0,1,2)-sector: V^{âŠ•6} with Sâ‚ƒ acting â†’ invariants = V (dim 2)  
- (1,1,1)-sector: V^{âŠ—3} with Sâ‚ƒ acting â†’ invariants = SymÂ³(V) (dim 4)
- Total: dim HÂ³(SymÂ³(E)) = 6, not 20.

The bâ‚ƒ = 20 is for the abelian variety EÂ³, not the symmetric product SymÂ³(E).

**Remaining candidate**: HÂ³(EÂ³) is 20-dimensional but carries Sâ‚ƒ action (permuting factors), not Sâ‚ˆ action (from matchings). There is no natural Sâ‚ˆ action on EÂ³, and the Galois representation is through â„š(âˆšâˆ’3), not â„š(âˆš43).

**Status**: TERMINATED. Three independent obstructions: (a) wrong dimension for symmetric products, (b) wrong group action, (c) wrong arithmetic (â„š(âˆšâˆ’3) vs â„š(âˆš43)).

---

## 4. Approach 3: Feynman Motives

**Hypothesis**: The matching chain's Euler-product structure is a Feynman motive, with the Gram matrix encoding a graph hypersurface period.

**Observation**: Feynman motives arise from Kirchhoff polynomials Ïˆ_Î“, which encode spanning trees. The Gram matrix encodes matching overlaps. These are different graph invariants:
- Ïˆ_{Kâ‚‡}: degree-15 polynomial in 21 variables, encoding spanning trees
- Gram matrix: 105Ã—105 integer matrix, encoding matching overlaps

The matching polynomial of Kâ‚ˆ has roots {0.29, 2.68, 7.85, 17.18} (in t = xÂ²), which are unrelated to the Gram eigenvalues {1.96, 3.19, 4.16, 6.83, 10.06, 17.80}.

No known construction relates Kirchhoff polynomials to matching overlap eigenvalues.

**Status**: TERMINATED. The Feynman motive framework requires a graph hypersurface whose period is the quantity of interest. The matching chain's spectral data does not arise from any graph hypersurface integral.

---

## 5. Approach 4: Artin Motive (Trivial Identification)

**Fact**: The vacuum sextic f(x) defines a number field K = â„š(Î»_vac) of degree 6 with Gal(K/â„š) = Câ‚‚ â‰€ Câ‚ƒ. This automatically defines an Artin motive h(Spec K) with:

- **L-function**: Î¶_K(s)/Î¶(s), which decomposes into irreducible Artin L-functions
- **Functional equation**: standard, from class field theory
- **Automorphy**: guaranteed (Câ‚‚ â‰€ Câ‚ƒ is solvable â†’ Langlands-Tunnell)

This IS a genuine motive with a genuine L-function. But it is **trivially so**: any irreducible polynomial defines an Artin motive. The identification captures the number field but not the matching chain's internal structure.

**What the Artin motive does NOT see**:
- The secular factorization char_poly(G_phys) = âˆâ‚– Sâ‚–(Î»)
- The budget identity Î»_mid Â· d_phys = 2(nâˆ’1) Â· Î»_max
- The self-focusing convergence Î·(n) â†’ 1
- The Ï†(2nâˆ’1) degree law across levels
- The Galois disjointness between levels

The Artin motive sees the output polynomial, not the constructive process. It's like knowing that 12 = 2Â² Ã— 3 defines a number field â„š(12^{1/n}) without knowing the prime factorization of 12.

**Status**: TERMINATED (as a non-trivial identification). The Artin motive exists but does not capture the spectral motive's distinguishing structure.

---

## 6. Terminal Diagnosis

### What the spectral motive IS NOT:
1. **Not a Grothendieck motive**: No algebraic variety X/â„š with H^i(X) â‰… Vâ‚â‚†,â‚‚â‚Ž
2. **Not a Feynman motive**: No graph hypersurface connection
3. **Not a CM motive**: Arithmetic is independent of the CM curve (â„š(âˆš43) âŠ¥ â„š(âˆšâˆ’3))
4. **Not a configuration space cohomology**: Dimension, group action, and arithmetic all mismatch

### What the spectral motive IS:
5. **An Artin motive** (trivially â€” any number field gives one, structure not captured)
6. **An integer matrix spectrum** with internal factorization structure (secular factors)
7. **Built FROM a CM elliptic curve** (Heawood map = dessin d'enfant on E with j = 0) **but producing independent arithmetic** (new primes 43, 71, 421; new field â„š(âˆš43))

### The creation mechanism:
```
CM curve E (j=0, Z[Ï‰]) 
  â†’ Heawood map (E[7] with Zâ‚‡) 
    â†’ 105 matchings of Kâ‚ˆ  [COMBINATORIAL STEP]
      â†’ Gram matrix (integer)  [LINEAR ALGEBRA STEP]
        â†’ vacuum sextic  [SPECTRAL STEP]
          â†’ Q(âˆš43)  [ARITHMETIC STEP]
```

The new arithmetic enters at the **combinatorial step**: enumerating perfect matchings on the 7-torsion. The matching complex M(E[7]) is not a subvariety of Eâ¿. It is a combinatorial object associated TO the torsion but not derivable FROM the cohomology.

### The gap:
There is no known mathematical framework that:
- Takes as input: an integer matrix from graph combinatorics on torsion points
- Produces as output: a motive whose L-function has Euler product matching the secular factorization
- Such that the budget identity becomes the functional equation

This gap identifies a **genuinely new structure**: more than an Artin motive (has internal factorization), less than a geometric motive (no underlying variety), parallel to Feynman motives (graph-defined, producing new arithmetic from combinatorial data) but distinct (different graph polynomial, different integral).

---

## 7. What This Means for the Program

### The spectral motive sits in a precise position:

**Below** the motivic landscape: it does not have Ã©tale/Betti/de Rham realizations that satisfy the standard comparison isomorphisms. Its five realizations (spectral, arithmetic, algebraic, evaluative, combinatorial) are a *new* set of "cohomology theories" that converge on the same invariants through a different mechanism.

**Parallel to** the motivic landscape: it exhibits the same wall structure (polynomial side / series side, Galois disjointness, Euler-type factorization, functional-equation-type constraint). The Ï†(2nâˆ’1) degree law gives the same cyclotomic mechanism that controls classical L-functions.

**Fed by** the motivic landscape: it is built from a genuine CM motive (the elliptic curve E with j = 0), using its torsion structure as scaffolding. The Heawood map IS a dessin d'enfant, which IS a motivic object. The spectral motive is a secondary construction on this motivic input.

### For RH:
The spectral motive's L-function (the Artin L-function of â„š(Î»_vac)) is a genuine L-function to which GRH applies. But GRH for this specific L-function is a statement about the *number field* â„š(Î»_vac), not about the matching chain's internal structure. The budget identity constrains eigenvalues at each finite level; GRH for the Artin L-function constrains something else (zeros of Î¶_K(s)/Î¶(s) on Re(s) = 1/2). These are different constraints on related but distinct objects.

The honest position: the matching chain exhibits the wall, produces pure numbers in the arithmetic universe, and generates genuine L-functions to which RH applies. But it has not crossed into the motivic landscape in a way that would make its internal structure (secular factorization, budget identity) equivalent to standard motivic structure (Euler product, functional equation).

### The open question, precisely stated:
Is there a categorical framework â€” extending Grothendieck's motives to include "combinatorial motives" built from matching complexes on torsion points â€” in which the spectral motive's secular factorization IS an Euler product, the budget identity IS a functional equation, and the Ï†(2nâˆ’1) degree law IS the degree of a cyclotomic Frobenius? If such a framework exists, the spectral motive enters the motivic landscape and the wall becomes the same wall. If not, the spectral motive remains a parallel construction that illuminates the wall from outside.
