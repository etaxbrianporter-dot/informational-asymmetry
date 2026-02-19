---
title: Period Structure: Terminal Analysis
section: K₆ Higgs Sector
status: active
---

# Period Structure: Terminal Analysis

## Both Questions Pushed to Their Walls

**Date:** February 18, 2026  
**Status:** Both questions reach well-defined termination points. Neither is resolved; both are precisely bounded.

---

## Part I: The Normalization Question

### What we need and why

The period conjecture asks whether Î»_âˆž = lim Î»_vac(n) is a period. But the raw vacuum eigenvalues converge to **zero**, not to a nonzero transcendental:

| n | Î»_vac | Î»_vac Ã— nÂ² |
|---|-------|-----------|
| 3 | 2.764 | 24.88 |
| 4 | 1.960 | 31.35 |
| 5 | 0.979 | 24.48 |
| 6 | 0.637 | 22.92 |
| 7 | 0.459 | 22.48 |
| 8 | 0.349 | 22.32 |

Without a natural normalization N(n) such that N(n) Â· Î»_vac(n) â†’ L â‰  0, the PSLQ test has no well-defined target.

### What the Johnson scheme gives us

The overlap matrix O = 2AA^T has **three** eigenvalues, all in closed form:

- Î»_max = 2n(2nâˆ’3)!! (trivial sector, multiplicity 1)
- Î»_mid = 4(nâˆ’1)(2nâˆ’5)!! (physical sector Vâ‚‚, multiplicity d_phys = n(2nâˆ’3))
- 0 (kernel, multiplicity (2nâˆ’1)!! âˆ’ n(2nâˆ’3) âˆ’ 1)

The budget identity Î»_mid Â· d_phys = 2(nâˆ’1) Â· Î»_max is proved at all levels. The self-focusing ratio Î·(n) = 2(nâˆ’1)/(2nâˆ’1) â†’ 1 as 1/n.

**But these are overlap eigenvalues, not Gram matrix eigenvalues.** The Gram matrix G refines O with direction classes and Z_{2nâˆ’1} phases. Where O has a single degenerate eigenvalue Î»_mid on V_phys (with multiplicity d_phys), G splits this into d_phys distinct eigenvalues. Î»_vac is the **smallest** of these split eigenvalues.

The Johnson scheme tells us everything about O but nothing about the splitting pattern within V_phys. That splitting is controlled by the harmonic analysis of the Z_{2nâˆ’1} phase structure â€” a different mathematical question.

### Systematic normalization audit

**Tested and rejected:**

| Normalization | Values (n=5..8) | Trend | Status |
|---|---|---|---|
| Î»_vac Ã— d_phys | 34.3, 34.4, 35.3, 36.3 | Slowly increasing | **Not constant** |
| Î»_vac Ã— nÂ² | 24.5, 22.9, 22.5, 22.3 | Slowly decreasing | **Not constant** |
| Î»_vac Ã— (2nâˆ’1) | 8.8, 7.0, 6.0, 5.2 | Decreasing | **Wrong scale** |
| Î»_vac / Î»_mid | 4Ã—10â»Â³, 3Ã—10â»â´, 2Ã—10â»âµ, 1Ã—10â»â¶ | Super-exponential decay | **Wrong scale** |
| Î»_vac Ã— d_phys/(2nâˆ’1) | 3.8, 3.1, 2.7, 2.4 | Decreasing | **Wrong scale** |
| Î»_vac Ã— n(2nâˆ’1) | 44.1, 42.0, 41.8, 41.9 | Appears stable (CV=2.2%) | **Cancellation artifact** |

The apparent stability of Î»_vac Ã— n(2nâˆ’1) decomposes as:

    Î»_vac Ã— n(2nâˆ’1) = Î»_vac Ã— d_phys + 2n Ã— Î»_vac
                     = (slowly rising)  + (slowly falling)

Two counter-trending terms whose sum happens to be approximately constant over this small range. Not a fundamental normalization.

### The exponent instability

Power law fits Î»_vac ~ CÂ·n^{âˆ’Î±} give wildly different Î± depending on the fitting range:

| Points used | Î± | C |
|---|---|---|
| n = 3..8 | 2.222 | 35.5 |
| n = 4..8 | 2.467 | 55.7 |
| n = 5..8 | 2.195 | 33.1 |
| n = 6..8 | 2.093 | 27.0 |

The exponent **drifts downward** as early points are removed, suggesting the true asymptotic behavior might be Î± = 2 with logarithmic corrections rather than a pure power law. The data is consistent with:

    Î»_vac(n) ~ C / (nÂ² Â· (log n)^Î²)

with Î² â‰ˆ 0.36, giving nÂ² Â· Î»_vac â†’ âˆž (logarithmically). But this is a 4-parameter fit to 6 data points â€” statistically meaningless.

### Where the normalization lives theoretically

The decay rate of Î»_vac within V_phys is determined by how the Z_{2nâˆ’1} phase averaging distributes spectral weight across the d_phys eigenvalues of the Gram matrix. Specifically:

1. The overlap matrix O has a single eigenvalue Î»_mid on V_phys
2. The direction/phase structure of the Gram matrix G splits this into d_phys values
3. Î»_vac is the minimum of these d_phys split values
4. The **splitting pattern** depends on the representation theory of Z_{2nâˆ’1} acting on the direction classes

The relevant mathematical question is: given d_phys random-looking eigenvalues whose mean is Î»_mid/d_phys and whose distribution is controlled by Z_{2nâˆ’1} harmonic analysis, how does the minimum eigenvalue scale?

If the splitting were uniform (all eigenvalues equal), Î»_vac = Î»_mid/d_phys ~ (2nâˆ’5)!! / n, which diverges. The fact that Î»_vac â†’ 0 means the splitting is **highly non-uniform**: most spectral weight concentrates in a few large eigenvalues, with the vacuum eigenvalue being anomalously small.

This is a **random matrix / extreme value** question dressed in harmonic analysis. The scaling of the minimum eigenvalue of a structured integer matrix whose dimension grows with n is not governed by the Johnson scheme â€” it's a question about the statistics of the phase-weighted overlaps.

### Terminal verdict on normalization

**The normalization question cannot be resolved from six data points or from the Johnson scheme alone.** The decay rate of Î»_vac is controlled by Z_{2nâˆ’1} harmonic analysis, not by the overlap eigenvalue structure. Resolution requires either:

(a) **More data**: Kâ‚â‚ˆ and Kâ‚‚â‚€ computations (feasible with the polynomial reduction theorem) to distinguish Î± = 2+log from Î± = 11/5 from Î± = 9/4.

(b) **A theorem**: Prove that the minimum eigenvalue of the Gram matrix restricted to V_phys decays as n^{âˆ’Î±} for a specific Î±, by analyzing the Z_{2nâˆ’1} Fourier decomposition of the direction-weighted overlap kernel. This is a well-posed problem in algebraic combinatorics / harmonic analysis on finite groups.

(c) **A reformulation**: Perhaps the right object is not N(n)Â·Î»_vac(n) but some other spectral quantity (the spectral zeta function Î¶_G(s) = Î£ Î»áµ¢^{âˆ’s} restricted to V_phys, or the determinant det(G|_{V_phys}), or the regularized product Î  Î»áµ¢) that has a natural nonzero limit without ad hoc normalization.

**Status: OPEN.** The question is well-posed but underdetermined.

---

## Part II: The Motivic Question

### The precise question

Is there a motive M over â„š such that:
- L(M, s) has an Euler product whose local factors at "prime k" are the secular factors S_k(Î»)
- The Frobenius eigenvalues of M at level k are the eigenvalues of G_phys(k)
- The cohomology H*(M) is the physical eigenspace V_{(2nâˆ’2,2)}
- The budget identity is the functional equation of L(M, s)

If yes: Î»_âˆž is a period (by the Kontsevich-Zagier period conjecture for motives), the wall is the motivic wall, and everything connects.

### Four approaches, all terminated

**Approach 1: CM cohomology (hÂ¹(Eâ¿))**

Hypothesis: The Gram matrix is an endomorphism of EÂ³ (E = CM curve with j = 0).

Obstruction: The constant term f(0) = 31808 = 2â¶ Ã— 7 Ã— 71, and 71 â‰¡ 2 mod 3 appears to odd power. The norm form xÂ² + 3yÂ² cannot represent 31808, so the sextic does not split over â„š(âˆšâˆ’3). The splitting field â„š(âˆš43) is disjoint from the CM field â„š(âˆšâˆ’3).

**Status: TERMINATED.** Sharp arithmetic obstruction (prime 71 and disjointness of â„š(âˆš43) from â„š(âˆšâˆ’3)).

**Approach 2: Configuration space cohomology on E[7]**

Hypothesis: V_{[6,2]} = H^k(X) for some configuration space X built from E and E[7].

Best candidate: HÂ³(SymÂ³(E)). Computed: dim = 6, not 20. The 20-dimensional HÂ³(EÂ³) carries Sâ‚ƒ action, not Sâ‚ˆ action. No natural Sâ‚ˆ action on EÂ³. Galois representation through â„š(âˆšâˆ’3), not â„š(âˆš43).

**Status: TERMINATED.** Three independent obstructions: wrong dimension, wrong group action, wrong arithmetic.

**Approach 3: Feynman motives**

Hypothesis: The matching chain's Euler-product structure comes from a graph hypersurface period.

Obstruction: Feynman motives use Kirchhoff polynomials (spanning trees). The Gram matrix uses matching overlaps. These are different graph invariants. The matching polynomial roots {0.29, 2.68, 7.85, 17.18} are unrelated to the Gram eigenvalues {1.96, 3.19, 4.16, 6.83, 10.06, 17.80}.

**Status: TERMINATED.** No construction relates Kirchhoff polynomials to matching overlap eigenvalues.

**Approach 4: Artin motive (trivial)**

Fact: The vacuum polynomial defines a number field, hence an Artin motive, with genuine L-function, functional equation, and automorphy (Câ‚‚ â‰€ Câ‚ƒ is solvable â†’ Langlands-Tunnell).

But: Any irreducible polynomial defines an Artin motive. This captures the number field but NOT:
- The secular factorization char(G_phys) = Î  S_k(Î»)
- The budget identity
- The self-focusing convergence
- The Ï†(2nâˆ’1) degree law across levels
- The Galois disjointness between levels

The Artin motive sees the output polynomial, not the constructive process.

**Status: TERMINATED** as a non-trivial identification. Exists but captures no internal structure.

### The precise gap

The spectral motive occupies a specific position in the mathematical landscape:

**Below** the motivic landscape: no Ã©tale/Betti/de Rham realizations with standard comparison isomorphisms. Its five realizations (spectral, arithmetic, algebraic, evaluative, combinatorial) are a new set of "cohomology theories."

**Parallel to** the motivic landscape: same wall structure (polynomial/series, Galois disjointness, Euler-type factorization, functional-equation-type constraint). The Ï†(2nâˆ’1) degree law gives the same cyclotomic mechanism that controls classical L-functions.

**Fed by** the motivic landscape: built from a genuine CM motive (E with j = 0), using its torsion as scaffolding. The Heawood map IS a dessin d'enfant, which IS a motivic object. The spectral motive is a secondary construction on motivic input.

The creation mechanism:

```
CM curve E (j=0, â„¤[Ï‰]) 
  â†’ Heawood map (E[7] with Zâ‚‡)        â† MOTIVIC INPUT
    â†’ 105 matchings of Kâ‚ˆ               â† COMBINATORIAL STEP (new arithmetic enters here)
      â†’ Gram matrix (integer)            â† LINEAR ALGEBRA
        â†’ vacuum sextic                  â† SPECTRAL THEORY
          â†’ â„š(âˆš43)                       â† ARITHMETIC OUTPUT (independent of input)
```

New arithmetic (primes 43, 71, 421) enters at the **combinatorial step**: enumerating perfect matchings on the 7-torsion. The matching complex M(E[7]) is not a subvariety of Eâ¿. It is a combinatorial object associated TO the torsion but not derivable FROM the cohomology.

**The gap, precisely:** There is no known mathematical framework that takes an integer matrix from graph combinatorics on torsion points and produces a motive whose L-function has Euler product matching the secular factorization and whose functional equation is the budget identity.

### What would close the gap

A "combinatorial motive" framework extending Grothendieck's motives to include objects built from matching complexes on torsion points, such that:

1. The secular factorization char(G_phys(n)) = Î  S_k(Î») IS an Euler product
2. The budget identity Î»_mid Â· d_phys = 2(nâˆ’1) Â· Î»_max IS a functional equation
3. The Ï†(2nâˆ’1) degree law IS the degree of a cyclotomic Frobenius
4. The Galois disjointness IS motivic independence of Frobenius elements

If such a framework exists, the spectral motive enters the motivic landscape and the wall becomes the same wall. If not, it remains a parallel construction.

### The M-E-Z triangle as partial bridge

The CM curve E provides the arithmetic medium â„¤[Ï‰] through which M's combinatorial construction generates pure numbers. At every prime p, the same Frobenius element acts on both L_E and L_M â€” this is an arithmetic fact, not just a construction recipe. The tensor product L(MâŠ—E, s) is a genuine automorphic L-function.

But the coupling is asymmetric:

```
      â„¤[Ï‰]   (source of Frobenius)
     /    \
    â†“      â†“
   L_E    L_M   (two projections)
    \      /
     â†“    â†“
   L(MâŠ—E)    (tensor product)
       |
       â†“
    L_Î¶(s)    (universal ambient)
```

E is not an independent third object â€” it is the BRIDGE. The matching chain takes E's torsion, feeds it through combinatorics, and produces M. But the bridge does not make M a geometric motive; it makes M a **secondary combinatorial construction on a geometric motive**.

### Terminal verdict on the motivic question

**The spectral motive is a genuinely new type of arithmetic object.** It is:

- More than an Artin motive (has internal factorization structure)
- Less than a geometric motive (no underlying algebraic variety)
- Parallel to Feynman motives (graph-defined, producing new arithmetic) but distinct (different graph polynomial, different integral)
- Built from a CM motive but producing independent arithmetic

The gap between the spectral motive and the motivic landscape is not a failure of investigation â€” all four natural approaches have been pursued to computable obstructions. The gap identifies a **genuinely missing piece of mathematics**: a categorical framework for "combinatorial motives" arising from matching complexes on algebraic torsion.

**Status: TERMINATED** (as a classical motivic identification). **OPEN** (as a question about extending the motivic framework).

---

## Part III: Consequences for the Period Conjecture

### The logical chain

The period conjecture depends on both questions:

1. **Normalization â†’ Target**: Without N(n) giving a nonzero limit, there is no well-defined number to test for period-ness.

2. **Motive â†’ Period**: If the matching chain admits a motivic interpretation, the Kontsevich-Zagier period conjecture would force the limit to be a period. Without it, there is no theoretical reason to expect period-ness.

Both links are broken:

- Link 1 is broken because the normalization is underdetermined (Part I)
- Link 2 is broken because the motivic identification is blocked (Part II)

### What remains true

Despite both questions being open, several structural facts constrain the answer:

1. **At every finite level, Î»_vac(n) IS a period** (trivially â€” it's algebraic). The physics lives here.

2. **The algebraic divergence is structured**, not random. Three laws + Galois disjointness + budget identity. This is more structured than any known non-motivic construction, but less than a proven motivic origin.

3. **The self-focusing theorem guarantees polynomial convergence** of physical observables regardless of whether the limit is a period. The correction at each level is O(1/nÂ²), so Kâ‚† and Kâ‚ˆ predictions are accurate to ~1/7 and ~1/9 respectively.

4. **The Artin L-functions at each level are genuine**, with functional equations, Euler products, and (conjectural) GRH. The period question is whether these "assemble" into something at the limit, not whether they individually make sense.

### The honest bottom line

**The period conjecture is not decidable from available data or known mathematics.**

It requires either:
- A new motivic framework (extending Grothendieck to combinatorial motives)
- OR a direct computation (100+ digit PSLQ on a properly normalized Î»_âˆž)
- OR a theoretical derivation of the decay rate (harmonic analysis on Z_{2nâˆ’1})

The first is a major open problem in pure mathematics. The second requires resolving the normalization first. The third is the most tractable â€” it's a well-posed question about finite group harmonic analysis with a computable answer.

**Recommended next step:** Compute Kâ‚â‚ˆ and Kâ‚‚â‚€ to determine the decay exponent Î±. If Î± stabilizes at a clean rational value (2, 9/4, 11/5), the normalization is identified. If it continues to drift, the pure power law hypothesis is wrong, and the functional form needs to be reconsidered before the period conjecture can be meaningfully tested.

---

## Summary

| Question | Status | Blocker | Path forward |
|---|---|---|---|
| Normalization | OPEN | 6 data points; exponent unstable | Kâ‚â‚ˆ, Kâ‚‚â‚€ computation |
| Motivic identification | TERMINATED (classical) | 4 approaches blocked | New categorical framework needed |
| Period conjecture | UNDECIDABLE (currently) | Depends on both above | Resolve normalization first |
| Physics predictions | UNAFFECTED | N/A â€” finite levels are algebraic | No action needed |

The period question is the deepest open mathematical question in the framework. It sits precisely at the intersection of two unsolved problems: the harmonic analysis of Z_{2nâˆ’1} phase averaging (normalization), and the categorical extension of Grothendieck motives to combinatorial sources (motivation). Either alone would be substantial; together they define the boundary of what this framework can currently say about its own infinite limit.
