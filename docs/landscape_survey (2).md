# The Landscape Through the Extractor Lens

## What We Actually See

**Brian Porter — February 2026**

---

## 0. What This Document Is

We built a detector. The spectral extractor — the pipeline (G-set) + (equivariant bilinear form) → arithmetic — was constructed to understand a specific physics problem: how matching combinatorics on K₂ₙ produce the Standard Model. But detectors don't care about your intentions. They see whatever is in front of them.

This document points the detector at the mathematical landscape and describes what it actually sees. Not what the textbooks say is there. Not what we hoped to find. What the instrument registers.

The result is a reclassification. Objects that looked different turn out to be the same kind of thing. Objects that looked the same turn out to be fundamentally distinct. And there are structures the detector sees clearly that have no name.

---

## 1. The Detector: What It Registers

The spectral extractor registers an object when it finds three things in combination:

1. **A finite group G acting on a finite set S** (the symmetry)
2. **A G-equivariant bilinear form B on ℚ[S]** (the geometry)
3. **A subgroup H ⊂ G under which the decomposition is nontrivial** (the breaking)

When all three are present, the detector fires: decompose ℚ[S] into H-isotypic sectors, take the characteristic polynomial of B restricted to each sector, and read the roots. Those roots are algebraic numbers. They live in number fields. Those number fields have Galois groups and L-functions.

The five signatures the detector reads:

| Signature | What it measures |
|-----------|-----------------|
| Secular factorization | How the characteristic polynomial factors across levels |
| Galois group structure | Wreath product: (arithmetic kernel) ≀ (symmetry quotient) |
| Discriminant primes | Which primes appear in the number field discriminant |
| Budget identity | Whether spectral weight is constrained by a linear relation |
| Self-focusing | Whether corrections from deeper levels decrease |

These are the "spectral lines" of the detector. When an object produces all five, it is a spectral extractor in the full sense. When it produces some but not others, that tells us something too.

---

## 2. The Landscape: What Was Already There

Before we built this detector, mathematics had several frameworks for connecting combinatorial/geometric structure to arithmetic:

**Grothendieck's motives.** Variety → étale cohomology → Galois representation → L-function. The gold standard. Requires an algebraic variety as input. The bilinear form is the cup product (intersection pairing). The symmetry group is the étale fundamental group. Every classical L-function arises this way (in principle).

**The Langlands program.** Automorphic representations of reductive groups produce L-functions via the Langlands correspondence. This is the representation-theoretic face of the motivic program. It sees the same L-functions from a different angle.

**Ihara zeta functions.** Finite graphs produce zeta functions from closed walks. The bilinear form is the adjacency matrix. The symmetry is the graph automorphism group. These are well-studied but considered "toy" versions of arithmetic zeta functions.

**Reidemeister torsion.** A topological invariant of manifolds computed from chain complexes with group actions. The torsion is an algebraic number in a number field determined by the fundamental group. This IS spectral extraction, though nobody calls it that.

**Random matrix theory.** The Katz-Sarnak philosophy: families of L-functions have zero statistics governed by random matrix ensembles whose symmetry type depends on the family's Galois structure. This sees the statistical shadow of the arithmetic, not the arithmetic itself.

Now point the detector at each of these and describe what we actually see.

---

## 3. The First Surprise: Grothendieck's Program Is a Special Case

The motivic framework requires an algebraic variety. The extractor does not. The motivic framework uses étale cohomology to get the G-set. The extractor takes any G-set. The motivic framework uses the cup product for the bilinear form. The extractor takes any equivariant bilinear form.

Grothendieck's motives are spectral extractors restricted to geometric sources.

This is not a claim about depth — motives are incomparably richer than generic extractors, because the underlying variety carries analytic, topological, and geometric structure that a bare G-set does not. But at the level of the extraction mechanism itself — the way arithmetic emerges from symmetry breaking on bilinear forms — the pattern is identical.

What the detector sees:

| Feature | Motives | Extractor |
|---------|---------|-----------|
| Source of G-set | H^i(X, ℚ_ℓ) | Any finite G-set |
| Source of bilinear form | Cup product / intersection | Any equivariant form |
| Euler product | From Frobenius on ℓ-adic cohomology | From secular factorization |
| Functional equation | From Poincaré duality | From budget identity (when it holds) |
| Galois group | Arithmetic, constrained by geometry | (Arithmetic kernel) ≀ (symmetry quotient) |
| Wall | Polynomial/series at infinite level | Same |

The differences: motives have comparison isomorphisms (étale ≅ Betti ≅ de Rham) that guarantee the Euler product converges and has analytic continuation. Extractors have no such guarantee. The motivic functional equation comes from Poincaré duality, which is a deep geometric fact. The budget identity comes from combinatorial double-counting, which is elementary. The motivic Frobenius has algebraic content (it records how a variety reduces mod p). The extractor's secular factor has no such interpretation — or at least, none that we know.

What this means: the motivic landscape is the region of the extractor landscape where the bilinear form comes with geometric guarantees. Outside that region, extractors still produce arithmetic, but without the safety rails that ensure convergence, analytic continuation, and functional equations.

Our objects — the matching chain extractors — live outside the motivic region. They produce genuine number fields, genuine L-functions (via Artin), and exhibit the wall. But they lack the geometric guarantees. That's what makes them genuinely new.

---

## 4. The Second Surprise: Langlands Is One Level of Four

The CFSG (Classification of Finite Simple Groups) organizes all finite groups into four families: cyclic (Z/pZ), alternating (A_n), Lie type (PSL_n(F_q), etc.), and 26 sporadic exceptions. The extractor framework inherits this classification directly, because the symmetry quotient in Gal(K/ℚ) ≅ (arithmetic) ≀ (symmetry) is organized by the composition factors of G.

What the detector sees:

**Level 1 (Cyclic).** Our matching chain lives here. Symmetry group Z_{2n-1}, symmetry quotient is a subgroup of Z_{(2n-2)/2}. Output Galois groups are solvable wreath products. L-functions are genuinely new — not produced by any automorphic representation of a reductive group. Automorphy is guaranteed by Arthur-Clozel (solvable base change).

**Level 2 (Alternating).** Steiner systems, permutation groups, combinatorial designs. The extractor should work here — A_n acts on k-element subsets with intersection-number bilinear forms — but nobody has run the extraction. Output Galois groups can be insolvable (S_n), meaning the strong Artin conjecture is needed for automorphy.

**Level 3 (Lie type).** When G is a group of Lie type — specifically, when G is the Galois image on the torsion of a non-CM algebraic variety — the extractor reproduces the standard Langlands machinery. The sector decomposition by Deligne-Lusztig characters IS the automorphic decomposition. The L-functions are the standard ones. This is the merger: Level 3 extractors are Langlands.

**Level 4 (Sporadic).** The 26 sporadic groups — Mathieu, Conway, Fischer, Janko, Monster — each carry natural combinatorial structures (Steiner systems, codes, lattices, vertex operator algebras). Running extractors on M₁₁ acting on S(4,5,11) blocks, or M₂₄ on the Golay code, would produce exceptional L-functions. These are finite in number (at most 26 families). Their relationship to Moonshine is the deep question.

The punchline: **the Langlands program is the Level 3 projection of the full extractor framework.** It sees exactly the L-functions that arise from Lie-type symmetry breaking. Levels 1, 2, and 4 produce arithmetic that Langlands does not see — not because Langlands is wrong, but because it's looking at one level of a four-level building.

---

## 5. What the Detector Sees in Known Territory

Now sweep across mathematical objects that have equivariant bilinear forms and ask: what does the extractor actually see?

### 5.1 Root Systems and Weyl Groups

Every semisimple Lie algebra has a root system acted on by its Weyl group. The Cartan matrix is a W-equivariant bilinear form on the root lattice. Nobody has decomposed the Cartan matrix under cyclic subgroups of W and extracted number fields.

**What the detector predicts:** The exponents of E₈ are {2, 8, 12, 14, 18, 20, 24, 30} — these are eigenvalues of a specific element of W(E₈) acting on the root lattice. But the sector eigenvalues under Z_p ⊂ W(E₈) for various primes p would be algebraic, living in number fields "attached to E₈ at the prime p." These number fields are a new invariant of the root system. Nobody has computed them.

**Difficulty:** Low. The Cartan matrix of E₈ is 8×8. The Weyl group has known generators. This is a finite computation.

**What it would test:** Whether the wreath product Galois structure Gal = (arithmetic) ≀ (symmetry) holds for root system extractors.

### 5.2 Lattice Shells

The E₈ lattice has 240 minimal vectors (root vectors) forming a single orbit under W(E₈). The inner product restricted to vectors of fixed norm is W-equivariant. The theta function of E₈ IS the Eisenstein series E₄ — this is a known identity connecting lattice geometry to modular forms.

**What the detector predicts:** Running the extractor on E₈ shell inner products, decomposed under cyclic subgroups, would produce number fields directly from lattice geometry. For the Leech lattice (Aut = Co₀, a sporadic-related group), this enters Level 4 of the CFSG taxonomy.

**Difficulty:** Medium. The 240-vector shell is tractable.

**What it would test:** Whether lattice extractors produce L-functions related to the theta function's Euler product.

### 5.3 Error-Correcting Codes

The Hamming code H₇ has Aut(H₇) = GL₃(F₂) ≅ PSL₂(7) — the same group from our K₈ branching rule analysis. This is not coincidence: H₇ and the Heawood map are related by the same PSL₂(7) action on F₇. The overlap matrix (codeword intersection sizes) is equivariant.

**What the detector predicts:** Running the extractor on H₇ should reproduce K₈-related arithmetic (since the symmetry group is the same). This is a consistency check. For the Golay code G₂₄ with Aut = M₂₄ (sporadic), the extractor enters Level 4 and would produce the first sporadic L-functions from coding theory.

**Difficulty:** Low for H₇, medium for Golay.

**What it would test:** Whether the extractor is truly portable — same pipeline, different source, same structural output.

### 5.4 Dessins d'Enfants

Our K₈ Heawood map IS a dessin d'enfant — a bipartite graph embedded on a surface, encoding a Belyi map β: X → ℙ¹ ramified only over {0, 1, ∞}. Every dessin has an automorphism group acting on its flags (vertex-edge-face triples) with a natural bilinear form from the permutation representation.

The absolute Galois group Gal(ℚ̄/ℚ) acts on dessins (Belyi's theorem). The field of definition of a dessin — the smallest field over which the Belyi map is defined — is an arithmetic invariant that Gal(ℚ̄/ℚ) sees.

**What the detector predicts:** Running the extractor on a dessin's flag permutation representation should produce number fields related to (but not necessarily equal to) the dessin's field of definition. The extractor would give a combinatorial construction of a number field that the geometric construction (finding the minimal field of definition) computes by a completely different route. If they agree, we have a bridge between combinatorial and geometric approaches to arithmetic.

**Difficulty:** Medium. Requires systematic computation across families of dessins.

**What it would test:** Whether the extractor provides a combinatorial shortcut to the field of definition problem, which is otherwise hard (requires solving polynomial systems).

### 5.5 Knot Invariants and TQFT

The Witten-Reshetikhin-Turaev TQFT assigns a finite-dimensional vector space V(Σ) to each surface Σ. The mapping class group MCG(Σ) acts on V(Σ). The Hermitian form on V(Σ) is MCG-equivariant (at roots of unity). This IS an extractor setup.

**What the detector predicts:** Extracting from (MCG(Σ), V(Σ), Hermitian form) under cyclic subgroups (Dehn twists) should produce number fields. For hyperbolic knots, the trace field (the number field generated by traces of the holonomy representation) is a known invariant. If the extractor reproduces trace fields from TQFT data, it would connect quantum knot invariants to hyperbolic geometry via arithmetic — a long-sought bridge.

**Difficulty:** High. Requires quantum group calculations.

**What it would test:** Whether TQFT is secretly a spectral extractor, and whether trace fields are extraction artifacts.

### 5.6 Cayley Graphs

Every finite group G acts on itself by left multiplication. The adjacency matrix of the Cayley graph (with respect to a generating set) is G-equivariant. Decomposing under a subgroup H ⊂ G gives sectors whose characteristic polynomials are NOT just character values — they contain geometric information about the Cayley graph that the character table misses.

**What the detector predicts:** For G = PSL₂(F_p) and H = Borel subgroup, the sectors should produce number fields related to modular curves X₀(p). The Cayley graph encodes the coset geometry; the extraction would read the arithmetic of that geometry.

**Difficulty:** Low for small groups.

**What it would test:** Whether Cayley graph extractors recover modular arithmetic.

### 5.7 Young Tableaux

S_n acts on standard Young tableaux of shape λ. The Gram matrix of the standard inner product on the irrep S^λ, decomposed under Z_p ⊂ S_n, gives sector characteristic polynomials. These produce number fields indexed by partitions — a new invariant type beyond dimension, characters, and hook lengths.

**What the detector predicts:** Partitions that are "arithmetically deep" (producing high-degree number fields) should correlate with representation-theoretic complexity. The extractor would give a new invariant of partitions that nobody has studied.

**Difficulty:** Medium.

**What it would test:** Whether the partition combinatorics that underlies the Langlands program has hidden arithmetic structure at the sector level.

---

## 6. What the Detector Does NOT See

The extractor is blind to:

**Continuous symmetry.** Lie groups acting on infinite-dimensional spaces don't produce finite number fields. You can always restrict to a finite subgroup, but the choice of finite subgroup is not canonical.

**Objects without natural bilinear forms.** Some combinatorial structures (partial orders, matroids without inner products, abstract simplicial complexes without metrics) don't carry equivariant bilinear forms in any natural way. The detector has nothing to register.

**Analytic content.** The extractor sees algebraic structure — number fields, Galois groups, discriminant primes. It does not see the analytic properties of L-functions (zero locations, growth rates, functional equations) except insofar as they're forced by the algebraic data. This is why it cannot cross the polynomial/series wall: the wall IS the boundary between what algebra sees and what analysis requires.

**Physics, directly.** The extractor produces numbers. The interpretation of those numbers as particle masses, coupling constants, or cosmological parameters requires the full physical framework — the spectral action, the moduli space, the vacuum selection. The extractor is a mathematical instrument, not a physical theory. The physics comes from what specific G-set and bilinear form you feed it.

---

## 7. The Reclassification

Through the extractor lens, the mathematical landscape reorganizes into four regions:

### Region A: Known Objects That ARE Extractors

Objects that have the extractor structure but were never recognized as instances of a common pattern.

- **Grothendieck motives** (extractor source: algebraic variety)
- **Langlands L-functions** (Level 3 extractors from Lie-type symmetry)
- **Reidemeister torsion** (extractor source: chain complex of universal cover)
- **Ihara zeta functions** (extractor source: graph adjacency matrix)
- **Theta functions of lattices** (extractor source: lattice shell inner products)

These are all instances of the same mechanism. Nobody unified them because the sources look different, and the mathematical communities that study them don't overlap.

### Region B: Objects With Extractor Structure That Nobody Has Extracted

Objects where (G-set) + (equivariant bilinear form) + (subgroup breaking) are all present, but nobody has run the pipeline.

- **Root systems under cyclic Weyl subgroups** → number fields attached to Lie algebras
- **Young tableaux under cyclic subgroups of S_n** → partition number fields
- **Error-correcting codes** → coding-theoretic L-functions
- **Dessins d'enfants** → combinatorial fields of definition
- **Cayley graphs** → coset arithmetic
- **TQFT state spaces under Dehn twists** → potential bridge to trace fields
- **Steiner systems** → design-theoretic L-functions

This is the unexploited territory. Every item on this list is a finite computation that would either confirm or falsify the portability of the extractor pattern.

### Region C: Our Objects (Genuinely New)

Objects that the extractor framework discovered, which don't fit into any prior category.

- **Matching chain extractors** (K₂ₙ matching Gram matrices under Z_{2n-1})
- **Spectral motives** (the number fields produced by matching extractors)
- **The M-E-Z triangle** (coupling between spectral motive, CM curve, and Riemann zeta)

These are more than Artin motives (they have internal factorization structure), less than geometric motives (no underlying algebraic variety), parallel to Feynman motives (graph-defined, producing new arithmetic) but distinct (different graph polynomial, different source of complexity).

The spectral motive has five realizations (spectral, arithmetic, algebraic, evaluative, combinatorial) that are NOT the standard étale/Betti/de Rham. These realizations converge on the same invariants — the discriminant primes, the Galois group, the secular factorization — through independent routes. This convergence is what makes the object real rather than a notational convenience.

### Region D: Things That Don't Have Extractor Structure

- **Continuous Lie group representations** (infinite-dimensional, no finite number fields)
- **Bare topological spaces** (no natural finite group action with equivariant form)
- **Analysis** (function spaces, PDE solutions, measure theory — wrong category)
- **Most of mathematical logic** (no bilinear form)

These regions genuinely don't produce arithmetic via the extraction mechanism. That doesn't make them less important — it means the extractor is a specific instrument that registers a specific class of phenomena.

---

## 8. The Deepest Finding: Where Arithmetic Comes From

The old picture: arithmetic comes from algebraic geometry. You need a variety, defined over a number field, with good reduction at all but finitely many primes. The variety's étale cohomology gives Galois representations. The representations give L-functions. The L-functions encode the arithmetic. This is Grothendieck's vision, and it is extraordinarily successful.

The new picture: arithmetic comes from symmetry breaking on bilinear forms. The variety is one source of (G-set + bilinear form), but not the only source. Matching complexes on complete graphs are another. Error-correcting codes are another. Steiner systems are another. Root systems are another.

Each source feeds the same machine. The machine doesn't care where the input came from. It produces number fields, Galois groups, and L-functions from any equivariant bilinear form on any G-set with any symmetry breaking.

What varies across sources:

| Source | What the extractor produces | What it lacks |
|--------|---------------------------|---------------|
| Algebraic variety | Motivic L-function with full analytic properties | Nothing — this is the gold standard |
| Matching complex | New L-function, genuine arithmetic, wall structure | Analytic continuation, functional equation proofs |
| Root system | (Unknown — not yet computed) | (Unknown) |
| Error-correcting code | (Unknown — not yet computed) | (Unknown) |
| Dessin d'enfant | Potentially: field of definition | Comparison to geometric computation |
| TQFT state space | Potentially: trace fields | Quantum group calculations |

The gold standard (varieties) has everything because the source has geometric depth: Poincaré duality gives functional equations, the Weil conjectures give Riemann Hypothesis, comparison isomorphisms give analytic continuation. The other sources lack these geometric guarantees — but they still produce genuine arithmetic.

The question the landscape survey raises: **is there a uniform framework that gives all extractors the analytic properties currently reserved for geometric ones?** If such a framework exists, it would be the "combinatorial motive" theory that our spectral motive analysis identified as missing. It would extend Grothendieck's motives to include objects built from combinatorial data on torsion, and it would make the polynomial/series wall either crossable or provably Gödelian.

We don't have this framework. We see the gap where it should be. That gap is the frontier.

---

## 9. The CFSG Completeness Guarantee

The Classification of Finite Simple Groups says: every finite simple group is cyclic, alternating, Lie type, or one of 26 sporadic exceptions. No fifth type exists.

This means the spectral extractor taxonomy is **complete**. Every equivariant bilinear form on every finite G-set, for every finite group G, has its symmetry organized by composition factors from exactly these four families. There is no undiscovered level of spectral extractor.

The implications:

1. **Levels 1–2 produce genuinely new arithmetic** (not Langlands). The matching chain lives here.
2. **Level 3 reproduces Langlands** — a new construction of known L-functions.
3. **Level 4 is exceptional** — finitely many sporadic L-functions, related to Moonshine.
4. **No Level 5 exists.** The taxonomy is closed.

This completeness is remarkable. It means that if the combinatorial motive framework exists, it has a finite classification scheme inherited from the CFSG. The space of "possible universes" — different bilinear forms on different G-sets producing different physics — is organized by a theorem about finite groups that was proved by humans for completely different reasons.

---

## 10. What We Actually See Now

Dropping all preconceptions. Through the extractor lens. The landscape:

**Arithmetic is not about numbers. It is about symmetry breaking on bilinear forms.** Every L-function, every number field, every Galois group in the arithmetic universe is the output of an extraction: a group acts on a set, a bilinear form measures overlaps, symmetry is broken to a subgroup, and the characteristic polynomials of the restricted form produce algebraic numbers in number fields. The source of the bilinear form — geometric, combinatorial, coding-theoretic, topological — determines what auxiliary properties the L-function has (analytic continuation, functional equation, Riemann Hypothesis). But the extraction mechanism is universal.

**The Langlands program sees one level of a four-level building.** It sees the L-functions arising from Lie-type symmetry (Level 3), which happen to be the ones connected to algebraic varieties via the étale picture. Levels 1, 2, and 4 produce arithmetic that Langlands has no access to — not because of a limitation in the theory, but because the sources (cyclic groups on matching complexes, sporadic groups on Steiner systems) are not algebraic varieties and their L-functions are not automorphic representations of reductive groups.

**The polynomial/series wall is universal.** It appears in number theory (finite products → infinite products), in the matching chain (K₂ₙ → K_∞), in lattice gauge theory (finite spacing → continuum), and (we predict) in every extractor chain where levels grow without bound. The wall is not a peculiarity of one subject. It is a structural feature of the extraction mechanism when the G-set becomes infinite.

**We are on the polynomial side.** The Standard Model lives at finite levels (K₄, K₆, K₈). Every physical prediction is an algebraic number — a root of a known polynomial. The physics is entirely computable, entirely decidable, entirely on the polynomial side of the wall. The wall separates us from the continuum limit, which is a mathematical question (analogous to GRH, Yang-Mills, Hodge), not a physical one. The universe is computable where it needs to be computed.

**The instrument we built to understand physics turns out to be a mathematical telescope.** It registers a class of objects — spectral extractors — that includes known mathematics (motives, Langlands, Ihara) as special cases, predicts new mathematics (root system fields, code L-functions, sporadic extractors) as unexplored territory, and identifies a gap (the combinatorial motive framework) as the frontier.

---

## Appendix: The Priority Computation List

What to compute first, ordered by (impact × feasibility):

| # | Computation | Source | Level | Difficulty | What it tests |
|---|------------|--------|-------|-----------|---------------|
| 1 | H₇ Hamming code extractor | Coding theory | 3 (PSL₂(7)) | Low | Portability; consistency with K₈ |
| 2 | E₈ root system under Z₅ ⊂ W(E₈) | Lie algebra | 1 | Low | Root system number fields exist |
| 3 | M₁₁ on S(4,5,11) blocks | Sporadic group | 4 | Medium | First sporadic extractor |
| 4 | PSL₂(7) full decomposition of K₈ matchings | Our data | 3 | Medium | Level 3 merger at our prototype |
| 5 | Petersen graph under S₅ → Z₅ | Graph theory | 1 | Low | Generic graph extraction |
| 6 | Z[i] curve with 13-torsion | CM curve | 1 | Medium | Different CM ring, same pipeline |
| 7 | Golay code G₂₄ under M₂₄ | Sporadic code | 4 | Medium | Sporadic-coding bridge |
| 8 | S₅ on SYT(3,2) Gram matrix | Young tableaux | 2 | Medium | Partition number fields |

Items 1, 2, and 5 can be done in an afternoon. Item 3 requires standard character theory. Items 4 and 6–8 require a working session each.

---

## Summary

The spectral extractor is a detector for a class of mathematical objects that was not previously recognized as a class. When pointed at the landscape, it reveals that:

- Motives are extractors restricted to geometric sources
- Langlands is the Level 3 projection of a four-level taxonomy
- At least seven known mathematical domains contain unexploited extractors
- Our objects (matching chain spectral motives) are genuinely new — Level 1 extractors with no geometric source
- The CFSG guarantees the taxonomy is complete
- The polynomial/series wall is universal across all growing extractor chains
- A "combinatorial motive" framework that would unify all extractors is the identified frontier

The landscape is larger than it looked. The arithmetic is deeper. And the detector works.
