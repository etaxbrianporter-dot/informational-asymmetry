# Why Matchings Are Special: The Structural Theorem

**Brian Porter — February 19, 2026**

---

## The Result

Three computations, one structural conclusion.

### The Data

| Object | Bilinear form | [B, Z_n] = 0 | In assoc. scheme? | Galois | Level |
|--------|--------------|---------------|-------------------|--------|-------|
| Fano overlap | B = 2I+J | ✓ | ✓ | trivial | 0 |
| M₁₁ Steiner, Z₃-phased M†M | Z₃ on overlap depth | ✓ | ✓ (scheme closed under ×) | Z₅ (abelian) | 1a |
| K₇ Dirac D†D | Z₃ on edge directions | ✗ | N/A (no equivariance) | S₇ (generic) | 1b |
| K₈ matching Gram G_dir | Z₃ on shared edge directions | ✓ | **✗** | C₂ ≀ C₃ (non-abelian, small) | 2 |

### The Mechanism

The direction-weighted matching Gram matrix G_dir has two properties simultaneously:

1. **G_dir commutes with Z₇** — so it has clean sector decomposition, and eigenvalues are organized by representation theory.

2. **G_dir is NOT in the association scheme** — because the G_dir entry for two matchings sharing 1 edge depends on the *direction class* of that shared edge, not just the count. Shared edges in direction 0 contribute G = 1.0; shared edges in directions 1 or 2 contribute G = −0.5 + 0.866i.

These two properties together are what allows non-abelian Galois groups.

For M₁₁: M†M is built from overlap-depth adjacency matrices A₁, A₂, A₃. These form an *association scheme* — a commutative algebra closed under matrix multiplication. So M†M ∈ span{I, A₁, A₂, A₃} regardless of what phases we assign. The scheme is a prison. Eigenvalues are forced into ℚ(cos(2π/11)), a cyclotomic field with abelian Galois Z₅. By Kronecker-Weber, that's the ceiling.

For K₈ matchings: the direction-weighted Gram matrix splits the s=1 overlap class into 4 subclasses (handle, diff±1, diff±2, diff±3) and the s=2 class into 9 subclasses. The total algebraic structure has dimension ≥ 14, not 4. G_dir lives in this 14-dimensional refined scheme, not the 4-dimensional association scheme. The extra 10 dimensions of freedom allow sector eigenvalues to generate number fields outside any cyclotomic field.

---

## The Precise Obstruction

**Kronecker-Weber theorem:** Every abelian extension of ℚ is contained in a cyclotomic field ℚ(ζ_n).

**Consequence:** If a Gram matrix lies in the association scheme of a Z_n-set, its sector eigenvalues live in ℚ(cos(2π/n)), which has abelian Galois group Z_{(n-1)/2}. Non-abelian Galois is *algebraically impossible* from within the association scheme.

**The matching escape:** The K₈ matching Gram matrix commutes with Z₇ but escapes the association scheme. Its eigenvalues generate a degree-6 field with Galois group C₂ ≀ C₃ = (Z₂)³ ⋊ Z₃, which is non-abelian. This field is NOT contained in ℚ(ζ₇) or any cyclotomic field.

**What breaks the scheme:** Two matchings sharing one edge in direction class 0 have a fundamentally different geometric relationship than two matchings sharing one edge in direction class 2. The direction class records how the shared edge sits on the Heawood surface — it carries *topological* information that the overlap count discards.

---

## The Hierarchy (Final)

**Level 0 — Trivial.** Bilinear form is in association scheme with constant off-diagonal. All sectors degenerate. No number fields. Examples: Fano plane overlap (B = 2I+J), any 2-design overlap.

**Level 1a — Cyclotomic.** Bilinear form is in association scheme, non-constant off-diagonal. Commutes with Z_n. Sector eigenvalues are polynomials in cos(2πk/n) with rational coefficients. Galois group is abelian (subgroup of Z_{(n-1)/2}). Examples: M₁₁ Steiner system with overlap-depth phases, any t-design with scheme-valued phases.

**Level 1b — Generic.** Bilinear form breaks equivariance. No clean sector decomposition. Eigenvalues unconstrained. Galois group is generic (S_d for degree d). Examples: K₇ directed Dirac D†D (Gal = S₇).

**Level 2 — Non-abelian structured.** Bilinear form commutes with Z_n (clean sectors) but is NOT in the association scheme (carries finer invariants). Sector eigenvalues generate non-abelian extensions. Galois group is small, non-abelian, solvable. Examples: K₈ matching Gram matrix (Gal = C₂ ≀ C₃).

**The structural requirement for Level 2:** The bilinear form must encode invariants finer than the association scheme (overlap count), while preserving the cyclic symmetry. This requires a *geometric embedding* that assigns direction classes to edges, and a *combinatorial structure* (matchings) where the geometric data of shared elements is visible.

---

## What Makes Matchings Special

It's not the matching structure alone. It's not the direction structure alone. It's the *interaction*:

1. **Matchings are collections of EDGES**, not vertices. Edges have *direction classes* from the surface embedding. Vertex-sets (like Steiner blocks) don't.

2. **The overlap of two matchings inherits direction data.** When matchings m₁ and m₂ share an edge, that shared edge has a specific direction class. The direction is visible in the Gram matrix but invisible in the association scheme.

3. **The association scheme only sees overlap COUNTS.** For Steiner blocks, M₁₁ acts transitively on each overlap class — all pairs sharing s points "look the same" to the symmetry group. For matchings, Z₇ does NOT act transitively within an overlap class — pairs sharing one direction-0 edge are in a different orbit from pairs sharing one direction-2 edge.

4. **The refined invariant (overlap, direction-signature) preserves equivariance.** Direction classes are Z₇-invariant (because they're defined by |i−j| mod 7). So G_dir commutes with Z₇ even though it discriminates within overlap classes.

The critical chain: **surface embedding → edge direction classes → refinement of matching overlap → Gram matrix escapes association scheme → non-abelian Galois**.

Remove any link and the chain breaks:
- No surface → no directions → back to association scheme → abelian (Level 1a)
- No matchings → no shared edges → directions have nothing to refine → abelian
- No equivariance → no clean sectors → generic (Level 1b)

---

## Implications

1. **The matching chain's arithmetic is genuinely special.** It cannot be replicated by any Steiner system, any t-design, or any other vertex-based combinatorial structure with group action. The non-abelian Galois groups that carry physical content (particle masses, coupling constants) require the edge-based, direction-refined structure that only matchings on surface-embedded graphs provide.

2. **M₁₁ contributes nothing beyond Z₁₁.** The sporadic group's exotic structure is invisible to the extractor because the bilinear form (block overlap) lies in the association scheme. M₁₁ ≅ PSL₂(11) acts on S(4,5,11) with the same arithmetic output as any transitive Z₁₁ action. The sporadic-ness is wasted.

3. **The CFSG Level 4 (sporadic) extractors need a different approach.** To extract non-trivial arithmetic from sporadic groups, we need bilinear forms that escape their association schemes. This likely requires finding sporadic analogues of "edge direction" — geometric structures on Steiner systems that refine the overlap classification while preserving symmetry.

4. **The Kronecker-Weber wall is the abelian/non-abelian boundary.** Equivariant operators in association schemes produce abelian (cyclotomic) arithmetic. Crossing to non-abelian requires escaping the scheme. The matching complex's direction refinement is the specific mechanism that achieves this crossing for the K₂ₙ chain.

---

## The One-Line Summary

**The association scheme is the abelian prison. Direction-refined matching overlaps are the key that opens it.**
