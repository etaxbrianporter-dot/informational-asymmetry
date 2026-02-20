# Inter-Level Coupling: The Rank-1 Embedding Theorem

## The Natural Map

There is a canonical embedding of M(K₂ₙ) into M(K₂ₙ₊₂): a matching M of K₂ₙ
maps to M ∪ {(2n, 2n+1)} in K₂ₙ₊₂. This pairs the two new vertices with each
other, leaving the original matching untouched. The embedded sector is exactly
the set of K₂ₙ₊₂ matchings containing the "top edge."

This is the unique structure-preserving embedding: it adds no new edges among
the original vertices, and the top edge is the only way to absorb both new
vertices without disturbing the existing pairing.

## The Rank-1 Perturbation Theorem

**Theorem.** Let O_nat be the overlap matrix of M(K₂ₙ) and O_emb be the overlap
matrix of the embedded sector in K₂ₙ₊₂, computed using K₂ₙ₊₂ overlaps. Then:

    O_emb = O_nat + 2J

where J is the (2n-1)!! × (2n-1)!! all-ones matrix.

**Proof.** Every embedded matching contains the top edge (2n, 2n+1). For any two
embedded matchings Mᵢ ∪ {top} and Mⱼ ∪ {top}, their edge overlap in K₂ₙ₊₂ is:

    |Mᵢ ∩ Mⱼ| + 1    (the +1 is the shared top edge)

The overlap matrix entry is 2|Mᵢ ∩ Mⱼ| + 2 when i ≠ j. When i = j, the diagonal
is (2n+2) = (2n) + 2. In both cases: O_emb[i,j] = O_nat[i,j] + 2.  ∎

**Corollary.** The eigenvalues of O_emb are:

    λ₀_emb = λ₀_nat + 2·(2n-1)!!    (ground state, uniform eigenvector)
    λᵢ_emb = λᵢ_nat                   (all other eigenvalues, i ≥ 1)

since J has eigenvalue (2n-1)!! on the uniform vector and 0 on its orthogonal
complement. The eigenvectors are identical.

## Computational Verification

| Embedding | n_matchings | ||O_emb - O_nat - 2J|| | Ground state shift | All other shifts |
|-----------|-------------|------------------------|--------------------|------------------|
| K₆ → K₈  | 15          | 0.0                    | +30 = 2×15         | exactly 0        |
| K₈ → K₁₀ | 105         | 0.0                    | +210 = 2×105       | exactly 0        |
| K₁₀→ K₁₂ | 945         | 0.0                    | +1890 = 2×945      | exactly 0        |

Machine-precision zeros. The theorem is exact, not approximate.

## Consequences

### 1. Topological Protection of Internal Geometry

The entire internal structure of each level — null space, half-rank, CP decomposition,
irrep content, all five laws — is invariant under inter-level embedding. The perturbation
2J acts only on the uniform (ground state) direction and leaves all relative structure
untouched. The vacuum geometry at K₂ₙ is the same whether computed natively or as a
subspace of K₂ₙ₊₂.

### 2. Geometric Realization of Galois Disjointness

The Galois disjointness theorem (number fields at different levels are algebraically
independent) has a geometric counterpart: the overlap spectra are rank-1 perturbations.
The only coupling is through the scalar ground-state shift. The spectral structure
that encodes physics (eigenvalue ratios, null space dimensions, parity selection) lives
entirely in the orthogonal complement of the uniform vector, where the embedding acts
as the identity.

### 3. The Aperture as the Inter-Level Coupling

The top edge (2n, 2n+1) is the aperture — the unique degree of freedom connecting
adjacent levels. In the NCG interpretation, this is the piece that connects the base
manifold to the internal space. The rank-1 theorem says this connection is maximally
simple: a single scalar coupling (the shared edge count) with no internal structure.

The aperture contributes to the cosmological constant term in the spectral action
(the a₀ Seeley-DeWitt coefficient, which is a trace). The rank-1 nature of the
inter-level coupling is consistent with this: traces see only the uniform direction.

### 4. Orbit Scrambling Under Torus Change

While the spectrum is protected, the orbit decomposition is not:

| Embedding | Vacuum → Vacuum fraction | Total orbits hit |
|-----------|--------------------------|------------------|
| K₆ → K₈  | 3/10 = 30%               | all 4            |
| K₈ → K₁₀ | 2/42 = 4.8%              | 11 of 12         |
| K₁₀→ K₁₂ | 1/54 = 1.9%              | 19 of 26         |

The K₂ₙ vacuum scatters across nearly all K₂ₙ₊₂ orbits, with the vacuum-to-vacuum
fraction dropping rapidly. This is because orbits depend on direction labels mod p,
and p changes between levels. The same edges get relabeled when the torus changes
from Z_p to Z_{p+2}.

The algebra (eigenvalues) remembers. The geometry (orbit classification) forgets.

### 5. Direction Relabeling

The same edge set acquires different torus directions at different levels:

    Edge set with K₆ dirs: (0, 0)     →  K₈ dirs: (0, 0, 0)
    Edge set with K₆ dirs: (0, 1)     →  K₈ dirs: (0, 1, 1)

This is not a bug — it's the mechanism by which each level contributes genuinely
independent geometric content. The torus Z_p at level K₂ₙ (p = 2n-1) is not a
subgroup of Z_{p+2} at level K₂ₙ₊₂. The direction classification starts fresh
at each level, and the rank-1 theorem guarantees this relabeling cannot corrupt
the internal spectral structure.

## Cross-Level Coupling Matrix

The overlap between the embedded K₂ₙ sector and the native K₂ₙ₊₂ vacuum was also
computed. This "cross-level coupling" has structure beyond rank-1:

| Coupling            | Matrix size    | Rank | Singular value pattern          |
|---------------------|---------------|------|---------------------------------|
| emb(K₆) × vac(K₈)  | 15 × 42      | 10   | 35.7, 12.4, 11.3(×5), 8.9(×2), 8.5 |
| emb(K₈) × vac(K₁₀) | 105 × 54     | 21   | 101.7, 30.3(×2), 28.4(×4), ...  |
| emb(K₁₀)× vac(K₁₂) | 945 × 110    | 36   | 408.2, 123.6(×2), 97.0(×2), ... |

The cross-level coupling rank is NOT trivial — it encodes how the embedded
sector projects onto the higher-level vacuum. This matrix carries information
about how Higgs-level physics (K₆ vacuum) couples to fermion-level physics
(K₈ vacuum), mediated through the shared matching structure. The singular value
multiplicities show remnants of the cyclic symmetry at the higher level.

## Summary

The inter-level embedding is exactly as simple as it could possibly be: a rank-1
perturbation that shifts the ground state by 2×(number of matchings) and leaves
all internal structure invariant. This is the matching-algebraic version of
"the levels don't interfere." Each K₂ₙ contributes independent physics because
the embedding preserves the spectral structure that encodes that physics.

The simplicity of this result — an exact rank-1 theorem with a one-line proof —
is itself significant. It means the tower of K₂ₙ levels is not a perturbation
series or an approximation scheme. Each level is an exact, self-contained algebraic
structure connected to adjacent levels by a single scalar. The inter-level coupling
is the aperture, and the aperture is a trace.
