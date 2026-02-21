# Pfaffian Orientation Probe (Sense 8): Results

## Summary

The Pfaffian sign — the ±1 orientation each matching carries in the Pfaffian expansion — reveals deep structure invisible to all previous probes. Six laws emerge, three of them exact.

---

## The Six Results

### Result 1: Global Net Sign = +1 (All Levels)

| Level | Total matchings | (+) | (−) | Net |
|-------|----------------|-----|-----|-----|
| K₆   | 15             | 8   | 7   | +1  |
| K₈   | 105            | 53  | 52  | +1  |
| K₁₀  | 945            | 473 | 472 | +1  |
| K₁₂  | 10395          | 5198| 5197| +1  |

Since (2n−1)!! is always odd, exact balance is impossible. The net sign is the minimum possible imbalance: exactly +1 at every level tested. The universe of matchings has a slight preference for one orientation.

### Result 2: CP-Even Fraction = 1/p (EXACT LAW)

Decompose the vacuum sign vector into CP-even and CP-odd components under the C₂ involution ×(−1) mod p:

| Level | p  | ‖s₊‖² | ‖s₋‖² | ‖s‖² | ‖s₊‖²/‖s‖² | 1/p     |
|-------|----|--------|--------|------|------------|---------|
| K₆   | 5  | 2      | 8      | 10   | 1/5        | 1/5     |
| K₈   | 7  | 6      | 36     | 42   | 1/7        | 1/7     |
| K₁₀  | 9  | 6      | 48     | 54   | 1/9        | 1/9     |
| K₁₂  | 11 | 10     | 100    | 110  | 1/11       | 1/11    |

**Law: The CP-even component of the Pfaffian sign vector carries exactly 1/p of the total norm-squared.**

The fermion sign is overwhelmingly CP-odd, and becomes more so at higher levels. At K₁₂, over 90% of the sign content is CP-odd.

### Result 3: Sign Preservation Under ×(−1) = 1/p (EXACT LAW)

The vertex map ×(−1) mod p preserves the Pfaffian sign of exactly (2n−1)!!/p matchings:

| Level | Preserved | Flipped | Total   | Fraction | 1/p  |
|-------|-----------|---------|---------|----------|------|
| K₆   | 3         | 12      | 15      | 1/5      | 1/5  |
| K₈   | 15        | 90      | 105     | 1/7      | 1/7  |
| K₁₀  | 105       | 840     | 945     | 1/9      | 1/9  |
| K₁₂  | 945       | 9450    | 10395   | 1/11     | 1/11 |

**Anatomy of preservation**: The preserved count decomposes as:

    preserved = (self-conjugate matchings) + 2 × (same-sign conjugate pairs)

| Level | Self-conj f(n−1) | Same-sign pairs | Total preserved |
|-------|-----------------|-----------------|-----------------|
| K₆   | 3               | 0               | 3               |
| K₈   | 7               | 4               | 15              |
| K₁₀  | 25              | 40              | 105             |
| K₁₂  | 81              | 432             | 945             |

The self-conjugate matchings follow the recurrence f(m) = f(m−1) + 2(m−1)·f(m−2) with f(0) = f(1) = 1, giving the sequence 1, 1, 3, 7, 25, 81, 331, ...

Note the recursive structure: the preserved count at K_{2n} equals the total matching count at K_{2n−2} (15 = total K₆, 105 = total K₈, 945 = total K₁₀). This is (2n−3)!! = (2n−1)!!/p.

### Result 4: Sign Projects Into Null Space at K₈

| Level | Null fraction | Live fraction | Interpretation                        |
|-------|--------------|---------------|---------------------------------------|
| K₆   | 0.00%        | 100.00%       | No null space (n=3 exception)         |
| K₈   | 64.60%       | 35.40%        | Sign is predominantly NULL            |
| K₁₀  | 67.34%       | 32.66%        | Sign is predominantly NULL            |
| K₁₂  | 25.21%       | 74.79%        | Sign is predominantly LIVE            |

At the fermion level K₈, nearly 2/3 of the Pfaffian sign lives in the null space — the "silent" sector that the overlap Hamiltonian cannot see. The fermionic orientation is mostly a hidden variable.

At K₁₂, this reverses: the sign is 3/4 live. The null-space fraction is not a simple function of p.

Eigenspace-resolved (K₈): O·s puts 91% of its weight into the λ = 24.49 eigenspace (the 2k = 6 generic multiplicity eigenspace). The overlap matrix acting on the sign vector concentrates into a single irrep.

### Result 5: Assortative Sign Mixing

Within vacuum direction blocks, same-sign matchings share more edges than opposite-sign matchings:

| Level | Block | ⟨O⟩(+,+) | ⟨O⟩(−,−) | ⟨O⟩(+,−) | Pattern        |
|-------|-------|-----------|-----------|-----------|----------------|
| K₈   | 0     | 1.29      | 1.33      | 1.17      | ASSORTATIVE    |
| K₈   | 1     | 1.20      | 1.21      | 1.25      | DISASSORTATIVE |
| K₈   | 2     | 1.36      | 1.47      | 1.08      | ASSORTATIVE    |
| K₁₂  | 0     | 3.45      | 3.73      | 2.20      | ASSORTATIVE    |

Majority pattern: assortative. Same-orientation fermions prefer to overlap. This is a combinatorial echo of the Pauli principle — matchings with the same Pfaffian sign cluster together in the overlap landscape.

### Result 6: No Orbit Is Sign-Balanced at K₈

| Level | Balanced orbits | Total orbits |
|-------|----------------|--------------|
| K₆   | 1/2            | Vacuum balanced, FP not |
| K₈   | 0/4            | **None** balanced       |
| K₁₀  | 4/12           | Mixed                  |
| K₁₂  | 0/26           | **None** balanced       |

At K₈ and K₁₂ (prime p), every orbit carries a net Pfaffian contribution. No orbit is "sign-neutral." At K₁₀ (composite p = 9), four orbits achieve exact balance — consistent with the pattern that composite p allows degeneracies prime p does not.

---

## Structural Identification

### Results 2 and 3 are the same law seen from two sides:

- **Result 2** (algebraic): decompose the sign vector into C₂ eigenspaces → CP-even fraction = 1/p
- **Result 3** (combinatorial): count sign-preserved matchings under ×(−1) → preserved fraction = 1/p

Both measure the same quantity: the "self-conjugate content" of the Pfaffian orientation. The vertex negation involution has exactly 1/p worth of sign-preserving matchings, and the sign vector's CP-even projection carries exactly 1/p of the norm-squared.

### The 1/p law connects to the tower structure:

The preserved count at K_{2n} = (2n−3)!! = total matchings at K_{2n−2}. This means **the sign-neutral content of each level equals the full content of the level below.** The fermion level K₈ has sign-neutral content equal to the Higgs level K₆. Each level's "symmetric core" is literally the previous level.

---

## Physics Interpretation

1. **CP violation is structural, not parametric.** The Pfaffian sign is (p−1)/p CP-odd at every level. This is not a tunable parameter — it's forced by the combinatorics. At K₈ (the fermion level), the CP-odd fraction is 6/7 ≈ 85.7%.

2. **The fermion sign is mostly hidden.** At K₈, 65% of the sign content lives in the null space — the sector invisible to the overlap Hamiltonian. The "observed" fermion physics (live sector) sees only 35% of the orientation structure. The majority of the fermionic nature is "dark."

3. **Sign assortativity = Pauli clustering.** Same-sign matchings preferentially overlap, creating clusters in the matching landscape. This is a pre-quantum version of the exclusion principle: the Pfaffian sign creates a natural partition of the matching space that correlates with the overlap geometry.

4. **The recursive 1/p structure** means each K₂ₙ level encodes the full matching structure of K_{2n−2} in its sign-neutral sector. The chain K₄ → K₆ → K₈ → K₁₀ → ... is not just a sequence of independent levels — each level's orientation structure contains the previous level as its symmetric core.

---

## Computational Infrastructure

All results computed by `pfaffian_probe.py`, which runs K₆ through K₁₂ in ~30 seconds. The probe implements seven sub-analyses:
- 8a: Orbit sign distribution
- 8b: Null/live eigenspace projection
- 8c: Sign × CP correlation
- 8d: Sign-stratified overlap matrices
- 8e: Signed overlap (O·s analysis)
- 8f: Sign flow under vertex maps
- 8g: Direction × sign decomposition

---

## Open Questions

A. **Why is the null fraction ~2/3 at K₈ and K₁₀ but ~1/4 at K₁₂?** The null fraction rational forms are 104/161 (K₈) and 200/297 (K₁₀). Is there a closed-form expression in terms of p and n?

B. **What determines the same-sign conjugate pair count?** We have 0, 4, 40, 432 same-sign pairs at K₆, K₈, K₁₀, K₁₂. Is there a combinatorial formula?

C. **Does the global net sign remain +1 for all K₂ₙ?** If so, this constrains the Pfaffian of the complete graph's adjacency matrix.

D. **Sense 5 (swap graph) × Sense 8 (Pfaffian sign)**: Does the swap graph respect the sign partition? If sign-preserving swaps connect same-sign matchings, the swap graph has a natural Z₂-grading that would be the dynamical analogue of the Pfaffian structure.
