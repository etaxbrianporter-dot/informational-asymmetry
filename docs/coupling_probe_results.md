# Coupling Probe: K₈ Results

## Setup
- K₈: 105 matchings, 4 orbits, vacuum = Orbit 2 (bs=14 per block, 42 total across 3 blocks)
- Non-vacuum orbits: 0 (pure, bs=7×3=21), 1 (mixed, bs=7×3=21), 3 (fixed pt, bs=21×1=21)

---

## Sense 2: Vacuum Block Spectral Structure

Full overlap spectrum O|_vacuum (42×42 matrix):

| Eigenvalue | Algebraic form | Multiplicity |
|---|---|---|
| 48 | 6·nv | 1 |
| 24.485 | 16 + 6√2 | 6 |
| 12 | 12 | 8 |
| 7.515 | 16 − 6√2 | 6 |
| 0 | 0 | 21 |

**Key findings:**
- Only 5 distinct eigenvalues. The irrational pair satisfies t² − 32t + 184 = 0
- The **√2 is the Pfaffian signature** — overlap entries are 2×(shared edges), so the spectral algebra lives over ℚ(√2)
- **21-dimensional null space** = exactly half the vacuum orbit = C(7,2) = size of fixed-point orbit
- Degeneracy pattern: **1, 6, 8, 6, 21** — the 6-fold comes from wreath product (C₃ × C₂ action), the 8-fold matches nv
- K₆ comparison: eigenvalues {12, 8(×5), 2(×4)} — full rank (no null space), simpler structure

---

## Sense 1: Inter-Orbit Coupling

### Universal leading coupling
**σ_max = 24√2 ≈ 33.94 for ALL three non-vacuum orbits.**

This is a topological invariant — the peak coupling amplitude from the vacuum to any excitation sector is identical. Discrimination between sectors lives entirely in the sub-leading structure.

### Orbit 1 coupling (mixed orbit): σ² values reveal same algebraic family

| σ² | Algebraic form | Multiplicity |
|---|---|---|
| 1152 | 128·9 = (24√2)² | 1 |
| 102.63 | 80 + 16√2 | 6 |
| 72 | 8·9 | 8 |
| 57.37 | 80 − 16√2 | 6 |

**Same degeneracy pattern (1, 6, 8, 6) as the vacuum eigenvalues!** The √2-paired structure propagates from the vacuum to the coupling. No null space here (rank = 21 = full).

### Frobenius norms (total coupling strength)

| Coupling | ||O_cross||²_F | Factored | Per matching |
|---|---|---|---|
| Vac ↔ Orbit 0 | 3024 | 42 × 72 = 42 × 8 × 9 | 72 |
| Vac ↔ Orbit 1 | 2688 | 42 × 64 = 42 × 8² | 64 |
| Vac ↔ Orbit 3 (fp) | 3024 | 42 × 72 = 42 × 8 × 9 | 72 |

- Orbits 0 and 3 couple with identical total strength despite having completely different structure (3 blocks × 7 vs 1 block × 21)
- Orbit 1 is weaker by factor 64/72 = **8/9**
- Per-matching coupling: nv² = 64 (Orbit 1) vs nv(nv+1) = 72 (Orbits 0, 3)

### Rank deficiency at fixed point
- Orbits 0, 1: rank 21 (full rank)
- Orbit 3 (fixed pt): **rank 19** — 2 coupling modes are silent

---

## Structural Observations

1. **√2 is the coupling field**: All eigenvalues and singular values live in ℚ(√2). This is the Pfaffian extension — the overlap matrix counts shared edges (integer), but its spectral geometry requires √2.

2. **Universal σ_max vs discriminating sub-structure**: The vacuum couples to all excitation sectors with the same peak amplitude 24√2. The physics (which modes are heavy, which are light) comes from the sub-leading σ distribution. This is analogous to universal gauge coupling at tree level with symmetry breaking at loop level.

3. **Null space = fixed-point sector**: The 21-dimensional null space of O|_vacuum has dimension equal to the fixed-point orbit. This suggests these "silent modes" in the vacuum are precisely the directions that point toward the fixed-point sector (the "all-directions-mixed" configurations that have no coupling to pure vacuum excitations).

4. **8/9 ratio for Orbit 1**: The reduced coupling 64/72 = 8/9 = nv/(nv+1). This is the only orbit with a different Frobenius norm, suggesting a structural distinction between the three non-vacuum sectors.

---

## Next probes

- **C₃ decomposition**: Decompose the 42×21 coupling matrices into irreps of the wreath product C₃ action. The 6-fold degeneracies should split into trivial + nontrivial irreps, revealing individual "channel" couplings.
- **K₁₀ / K₁₂ comparison**: Does the universal σ_max persist? Does the null space dimension always match the fixed-point orbit?
- **Single-block coupling**: Compute O(block_i, block_j) for individual direction blocks within and across orbits to see if the hierarchy emerges at the block level before orbit averaging.
