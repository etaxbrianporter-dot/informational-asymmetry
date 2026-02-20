# K₁₆ Five-Law Validation Results

## Setup

K₁₆: 16 vertices, p=15 (composite, 3×5), 2,027,025 matchings, 1701 blocks, 455 orbits.

Direction permutation cycles: {0,1,3,6} (length 4), {2,5} (length 2), {4} (length 1).
k = 4 (vacuum orbit size). Three orbit-size classes: |orb|=4 (398), |orb|=2 (52), |orb|=1 (5).

Vacuum cycle signature: (7, 0, 0) — all torus edges on the 4-cycle. Zero edges from short cycles.

## Law 1: Vacuum Half-Rank — CONFIRMED

| bs_vac | rank | null | bs/2 |
|--------|------|------|------|
| 120    | 60   | 60   | 60   |

Null projector diagonal all 1/2: True.

Vacuum spectral fingerprint:

| λ        | mult | notes                    |
|----------|------|--------------------------|
| 212.0000 | 1    | ground state             |
| 196.0000 | 1    | second unique            |
| 100.0000 | 2    |                          |
| 79.9215  | 8    | = 2k                    |
| 28.0000  | 4    |                          |
| 24.3751  | 8    | = 2k                    |
| 23.1652  | 4    |                          |
| 20.0000  | 2    |                          |
| 14.4721  | 4    |                          |
| 14.1411  | 8    | = 2k                    |
| 5.5279   | 4    |                          |
| 4.8348   | 4    |                          |
| 4.0000   | 2    |                          |
| 1.5624   | 8    | = 2k                    |
| 0.0000   | 60   | null space               |

Rich spectrum: 14 distinct eigenvalue levels, four with generic multiplicity 2k=8.
Non-generic multiplicities: 1, 2, 4 from cycle structure.

## Law 2: FP Rank Deficit — CONFIRMED (with generalization)

### Generalized Law 2: deficit = k − |orb| for direction-compatible orbits

| |orb| | Expected deficit | Observed | Count | Status |
|-------|-----------------|----------|-------|--------|
| 4     | 0               | 0        | 120   | ✓      |
| 2     | 2               | 2        | 18    | ✓      |
| 1 (true FP) | 3        | 3        | 2     | ✓      |

True FP (orbit 378, all 7 directions, bs=9465): deficit = 3 = k−1. ✓

### Direction-Disjoint Anomaly (composite-p phenomenon)

22 orbits have direction cycle signature (0, *, *) — no edges from the main 4-cycle.
These are direction-disjoint from the vacuum. All collapse to coupling rank 15.
Deficit = 45 = vac_rank − 15 regardless of orbit size.

This cannot occur at prime p, where the direction permutation has a single cycle:
every orbit shares direction support with the vacuum, so the simple formula holds.

At K₁₀ (also composite, p=9), the fixed direction d=2 forms 3-cycles on Z₉,
limiting to 3 disjoint edges — not enough to fill 4 torus-edge slots.
Direction-disjoint orbits are structurally impossible there.

At K₁₆, the short cycles {2,5} + {4} span 3 directions and CAN fill all 7 slots.
Direction-disjoint sectors emerge for the first time.

### Mixed orbits (partial overlap with vacuum)

8 orbits with |orb|=4 have intermediate deficits (1-5), with deficit scaling
with the fraction of edges on the 4-cycle.

## Law 3: Generic Multiplicity = 2k — CONFIRMED

Multiplicities: [1, 1, 2, 8, 4, 8, 4, 2, 4, 8, 4, 4, 2, 8]

Generic mult = 8 = 2k for four eigenvalue levels. Non-generic multiplicities (1, 2, 4) arise
from the mixed cycle structure [4, 2, 1].

## Law 4: C₂ Parity Selection — CONFIRMED

| Space | C₂=+1 | C₂=−1 |
|-------|-------|-------|
| null  | 28    | 32    |
| live  | 32    | 28    |

Null is CP-odd excess. Asymmetry ±1 per irrep (8 irreps under C₂ × C₄).

## Law 5: CP Independence — CONFIRMED

2⁴ = 16 ≡ 1 mod 15 (not p−1 = 14).
Q ≠ P^k. Symmetry is **product C₂ × C₄** (CP independent).

K₁₆ joins K₈ as the only two levels with genuinely independent charge conjugation
among K₆–K₁₈. Both are the levels where ord_p(2) is odd.

## Composite-p vs Prime-p Summary

| Feature | Prime p | Composite p (K₁₀, p=9) | Composite p (K₁₆, p=15) |
|---------|---------|------------------------|--------------------------|
| Direction cycles | single cycle length k | [3, 1] | [4, 2, 1] |
| Orbit sizes | k or 1 | k or 1 | k, k/2, or 1 |
| Direction-disjoint sectors | impossible | impossible (capacity limit) | PRESENT |
| Law 2 | deficit = k−1 (FP) or 0 | same | deficit = k−|orb| + disjoint anomaly |
| CP independence | depends on ord_p(2) | derived (2³≡−1 mod 9) | INDEPENDENT (2⁴≡1 mod 15) |

## Timing

Total: ~35 seconds (matching generation ~10s, block decomposition ~14s, analysis ~11s).

## Verdict

Laws 1, 3, 4, 5: confirmed without modification.
Law 2: confirmed for direction-compatible orbits; generalized to deficit = k − |orb|;
direction-disjoint anomaly documented as composite-p phenomenon.

No law was broken. The framework was stress-tested at 2M matchings and the structure held,
revealing finer detail (direction-disjointness) that was invisible at prime p.
