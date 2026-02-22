# K₈ Construction Reference
## Definitive Specification — DO NOT DEVIATE

**Date:** February 22, 2026  
**Purpose:** Stop chasing normalization ghosts. Every K₈ computation starts here.

---

## 1. The Three K₆ Matrices (They Are Different!)

There are THREE different 15×15 matrices that all involve "K₆ matchings."
Using the wrong one silently corrupts all downstream results.

| Matrix | Eigenvalue | Where it lives | Direction blocking | Use for |
|--------|-----------|----------------|-------------------|---------|
| **Standalone K₆ Gram** | λ = 3.306 | k6_tools.py | 3D scalar (d₀,d₁,d₂ on K₅ torus) | Higgs mass (R = a₄/a₂² = 0.3722) |
| **Hub sub-Gram** | λ = 2.764 | G₈[hub,hub] | 4D multiset (d₀,d₁,d₂,d₃ incl handle) | K₆→K₈ embedding, Yukawa |
| **Spectral extractor K₆** | λ = 2.764 | spectral_extractor_v2.py | Varies by extractor config | Frobenius, Galois, discriminants |

**RULE: For any computation involving K₈, use the hub sub-Gram vacuum.**
The standalone K₆ is ONLY for standalone Higgs mass computation.

---

## 2. Gram Matrix Construction (K₈)

### Matchings
- 105 perfect matchings of K₈ on vertices {0,1,2,3,4,5,6,7}
- Vertex 7 = handle/hub vertex (connects to all ring vertices)
- Vertices 0-6 = ring vertices (Z₇ cyclic structure)

### Edge directions (Heawood)
```
edge_direction(i, j):
    if i == 7 or j == 7: return 3      # handle edge
    diff = min(|i-j| mod 7, |j-i| mod 7)
    {1→0, 2→1, 3→2}[diff]             # torus directions
```

### Torus multiset
For each matching, collect sorted tuple of torus directions (exclude handle):
```
torus_multiset(m) = sorted([edge_direction(a,b) for (a,b) in m if dir < 3])
```
This gives 10 distinct multisets, each appearing 7 or 14 or 21 times.

### ℤ₃ phase
```
z3(m) = sum(torus_multiset(m)) mod 3
```

### Gram matrix
```
G[i,j] = overlap(i,j) × Re(ω^(z3[j] - z3[i])) × δ(multiset[i], multiset[j])
```
where:
- overlap(i,j) = 8 if i=j, else 2 × |shared edges|
- ω = exp(2πi/3)
- δ enforces same direction multiset (block-diagonal)

### Verified eigenspectrum
- 19 clusters, 105 total
- Null: 5 modes
- Vacuum: λ = 1.9595, multiplicity 6
- Maximum: λ = 24.000

---

## 3. Hub Sub-Gram Matrix

### Construction
```python
hub_indices = [i for i in range(105) if (0,1) in matchings[i]]  # 15 matchings
G_hub = G8[np.ix_(hub_indices, hub_indices)]                     # 15×15 submatrix
```

### Why it differs from standalone K₆
The standalone K₆ Gram matrix uses 3D scalar direction blocking on the K₅ torus.
The hub sub-Gram inherits K₈'s 4D multiset blocking (3 torus + 1 handle).
Every hub matching contains exactly one handle edge, but the handle direction
changes WHICH pairs of hub matchings share the same multiset.

### Verified eigenspectrum
Ground state: λ_hub = 2.7639 (NOT 3.306)

### Hub vacuum embedding
```python
ev_hub, ec_hub = eigh(G_hub)
v_hub = ec_hub[:, 0]              # ground state eigenvector

v_embedded = np.zeros(105)
for ki in range(15):
    v_embedded[hub_indices[ki]] = v_hub[ki]
```

---

## 4. K₆→K₈ Projection (Paper III Theorem 5.1)

### Verified numbers
| Quantity | Value | Source |
|----------|-------|--------|
| ‖Π_vac v_hub‖ | 0.2355 | 23.55% of hub vacuum |
| ρ₁ fraction | 100.0% | Pure ρ₁ |
| ρ₂ fraction | 0.0% | Selection rule |
| ρ₃ fraction | 0.0% | Pure ρ₁ |

### Yukawa hierarchy
| Construction | Hierarchy |
|-------------|-----------|
| Full projection (correct) | **415 : 135 : 1** |
| ρ₁-only (same, since 100% ρ₁) | **415 : 135 : 1** |
| WRONG (standalone K₆) | 2304 : 332 : 1 |

---

## 5. ℤ₇ Action

### Generator
```
vertex i → (i+1) mod 7  for i < 7
vertex 7 → 7             (fixed)
```

### Irrep projectors
```python
zeta7 = exp(2πi/7)
P_k(v) = (1/7) Σ_{g=0}^{6} zeta7^{-kg} × v[σ^g]
```
where σ^g permutes matching indices by the g-fold rotation.

### Conjugacy: ρ_k and ρ_{7-k} are conjugate pairs
- (ρ₁, ρ₆): quadratic residues mod 7
- (ρ₂, ρ₅): blocked by selection rule
- (ρ₃, ρ₄): non-residues

---

## 6. Yukawa Matrix

### Vertex assignments
- Doublet: {0, 1} (hub edge)
- Generation 0: {2, 3}
- Generation 1: {4, 5}
- Generation 2: {6, 7} — contains handle vertex 7 (topologically distinguished)

### Construction
Non-hub matchings (90 of 105) pair doublet vertices with generation vertices:
```
Y[α,β] = Σ v_M × ω^{z3_M}
```
summed over non-hub matchings M coupling generation α to generation β.

### Counting
- Diagonal Y[α,α]: 6 matchings each
- Off-diagonal Y[α,β]: 12 matchings each
- Total: 3×6 + 6×12 = 90 ✓

---

## 7. Error Log (Don't Repeat These)

| Date | Error | Root cause | Time wasted |
|------|-------|-----------|-------------|
| Feb 21 evening | 2304:332:1 instead of 415:135:1 | Used standalone K₆ vacuum (λ=3.306) instead of hub sub-Gram vacuum (λ=2.764) | ~3 hours |
| Feb 19 | 29.88% projection instead of 23.55% | Same root cause — standalone K₆ in extractor | ~2 hours |
| Feb 18 | λ_hub = 2.764 noted as "anomaly" | Didn't recognize it IS the correct K₆ vacuum for embedding | ~1 hour |

**Pattern:** Every time we compute K₆→K₈ embedding, we default to standalone K₆.
**Fix:** Always extract hub sub-Gram from G₈. Never import from k6_tools for embedding.

---

## 8. Quick Self-Test

Any K₈ computation should pass these checks before proceeding:
1. G₈ has 5 null modes, λ_vac = 1.9595 (×6), λ_max = 24
2. Hub sub-Gram ground state: λ = 2.764 (NOT 3.306)
3. Hub vacuum → K₈ vacuum projection: 23.55% (NOT 29.9%)
4. ρ₁ = 100% at vacuum eigenspace (NOT 74.7%)
5. Yukawa: 415:135:1 (NOT 2304:332:1)