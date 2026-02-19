---
title: Q3 Resolution: Algebraic Status of R₈ = 0.196975716486...
section: K₆ Higgs Sector
status: active
---

# Q3 Resolution: Algebraic Status of R₈ = 0.196975716486...

## Summary

**R₈ is an algebraic number of degree 6 over Q, but its minimal polynomial 
has impractically large integer coefficients. There is no clean closed-form 
expression.**

## What We Established

### 1. The Gram Matrix is Integer-Valued
The 105×105 Gram matrix G has entries in {0, 2, 4, 8} only. All eigenvalues 
are algebraic numbers.

### 2. The Vacuum Eigenvalue Has a Known Minimal Polynomial
λ_vac satisfies the degree-6 integer polynomial:

  **λ⁶ − 44λ⁵ + 720λ⁴ − 5648λ³ + 22512λ² − 43456λ + 31808 = 0**

The 6 roots are {1.9595, 3.1878, 4.1631, 6.8327, 10.0567, 17.8001}, which 
form a single Galois orbit. These appear as eigenvalues of the full Gram 
matrix with multiplicity 6 each (the three ℤ₇ irreps × 2D real embedding).

### 3. The Active Subspace Has Integer Structure  
At any pure ρ_k direction, only 14 matchings are active, all sharing the 
same direction vector d = (1,−1,1,0) and phase ζ = −1/2 − i√3/2.

The 14×14 sub-Gram matrix is **integer-valued** with entries in {0, 2, 4, 8}.
Its characteristic polynomial factors as:

  **x · (x − 24) · [λ⁶ − 44λ⁵ + 720λ⁴ − 5648λ³ + 22512λ² − 43456λ + 31808]²**

confirming that λ_vac is an eigenvalue of this integer matrix with multiplicity 2.

### 4. R₈ Reduces to Spectral Kurtosis of an 8×8 Matrix
Since all active matchings have the same direction and phase, BZ averaging 
is trivial. The Dirac operator simplifies to D = ζ·e^{ik·d}·S where 
S = Σᵢ vᵢ Mᵢ is a fixed 8×8 real symmetric matrix.

Then:
- a₂ = Tr(S²) = Σ μᵢ² = λ_vac
- a₄ = Tr(S⁴) = Σ μᵢ⁴  
- **R₈ = Σμᵢ⁴ / (Σμᵢ²)²** = spectral kurtosis of S

### 5. R₈ Has Degree 6 Over Q
Computing R₈ at all 6 Galois conjugates of λ_vac gives 6 distinct values:

| λ_k        | R(λ_k)         |
|-------------|-----------------|
| 1.9595...   | 0.196975716486  |
| 3.1878...   | 0.368917445321  |
| 4.1631...   | 0.240731354047  |
| 6.8327...   | 0.252346587873  |
| 10.0567...  | 0.280403342188  |
| 17.8001...  | 0.183703754702  |

All distinct → R₈ has degree 6 over Q (same as λ_vac).

### 6. The Minimal Polynomial Has Very Large Coefficients
PSLQ searches with maxcoeff up to 10²⁰ failed to find the degree-6 
minimal polynomial of R₈. All purported relations at lower degree 
(4 or 5) are spurious.

The symmetric functions of the 6 R values give approximate coefficients:
- e₁ = Σ Rᵢ ≈ 1.52308
- e₂ = Σ RᵢRⱼ ≈ 0.95543
- e₃ ≈ 0.31618
- e₄ ≈ 0.05824
- e₅ ≈ 0.00567
- e₆ = Π Rᵢ ≈ 0.000227

These are **not** small rationals. The common denominator exceeds 10²⁰, 
making the integer-coefficient minimal polynomial impractical to write.

## Why the Coefficients Are Large

The eigenvector of the 14×14 integer Gram matrix at λ_vac involves the 
adjugate of (G_sub − λI), which is a 14×14 matrix of degree-13 polynomials 
in λ. The quartic form a₄ = Σ vᵢvⱼvₖvₗ Tr(MᵢMⱼMₖMₗ) then involves 
4th powers of these polynomials, producing rational functions of λ with 
degree ~52 numerator/denominator (before reduction mod the degree-6 
minimal polynomial). The resulting expression for R in Q(λ) has rational 
coefficients whose numerators and denominators reflect these high-degree 
intermediate expressions.

## High-Precision Value

To 50 significant digits:

  **R₈ = 0.19697571648572320585882537334051356582159690742413**

Computed from:
- a₂ = 1.959511692168087010271145404248703875303
- a₄ = 0.756324915061916048768917223800601492499

Using:
- 14×14 integer sub-Gram matrix
- Eigenvector at λ_vac via inverse iteration (80-digit mpmath)
- S matrix eigenvalues via mpmath.eigsy

## Physical Interpretation

R₈ = 0.19698... is **not** a simple fraction, **not** expressible via √7 
or cos(2π/7) in any clean way, and **not** related to R₆ = 0.37221... by 
a simple algebraic factor. The spectral kurtosis of the vacuum Dirac 
operator is an irreducible degree-6 algebraic number determined by the 
combinatorial structure of K₈ perfect matchings.

This confirms the picture from the main analysis: K₈ spectral data lives 
in a genuinely different algebraic extension from K₆, and the two systems 
compute independent physical parameters (Yukawa couplings vs Higgs mass).

## Files

- k8_gram.npz: Saved Gram matrix, eigenvectors, phases, directions
- k8_final_R.py: High-precision R computation via 14×14 integer matrix
- k8_R_all6.py: R at all 6 Galois conjugates
- k8_factor.py: Factorization of degree-15 characteristic polynomial
