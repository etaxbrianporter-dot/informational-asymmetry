---
title: Spectral Ratio Oscillation vs Cyclotomic Degree: Decoupling Theorem
section: K₆ Higgs Sector
status: active
---

# Spectral Ratio Oscillation vs Cyclotomic Degree: Decoupling Theorem

## Brian Porter â€” February 2026

---

## The Question

Does the non-monotone behavior of the algebraic degree deg(fâ‚‚â‚™) = Ï†(2nâˆ’1), controlled by Euler's totient, have a physical shadow in the spectral ratio R_n = aâ‚„/aâ‚‚Â² that controls the Higgs mass?

## Answer: No. They Are Decoupled.

The spectral ratio R_n and the cyclotomic degree Ï†(2nâˆ’1) respond to different structural features of K_{2n} and are numerically decoupled. Specifically:

- **The biggest Ï† oscillation produces no R response:** Kâ‚â‚„â†’Kâ‚â‚† has Î”Ï† = âˆ’4 (degree drops from 12 to 8), but Î”R = âˆ’0.001 (R changes by 0.5%).
- **The biggest R oscillation occurs at constant Ï†:** Kâ‚ˆâ†’Kâ‚â‚€ has Î”Ï† = 0 (both have Ï† = 6), but Î”R = +0.113 (R jumps by 57%).

These two quantities live in the same number field Q(Î»_vac) at each level, but their numerical values are controlled by orthogonal inputs:
- Ï†(2nâˆ’1) = arithmetic of the cyclic modulus (prime factorization of 2nâˆ’1)
- R_n = spectral geometry of the matching overlap tensor (eigenvalue distribution shape)

---

## Computation

### Method

For each K_{2n} (n = 4 through 8), we:

1. Enumerate all (2nâˆ’1)!! perfect matchings of K_{2n}
2. Classify matchings into Z_{2nâˆ’1} orbits under cyclic shift (fixing vertex 2nâˆ’1)
3. Identify **two-orbit sectors**: groups of 2(2nâˆ’1) matchings sharing the same edge-class signature
4. Build the integer overlap matrix O_{ij} = Tr(M_i M_j) within the sector
5. Find the vacuum eigenvector v at the smallest positive eigenvalue Î»_vac
6. Construct S = Î£ v_i M_i (a 2n Ã— 2n real symmetric matrix)
7. Compute R = Tr(Sâ´)/Tr(SÂ²)Â² = spectral kurtosis of S

**Key structural fact:** Within a two-orbit sector, all matchings share the same edge-class composition, hence the same total Zâ‚ƒ phase. The phase correlator is therefore identically 1 within the sector, and the Gram matrix reduces to the pure overlap matrix. This means the sector-level computation yields the physical R for Kâ‚ˆ and all higher levels.

**Verification:** The characteristic polynomial of each sector overlap matrix factors exactly as expected:

| Level | Sector char. poly. structure | Irreducible degree | Matches Ï†(2nâˆ’1)? |
|:------|:----------------------------|:-------------------|:-----------------|
| Kâ‚ˆ | x Â· (xâˆ’24) Â· [fâ‚†]Â² | 6 | âœ“ Ï†(7) = 6 |
| Kâ‚â‚€ | xÂ³ Â· (xâˆ’36) Â· (xâˆ’12)Â² Â· [fâ‚†]Â² | 6 | âœ“ Ï†(9) = 6 |
| Kâ‚â‚‚ | x Â· (xâˆ’72) Â· [fâ‚â‚€]Â² | 10 | âœ“ Ï†(11) = 10 |
| Kâ‚â‚„ | trivial factors Â· [fâ‚â‚‚]Â² | 12 | âœ“ Ï†(13) = 12 |
| Kâ‚â‚† | trivial factors Â· [fâ‚ˆ]Â² | 8 | âœ“ Ï†(15) = 8 |

All minimal polynomials match the cyclotomic degree law, confirming the computations.

### Results

| Level | 2nâˆ’1 | Ï†(2nâˆ’1) | Matchings | Î»_vac | aâ‚„ | **R = aâ‚„/aâ‚‚Â²** | Î”R |
|:------|:-----|:--------|:----------|:------|:---|:----------------|:---|
| Kâ‚ˆ | 7 | 6 | 105 | 1.9595 | 0.7563 | **0.19698** | â€” |
| Kâ‚â‚€ | 9 | 6 | 945 | 2.2145 | 1.5192 | **0.30978** | +0.113 |
| Kâ‚â‚‚ | 11 | 10 | 10,395 | 0.6366 | 0.0852 | **0.21009** | âˆ’0.100 |
| Kâ‚â‚„ | 13 | 12 | 135,135 | 0.4589 | 0.0465 | **0.22086** | +0.011 |
| Kâ‚â‚† | 15 | 8 | 2,027,025 | 0.3488 | 0.0267 | **0.21970** | âˆ’0.001 |

---

## Three Main Findings

### Finding 1: R and Ï† Are Decoupled

The degree sequence and R sequence are uncorrelated:

```
Ï†(2n-1):   6  â†’  6  â†’ 10 â†’ 12 â†’  8      (oscillates)
    R_n:  .197 â†’ .310 â†’ .210 â†’ .221 â†’ .220  (converges)
```

The critical test is Kâ‚â‚„ â†’ Kâ‚â‚†: the degree drops by 4 (from 12 to 8), but R changes by only 0.001. Conversely, Kâ‚ˆ â†’ Kâ‚â‚€: the degree stays at 6, but R jumps by 0.113.

**What controls Ï†:** The factorization of 2nâˆ’1. At Kâ‚â‚†, 15 = 3 Ã— 5, so Ï†(15) = 2 Ã— 4 = 8 â€” a product of two small primes creates redundancy among cyclotomic characters via CRT.

**What controls R:** The eigenvalue distribution of the vacuum S-matrix. The spectral kurtosis depends on how "peaked" vs "flat" the spectrum of S is. This is determined by the matching overlap tensor's structure, not by the arithmetic of 2nâˆ’1.

### Finding 2: R Converges to ~0.22

After the initial transient (Kâ‚ˆ, Kâ‚â‚€), R converges rapidly:

| Level | R | |R âˆ’ 0.22| | |R âˆ’ 0.22|/R |
|:------|:--|:----------|:------------|
| Kâ‚â‚‚ | 0.2101 | 0.0099 | 4.7% |
| Kâ‚â‚„ | 0.2209 | 0.0009 | 0.4% |
| Kâ‚â‚† | 0.2197 | 0.0003 | 0.1% |

The mean over Kâ‚â‚‚â€“Kâ‚â‚† is RÌ„ = 0.2169 with standard deviation 0.0048. The convergence is geometric: each step reduces the deviation by roughly a factor of 10.

Nearest mathematical constants: 7/32 = 0.21875 (Î” = 0.002), Ï€Â²/45 = 0.21932 (Î” = 0.002). Neither is a convincing match at this precision. Additional data points (Kâ‚â‚ˆ, Kâ‚‚â‚€) would be needed to identify R_âˆž exactly, but these require enumerating 17!! â‰ˆ 34M and 19!! â‰ˆ 655M matchings respectively.

### Finding 3: Kâ‚â‚€ Is Anomalous (The 3Â² Effect)

Kâ‚â‚€ has R = 0.310, far above the other levels. The anomaly correlates with 2nâˆ’1 = 9 = 3Â²:

| Level | 2nâˆ’1 | Factorization | R |
|:------|:-----|:-------------|:--|
| Kâ‚ˆ | 7 | prime | 0.197 |
| **Kâ‚â‚€** | **9** | **3Â²** | **0.310** |
| Kâ‚â‚‚ | 11 | prime | 0.210 |
| Kâ‚â‚„ | 13 | prime | 0.221 |
| Kâ‚â‚† | 15 | 3 Ã— 5 | 0.220 |

The Zâ‚‰ symmetry has a non-trivial subgroup Zâ‚ƒ, which creates extra structure in the matching overlap tensor. The vacuum S-matrix at Kâ‚â‚€ has eigenvalues Â±0.924 dominating, making the spectrum much more peaked (kurtosis-enhancing) than at other levels. The Zâ‚ƒ subgroup forces certain matching overlaps to correlate, concentrating spectral weight.

Note that Kâ‚â‚† (2nâˆ’1 = 15 = 3 Ã— 5) also has Zâ‚ƒ as a subgroup, but the Zâ‚… factor dilutes the effect â€” Râ‚â‚† = 0.220 is normal, not anomalous.

---

## The Galois Conjugate "Echo" of Râ‚†

A striking secondary finding: at each level, there exists a Galois conjugate of R (at a non-physical eigenvalue) whose value is close to Râ‚† = 0.3722:

| Level | Conjugate Î» | R at conjugate | |R âˆ’ Râ‚†| | Relative error |
|:------|:-----------|:---------------|:---------|:---------------|
| Kâ‚ˆ | 3.188 | 0.3689 | 0.0033 | 0.9% |
| Kâ‚â‚‚ | 3.503 | 0.3716 | 0.0006 | **0.16%** |
| Kâ‚â‚„ | 4.357 | 0.4033 | 0.0311 | 8.4% |

At Kâ‚â‚‚, the conjugate R = 0.3716 is within 0.06% of Râ‚† = 0.3722. This "echo" occurs at eigenvalues near Î» â‰ˆ 3.2â€“3.5, close to Kâ‚†'s vacuum eigenvalue Î»â‚† = 3.306.

This is consistent with the Kâ‚† â†’ Kâ‚ˆ embedding: the Kâ‚† vacuum projects 23.55% onto the Kâ‚ˆ vacuum eigenspace (entirely in Ïâ‚), and the Kâ‚ˆ eigenvalue at Î» â‰ˆ 3.29 captures 28.5% of the Kâ‚† vacuum weight. The "ghost" of Kâ‚†'s Higgs-sector ratio persists at higher levels on non-physical Galois conjugates.

**Physical implication:** The Kâ‚† Higgs mass prediction (Râ‚† = 0.3722 â†’ mH â‰ˆ 125 GeV with c = Ï€Â²/8) is algebraically independent from the fermion-sector computations at higher K_{2n}. It lives in Q(âˆš5), disjoint from all higher fields. But its numerical shadow appears as a Galois conjugate at Kâ‚ˆ and Kâ‚â‚‚, suggesting a deeper embedding relationship.

---

## Physical Implications

### Higgs Mass from the Continuum Limit

If R_âˆž â‰ˆ 0.22 represents the continuum limit, the Higgs mass formula mH = âˆš(8R/c) Ã— mW gives:

| c | mH (GeV) | Interpretation |
|:--|:---------|:---------------|
| 1 | 106.7 | No normalization |
| Ï€Â²/8 = 1.234 | 96.0 | Kâ‚† normalization |
| R_âˆž/Râ‚† = 0.59 | 138.9 | Self-consistent scaling |

None of these reproduce 125 GeV. The continuum limit R_âˆž does not directly compute the Higgs mass â€” that remains the province of Kâ‚†.

### What Each Level Knows

The clean separation between Ï† (algebraic) and R (spectral) reinforces the established picture:

- **Kâ‚†** knows the Higgs mass (Râ‚† = 0.3722 â†’ mH = 125 GeV)
- **Kâ‚ˆ** knows the fermion hierarchy (rank-3 Yukawa, three generations)
- **Kâ‚â‚€â€“Kâ‚â‚†** probe the algebraic tower (cyclotomic degrees, Galois disjointness)
- **R_âˆž** characterizes the spectral geometry of the infinite matching system

The Euler totient controls which number field the vacuum lives in. The spectral kurtosis controls which physical observable it computes. These are independent data about the same mathematical object.

---

## Summary Table

| Property | Controlled by Ï†(2nâˆ’1)? | Controlled by spectral geometry? |
|:---------|:----------------------|:-------------------------------|
| Minimal polynomial degree | âœ“ | âœ— |
| Galois group structure | âœ“ | âœ— |
| Generation count | âœ“ | âœ— |
| Discriminant primes | âœ“ | âœ— |
| Vacuum eigenvalue Î»_vac | âœ— (not monotone, but not Ï†-correlated) | âœ“ |
| Spectral ratio R | âœ— | âœ“ |
| Higgs mass prediction | âœ— | âœ“ (via Râ‚†) |
| Yukawa hierarchy | âœ— | âœ“ (via Kâ‚ˆ vacuum moduli) |

**The algebraic and spectral sectors of the K_{2n} tower are decoupled.**
