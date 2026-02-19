# Informational Asymmetry Program: Documentation & Verification Plan

## Brian Porter — February 2026

---

## I. Documentation Architecture

### Paper Series

| Paper | Title | Core Claim | Status |
|-------|-------|-----------|--------|
| **I** | Lorentzian Signature from Spectral Action on K₄ | Three axioms → K₄, V₄, Z₃ flux, signature (1,3) | Draft complete |
| **II** | Cosmic Birefringence from K₄ Torsion | Five parameter-free CMB predictions | Draft complete |
| **III** | Fermion Mass Hierarchy from K₈ on Genus-2 | Three generations, 415:135:1 hierarchy, ρ₂ selection rule | Draft complete |

### Supporting Documents

| Document | Content | Audience |
|----------|---------|----------|
| **Executive Summary** | One-paragraph claim → five-axiom chain → status table | General physics |
| **Spectral Budget Identity** | Johnson scheme → overlap eigenvalues → polynomial reduction | Mathematical |
| **Galois Disjointness** | Number fields at K₆, K₈, K₁₀ are linearly disjoint over ℚ | Mathematical |
| **Classification Theorem** | FORCED/SILENT/REDUCED taxonomy of open problems | Program-internal |

### Companion Verification Repository

**This document.** A self-contained Python package that independently verifies every computational claim in Papers I–III and the supporting documents.

---

## II. Verification Suite Structure

```
verification_suite/
├── PLAN.md                          # This file
├── README.md                        # How to run, what each test checks
├── run_all.py                       # Master runner with PASS/FAIL summary
├── requirements.txt                 # numpy, scipy, sympy (standard only)
│
├── tests/
│   ├── test_01_saturation.py        # (2n-1)!! = I(n), n=3 uniqueness
│   ├── test_02_k4_matching.py       # V₄ algebra, six-face equivalence, order-5
│   ├── test_03_k4_pfaffian.py       # 15 subgraphs, Pf classification, 13/15
│   ├── test_04_k4_eigenvalues.py    # D² spectrum, R = 7+4√3, s_crit
│   ├── test_05_k4_hessian.py        # Lorentzian signature emergence
│   ├── test_06_k4_heat_kernel.py    # a₀=20, b₀=-4, parity ratio -1/5
│   ├── test_07_k4_commutator.py     # Arrow of time: D≠D*, Tr(C)=0
│   ├── test_08_k6_matchings.py      # 15 matchings, 120:90:15 overlaps
│   ├── test_09_k6_gram.py           # Gram matrix, eigenspectrum, a₂, a₄, R
│   ├── test_10_k6_higgs.py          # Higgs mass formula, moduli distribution
│   ├── test_11_k6_bz.py            # BZ k-independence, gauge kinetic = a₂
│   ├── test_12_birefringence.py     # A₀*, b₀/a₀, α₀ bound
│   ├── test_13_johnson_scheme.py    # Overlap eigenvalues, all levels K₆–K₁₆
│   ├── test_14_budget_identity.py   # λ_mid · d_phys = 2(n-1) · λ_max
│   ├── test_15_physical_projector.py # Polynomial projector, hub memory
│   ├── test_16_k8_structure.py      # 105 matchings, Z₇, genus-2
│   ├── test_17_k8_vacuum.py         # 6D vacuum, ρ₁⊕ρ₂⊕ρ₃ decomposition
│   ├── test_18_k8_yukawa.py         # 415:135:1 hierarchy, ρ₂ selection rule
│   ├── test_19_k10_k12.py          # Higher levels, degree stabilization
│   ├── test_20_galois.py           # Number field disjointness
│   ├── test_21_d4_uniqueness.py     # Sub-Pfaffian frustration theorem
│   └── test_22_bug_demonstration.py # The .real bug: showing both bugs needed
│
└── lib/
    ├── matching_tools.py            # Perfect matching enumeration for K_{2n}
    ├── gram_tools.py                # Gram matrix construction with phases
    ├── k4_tools.py                  # K₄-specific: Pfaffian, Hessian, Fisher
    ├── k6_tools.py                  # K₆-specific: verified from project
    ├── k8_tools.py                  # K₈-specific: Heawood, Z₇, Yukawa
    └── johnson_tools.py             # Johnson scheme: closed-form eigenvalues
```

---

## III. What Each Test Verifies

### Foundational (Tests 1–2)

**Test 01 — Saturation Equation**
- Claim: (2n−1)!! = n(2n−1) has exactly two solutions: n=1 (trivial) and n=3 (K₆)
- Method: Direct evaluation for n=1..100, symbolic proof that f(n) = (2n-1)!! - n(2n-1) is strictly increasing for n≥4

**Test 02 — K₄ Matching Algebra**
- Claim: V₄ = Z₂×Z₂, abelian, six-face equivalence
- Method: Enumerate 3 matchings of K₄, compute products, verify group table
- Claim: Order-5 obstruction (no abelian group of order 5 with all involutions)
- Method: Exhaustive check of Z₅

### K₄ Sector (Tests 3–7)

**Test 03 — Pfaffian Classification**
- Claim: 15 four-edge subgraphs, classified as 12 hub-spoke + 1 sequential HC + 2 scrambled HC
- Claim: Pf(D) ≠ 0 for exactly 13/15
- Method: Enumerate all C(6,4)=15 subgraphs, compute Pf for each

**Test 04 — D² Eigenvalues**
- Claim: Hub-spoke eigenvalues s²(2±√3), ratio R = 7+4√3
- Claim: A₀* = (2−√3)² ≈ 0.0718
- Claim: s_crit = √(2ln(2+√3)/√3) ≈ 1.233
- Method: Direct eigenvalue computation for all 15 subgraphs

**Test 05 — Hessian Signature**
- Claim: Hessian of I[w] has signature (1,3) for 13/15 subgraphs at s > s_crit
- Method: Compute Hessian numerically, count positive/negative eigenvalues

**Test 06 — Heat Kernel Coefficients**
- Claim: a₀ = Tr(1) = 20, b₀ = Tr(γ₅) = −4, b₀/a₀ = −1/5
- Method: Direct trace computation with K₄ internal spectral triple

**Test 07 — Arrow of Time**
- Claim: [D, D*] ≠ 0 (informational asymmetry)
- Claim: Tr(C) = 0 for commutator C
- Method: Construct D with Z₃ flux, compute commutator

### K₆ Sector (Tests 8–11)

**Test 08 — K₆ Matchings**
- Claim: 15 = 5!! perfect matchings of K₆
- Claim: Overlap distribution 120:90:15 (pairs sharing 0,1,2 edges → overlaps 0,2,4)
- Method: Enumerate all matchings, compute pairwise overlaps

**Test 09 — K₆ Gram Matrix**
- Claim: With sorted assignment, a₂ = 3.3060051829, a₄ = 4.0681621736, R = 0.3722127085
- Claim: Full 15-eigenvalue spectrum matches published values
- Method: Construct Gram matrix with Z₃ phases, diagonalize

**Test 10 — Higgs Mass**
- Claim: mH² = 2(a₄/a₂)mW² = 126.1 GeV (c = 4/a₂) or mH² = 8(R/c)mW²
- Claim: Over 50,000 random assignments, mean mH ≈ 122 GeV, std ≈ 10
- Claim: R ∈ [0.17, 0.50] for all assignments
- Method: Monte Carlo over random direction/phase assignments

**Test 11 — BZ Independence**
- Claim: At sorted vacuum, Var(S₂(k)) = 0 (k-independent)
- Claim: Gauge kinetic coefficient = a₂
- Method: Sample BZ zone, verify variance

### Birefringence (Test 12)

**Test 12 — CMB Predictions**
- Claim: α₀_max = 0.41° from (2−√3)²
- Claim: Five falsifiable predictions (stated, not derived here)
- Method: Compute from K₄ eigenvalue ratio and heat kernel coefficients

### Johnson Scheme / Polynomial Reduction (Tests 13–15)

**Test 13 — Johnson Scheme Eigenvalues**
- Claim: O has exactly 3 distinct eigenvalues for all K_{2n}
- Claim: λ_max = 2n(2n−3)!!, λ_mid = 4(n−1)(2n−5)!!, d_phys = n(2n−3)
- Method: Construct overlap matrices for K₆, K₈, K₁₀, K₁₂; verify eigenvalues against formulas

**Test 14 — Budget Identity**
- Claim: λ_mid · d_phys = 2(n−1) · λ_max for all n
- Claim: Physical weight fraction = (2n−2)/(2n−1)
- Method: Verify algebraically and computationally for n=3..8

**Test 15 — Physical Projector**
- Claim: Π_phys = O(O − λ_max I) / [λ_mid(λ_mid − λ_max)] is exact projector
- Claim: Hub memory in physical sector = (2n−2)/(2n−1)
- Method: Construct projector, verify Π² = Π, verify hub projection fraction

### K₈ Sector (Tests 16–18)

**Test 16 — K₈ Structure**
- Claim: 105 = 7!! perfect matchings
- Claim: Genus γ(K₈) = ⌈(5)(4)/12⌉ = 2
- Claim: 28 edges decompose into 4 classes of 7 via Heawood embedding
- Method: Enumerate matchings, verify Ringel-Youngs formula

**Test 17 — K₈ Vacuum**
- Claim: 105×105 Gram matrix has 6D vacuum eigenspace at λ_vac ≈ 1.9595
- Claim: Vacuum decomposes as ρ₁ ⊕ ρ₂ ⊕ ρ₃ under Z₇ (no trivial component)
- Method: Construct Gram matrix, diagonalize, verify Z₇ sector decomposition

**Test 18 — Yukawa Hierarchy**
- Claim: K₆ vacuum projects purely into ρ₁
- Claim: ρ₁ Yukawa hierarchy = 415:135:1
- Claim: No K₆ eigenvector reaches ρ₂ (selection rule)
- Method: Project K₆ vacuum into K₈ space, compute Yukawa eigenvalues

### Higher Levels (Tests 19–20)

**Test 19 — K₁₀/K₁₂ Computation**
- Claim: K₁₀ vacuum has degree-6 minimal polynomial with Gal = C₂≀C₃
- Claim: Johnson scheme predictions verified at K₁₀, K₁₂
- Method: Construct Gram matrices, verify eigenvalues, compute minimal polynomials

**Test 20 — Galois Disjointness**
- Claim: ℚ(λ₆) ∩ ℚ(λ₈) = ℚ
- Claim: ℚ(λ₆) = ℚ(√5), degree 2
- Method: Compute minimal polynomials, verify field structure (requires sympy)

### Structural (Tests 21–22)

**Test 21 — d=4 Uniqueness**
- Claim: Tr(D⁴) = ½[Tr(D²)]² − 4·Σ Pf(D_S)² with C(N,4) sub-Pfaffians
- Claim: Unique irreducible Pfaffian invariant iff N=4
- Method: Verify Newton identity for K₃, K₄, K₅, K₆

**Test 22 — Bug Demonstration**
- Claim: Original mH = 125 GeV required BOTH .real bug AND unsorted enumeration
- Claim: Fix either bug alone → mH shifts 7–16 GeV
- Method: Four-way computation (sorted/unsorted × .real/Hermitian)

---

## IV. Documentation Deliverables

### For Each Paper

1. **LaTeX source** (already drafted)
2. **Verification script** mapping every theorem/proposition to a computational check
3. **Cross-reference table**: Theorem number → test file → specific assertion

### Master Verification Report

`run_all.py` produces a single summary:

```
=== INFORMATIONAL ASYMMETRY PROGRAM: VERIFICATION REPORT ===
Date: 2026-02-19
Python: 3.x.x, NumPy: x.x.x, SciPy: x.x.x

Test 01: Saturation equation ........................ PASS (2 claims, 2 verified)
Test 02: K₄ matching algebra ....................... PASS (4 claims, 4 verified)
...
Test 22: Bug demonstration ......................... PASS (3 claims, 3 verified)

TOTAL: XX/XX claims verified. 0 failures.
```

---

## V. Principles

1. **Every numerical claim has a test.** No exceptions.
2. **Tests use only standard libraries** (numpy, scipy, sympy). No custom physics packages.
3. **Tests are self-contained.** Each can run independently.
4. **Tests report what they check.** Each prints the claim, the computed value, the expected value, and PASS/FAIL.
5. **Tolerances are explicit.** Exact integer results use ==. Floating point uses stated tolerance.
6. **The .real bug is preserved as a test**, not hidden. Demonstrating you found and fixed your own bugs is stronger than pretending they didn't exist.
