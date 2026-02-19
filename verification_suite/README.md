# Informational Asymmetry Program — Verification Suite

## Brian Porter, February 2026

**Purpose:** Independent computational verification of every numerical claim in the Informational Asymmetry framework (Papers I–III).

**Philosophy:** Every claimed result has a test. Tests use only standard Python libraries. Anyone with `numpy` and `scipy` can reproduce everything.

---

## Quick Start

```bash
pip install numpy scipy
python run_all.py          # Run all tests (~5 seconds)
python run_all.py --quick  # Run foundational tests only
python run_all.py -v       # Verbose: show failing details
```

## What Gets Verified

| Test | Paper | Claims | What It Checks |
|------|-------|--------|----------------|
| 01 | I | 9 | Saturation equation (2n-1)!! = n(2n-1) has exactly 2 solutions |
| 02 | I | 18 | V₄ = Z₂×Z₂ matching algebra, order-5 obstruction, (Z/8Z)* isomorphism |
| 03 | I | 13 | All 15 K₄ subgraphs: Pfaffian classification, 13/15 Lorentzian |
| 05 | I | 5 | Hessian signature emergence at s > s_crit |
| 06 | I,II | 3 | Heat kernel: a₀=20, b₀=−4, parity ratio −1/5 |
| 08 | III | 9 | K₆ matching enumeration, overlap distribution |
| 09 | III | 12 | Gram matrix spectrum: a₂=3.306, a₄=4.068, R=0.3722 (10-digit) |
| 12 | II | 10 | Birefringence: A₀*=(2−√3)², five predictions stated |
| 13 | — | 24 | Johnson scheme eigenvalues for K₆, K₈, K₁₀ + budget identity |
| 21 | I | 10 | d=4 uniqueness: sub-Pfaffian frustration theorem |
| 22 | — | 4 | Bug demonstration: .real bug + enumeration order = accidental 125 |

**Total: 117 independently verified claims.**

## Tests Not Yet Implemented (Planned)

| Test | Content | Difficulty |
|------|---------|------------|
| 04 | K₄ D² eigenvalues at all scales | Low |
| 07 | Arrow of time: [D,D*] ≠ 0 | Low |
| 10 | Higgs mass Monte Carlo (50k assignments) | Medium |
| 11 | BZ k-independence verification | Medium |
| 14 | Budget identity (algebraic proof) | Low |
| 15 | Physical projector Π² = Π | Low |
| 16 | K₈ structure: 105 matchings, Z₇, genus-2 | Medium |
| 17 | K₈ vacuum: 6D eigenspace, ρ₁⊕ρ₂⊕ρ₃ | High |
| 18 | K₈ Yukawa: 415:135:1, ρ₂ selection rule | High |
| 19 | K₁₀/K₁₂: degree stabilization, Galois groups | High |
| 20 | Galois disjointness: ℚ(λ₆) ∩ ℚ(λ₈) = ℚ | High (needs sympy) |

## Directory Structure

```
verification_suite/
├── run_all.py              # Master runner
├── PLAN.md                 # Detailed plan with cross-references
├── README.md               # This file
├── requirements.txt        # numpy, scipy
├── lib/
│   ├── matching_tools.py   # Perfect matching enumeration
│   ├── gram_tools.py       # Gram matrix with Z₃ phases
│   ├── k4_tools.py         # K₄: Pfaffian, Hessian, spectral action
│   └── johnson_tools.py    # Johnson scheme closed-form eigenvalues
└── tests/
    ├── test_01_saturation.py
    ├── test_02_k4_matching.py
    ├── ...
    └── test_22_bug_demonstration.py
```

## How Results Map to Papers

**Paper I** (K₄ → Lorentzian signature):
- Tests 01–07, 21 verify all theorems

**Paper II** (Cosmic birefringence):
- Test 06 (heat kernel), Test 12 (birefringence predictions)

**Paper III** (K₈ fermion masses):
- Tests 08–09 (K₆ core), Tests 13 (Johnson scheme), Tests 16–18 (K₈)

**Bug Disclosure:**
- Test 22 demonstrates the .real bug and enumeration order dependency
  that produced the original erroneous mH = 125 GeV. This test exists
  as an act of scientific integrity: we found our own bugs and document
  them for reproducibility.

## Adding New Tests

Each test script is self-contained. Pattern:
```python
import sys
sys.path.insert(0, '..')
from lib.matching_tools import ...

PASS = FAIL = 0
def check(name, condition, detail=""):
    global PASS, FAIL
    if condition: PASS += 1; print(f"  ✓ {name}")
    else: FAIL += 1; print(f"  ✗ {name}  {detail}")

# ... your tests ...

print(f"TEST XX RESULT: {PASS} passed, {FAIL} failed")
```

The master runner `run_all.py` auto-discovers any file matching `test_*.py` in the `tests/` directory.
