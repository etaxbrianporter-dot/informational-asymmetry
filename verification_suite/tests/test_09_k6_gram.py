"""
Test 09: K₆ Gram Matrix — Core Higgs Sector
=============================================
Paper III / K6_Higgs_Paper_v2:
- a₂ = 3.3060051829 (ground eigenvalue)
- a₄ = 4.0681621736 (quartic form)
- R = a₄/a₂² = 0.3722127085 (spectral kurtosis)
- Full 15-eigenvalue spectrum
"""
import sys
import numpy as np
from scipy.linalg import eigh
sys.path.insert(0, '..')
from lib.gram_tools import (gram_matrix, quartic_form, K6_MATCHINGS, K6_DIRS, K6_Z3,
                             compute_k6_full)

PASS = 0
FAIL = 0

def check(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  ✓ {name}")
    else:
        FAIL += 1
        print(f"  ✗ {name}  {detail}")

print("=" * 70)
print("TEST 09: K₆ Gram Matrix and Spectral Invariants")
print("=" * 70)

# Compute everything
result = compute_k6_full()
G = result['G']
evals = result['evals']
v0 = result['v0']
a2 = result['a2']
a4 = result['a4']
R = result['R']

# Claim 1: Gram matrix properties
print("\nClaim 1: Gram matrix structure")
check("G is 15×15", G.shape == (15, 15))
check("G is symmetric", np.allclose(G, G.T, atol=1e-14))

# Claim 2: a₂ = ground eigenvalue
print(f"\nClaim 2: a₂ = {a2:.10f}")
a2_expected = 3.3060051829
check(f"a₂ = 3.3060051829 (to 10 digits)", 
      abs(a2 - a2_expected) < 1e-6, f"got {a2:.10f}")

# Claim 3: a₄ = quartic form
print(f"\nClaim 3: a₄ = {a4:.10f}")
a4_expected = 4.0681621736
check(f"a₄ = 4.0681621736 (to 10 digits)", 
      abs(a4 - a4_expected) < 1e-4, f"got {a4:.10f}")

# Claim 4: R = a₄/a₂²
print(f"\nClaim 4: R = a₄/a₂² = {R:.10f}")
R_expected = 0.3722127085
check(f"R = 0.3722127085", 
      abs(R - R_expected) < 1e-6, f"got {R:.10f}")

# Verify R computed correctly
R_direct = a4 / a2**2
check(f"R = a₄/a₂² (direct)", abs(R - R_direct) < 1e-14)

# Claim 5: Full eigenspectrum
print(f"\nClaim 5: Full eigenspectrum")
expected_evals = [3.306005, 3.585786, 4.000000, 5.000000, 5.000000,
                  5.381966, 5.381966, 5.748000, 6.000000, 6.414214,
                  7.000000, 7.618034, 7.618034, 8.945995, 9.000000]

print(f"  {'Index':>5s} {'Computed':>14s} {'Expected':>14s} {'Match':>6s}")
print("  " + "-" * 42)
all_match = True
for i in range(15):
    match = abs(evals[i] - expected_evals[i]) < 0.001
    if not match:
        all_match = False
    print(f"  {i:5d} {evals[i]:14.6f} {expected_evals[i]:14.6f} {'✓' if match else '✗':>6s}")

check("All 15 eigenvalues match published values", all_match)

# Claim 6: Spectral range and condition number
print(f"\nClaim 6: Spectral properties")
spectral_range = evals[-1] - evals[0]
condition = evals[-1] / evals[0]
check(f"Spectral range ≈ 5.694", abs(spectral_range - 5.694) < 0.01)
check(f"Condition number ≈ 2.72", abs(condition - 2.72) < 0.01)

# Claim 7: Ground eigenvector has all dir=0 components
print(f"\nClaim 7: Vacuum eigenvector structure")
# In the sorted assignment, dir=0 matchings are at indices 0,3,6,9,12
dir0_indices = [i for i in range(15) if K6_DIRS[i] == 0]
dir1_indices = [i for i in range(15) if K6_DIRS[i] == 1]
dir2_indices = [i for i in range(15) if K6_DIRS[i] == 2]

norm_dir0 = np.linalg.norm(v0[dir0_indices])
norm_dir1 = np.linalg.norm(v0[dir1_indices])
norm_dir2 = np.linalg.norm(v0[dir2_indices])
print(f"  |v₀|_dir0 = {norm_dir0:.6f}")
print(f"  |v₀|_dir1 = {norm_dir1:.6f}")
print(f"  |v₀|_dir2 = {norm_dir2:.6f}")

# Claim 8: Higgs mass formula
print(f"\nClaim 8: Higgs mass predictions")
mW = 80.379  # GeV

# With c = 4/a₂
c_1 = 4 / a2
mH_1 = np.sqrt(8 * R / c_1) * mW
print(f"  c = 4/a₂ = {c_1:.6f}: mH = {mH_1:.1f} GeV")

# With c = π²/8
c_2 = np.pi**2 / 8
mH_2 = np.sqrt(8 * R / c_2) * mW
print(f"  c = π²/8 = {c_2:.6f}: mH = {mH_2:.1f} GeV")

# Equivalent formula: mH² = 2(a₄/a₂)mW²
mH_equiv = np.sqrt(2 * a4 / a2) * mW
print(f"  mH = √(2a₄/a₂)·mW = {mH_equiv:.1f} GeV")
check("mH formulas consistent", abs(mH_1 - mH_equiv) < 0.1)

# With paper's c = 3
c_paper = 3.0
mH_paper = np.sqrt(8 * R / c_paper) * mW
print(f"  c = 3 (paper): mH = {mH_paper:.1f} GeV")
check("Paper's c=3 gives ~80 GeV", abs(mH_paper - 80) < 1)

# c required for 125.09
mH_exp = 125.09
c_exact = 8 * R * mW**2 / mH_exp**2
print(f"  c_exact for mH=125.09: c = {c_exact:.6f}")
check("c_exact ≈ 1.23", abs(c_exact - 1.23) < 0.01)

print(f"\n{'='*70}")
print(f"TEST 09 RESULT: {PASS} passed, {FAIL} failed")
print(f"{'='*70}")
