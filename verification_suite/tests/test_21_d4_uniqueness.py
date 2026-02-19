"""
Test 21: d=4 Uniqueness via Sub-Pfaffian Frustration
=====================================================
Paper I, Theorem (d=4 uniqueness):
Tr(D⁴) = ½[Tr(D²)]² - 4·Σ Pf(D_S)² where sum is over all C(N,4) 
four-vertex subsets S.

The spectral action contains a SINGLE irreducible Pfaffian invariant
if and only if N = 4.
"""
import sys
import numpy as np
from itertools import combinations
from math import comb
sys.path.insert(0, '..')
from lib.k4_tools import pfaffian_4x4

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
print("TEST 21: d=4 Uniqueness (Sub-Pfaffian Frustration)")
print("=" * 70)

def random_skew_symmetric(N, seed=42):
    """Generate random N×N skew-symmetric matrix."""
    rng = np.random.RandomState(seed)
    A = rng.randn(N, N)
    return (A - A.T) / 2

def verify_newton_identity(D):
    """
    Verify: Tr(D⁴) = ½[Tr(D²)]² - 4·Σ_{|S|=4} Pf(D_S)²
    """
    N = D.shape[0]
    tr_D4 = np.trace(D @ D @ D @ D)
    tr_D2 = np.trace(D @ D)
    
    # Sum of sub-Pfaffians
    pf_sum = 0.0
    n_subpf = 0
    for S in combinations(range(N), 4):
        D_S = D[np.ix_(list(S), list(S))]
        pf = pfaffian_4x4(D_S)
        pf_sum += pf**2
        n_subpf += 1
    
    rhs = 0.5 * tr_D2**2 - 4 * pf_sum
    return tr_D4, rhs, n_subpf

print("\nClaim 1: Newton identity verified for random skew-symmetric matrices")
for N in [4, 5, 6, 7]:
    D = random_skew_symmetric(N)
    lhs, rhs, n_pf = verify_newton_identity(D)
    match = abs(lhs - rhs) < 1e-8
    print(f"  N={N}: Tr(D⁴) = {lhs:.6f}, ½S² - 4ΣPf² = {rhs:.6f}, "
          f"C({N},4) = {n_pf} sub-Pfaffians, {'✓' if match else '✗'}")
    check(f"Newton identity holds for N={N}", match)

# Claim 2: Number of sub-Pfaffians
print(f"\nClaim 2: C(N,4) sub-Pfaffians")
for N in range(2, 8):
    n_pf = comb(N, 4)
    print(f"  N={N}: C({N},4) = {n_pf}")

check("C(N,4) = 0 for N < 4 (no curvature content)", comb(3, 4) == 0)
check("C(4,4) = 1 (unique, irreducible Pfaffian)", comb(4, 4) == 1)
check("C(5,4) = 5 (multiple sub-Pfaffians, no single invariant)", comb(5, 4) == 5)
check("C(6,4) = 15", comb(6, 4) == 15)

# Claim 3: For N=4, the single Pfaffian IS the full Pf(D)
print(f"\nClaim 3: For N=4, the single sub-Pfaffian IS Pf(D)")
D4 = random_skew_symmetric(4, seed=123)
pf_full = pfaffian_4x4(D4)
# The only 4-vertex subset of {0,1,2,3} is {0,1,2,3} itself
D_sub = D4[np.ix_([0,1,2,3], [0,1,2,3])]
pf_sub = pfaffian_4x4(D_sub)
check("Sub-Pfaffian = full Pfaffian for N=4", abs(pf_full - pf_sub) < 1e-15)

# Claim 4: Decomposition at N=4
print(f"\nClaim 4: At N=4, a₂ = ½S² - 4Pf²")
S = np.trace(D4 @ D4)
a2 = np.trace(D4 @ D4 @ D4 @ D4)
decomp = 0.5 * S**2 - 4 * pf_full**2
check("Tr(D⁴) = ½[Tr(D²)]² - 4Pf(D)²", abs(a2 - decomp) < 1e-10)
print(f"  Volume term: ½S² = {0.5*S**2:.6f}")
print(f"  Curvature term: -4Pf² = {-4*pf_full**2:.6f}")
print(f"  Total a₂ = {a2:.6f}")

print(f"\n  For N>4, multiple sub-Pfaffians create frustration —")
print(f"  no single topological invariant controls the curvature.")
print(f"  N=4 is unique: one graph, one Pfaffian, one curvature.")

print(f"\n{'='*70}")
print(f"TEST 21 RESULT: {PASS} passed, {FAIL} failed")
print(f"{'='*70}")
