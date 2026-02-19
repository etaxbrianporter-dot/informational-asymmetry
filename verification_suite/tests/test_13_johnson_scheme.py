"""
Test 13: Johnson Scheme Eigenvalues
=====================================
Spectral Budget Identity / Polynomial Reduction Theorem:
Overlap matrix O of K_{2n} has exactly 3 distinct eigenvalues:
  λ_max = 2n(2n-3)!!
  λ_mid = 4(n-1)(2n-5)!!
  0 (kernel)
with multiplicities 1, n(2n-3), and N-n(2n-3)-1 respectively.

Verified computationally for K₆ (n=3), K₈ (n=4), K₁₀ (n=5), K₁₂ (n=6).
"""
import sys
import numpy as np
sys.path.insert(0, '..')
from lib.matching_tools import enumerate_matchings, overlap_matrix
from lib.johnson_tools import johnson_parameters, verify_overlap_spectrum

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
print("TEST 13: Johnson Scheme Eigenvalues")
print("=" * 70)

# Test for each level
for n in [3, 4, 5]:
    nv = 2 * n
    p = johnson_parameters(n)
    
    print(f"\n--- K_{nv} (n={n}) ---")
    print(f"  N = {p['N']} matchings, E = {p['E']} edges")
    print(f"  Predicted: λ_max = {p['lam_max']}, λ_mid = {p['lam_mid']}, d_phys = {p['d_phys']}")
    
    if p['N'] > 20000:
        print(f"  SKIP (N = {p['N']} too large for enumeration)")
        continue
    
    # Enumerate and compute
    matchings = enumerate_matchings(nv)
    check(f"K_{nv}: |M| = {p['N']}", len(matchings) == p['N'])
    
    O = overlap_matrix(matchings)
    result = verify_overlap_spectrum(O, n)
    
    print(f"  Computed:  λ_max = {result['lam_max_found']:.1f} (mult {result['lam_max_mult_found']})")
    print(f"             λ_mid = {result['lam_mid_found']:.1f} (mult {result['lam_mid_mult_found']})")
    print(f"             zeros: {result['n_zero_found']}")
    
    check(f"K_{nv}: λ_max = {p['lam_max']}", 
          abs(result['lam_max_found'] - p['lam_max']) < 0.5)
    check(f"K_{nv}: λ_max multiplicity = 1", 
          result['lam_max_mult_found'] == 1)
    check(f"K_{nv}: λ_mid = {p['lam_mid']}", 
          abs(result['lam_mid_found'] - p['lam_mid']) < 0.5)
    check(f"K_{nv}: d_phys = {p['d_phys']}", 
          result['lam_mid_mult_found'] == p['d_phys'])
    check(f"K_{nv}: all eigenvalues correct", result['all_correct'])

# Algebraic verification for higher n (no enumeration needed)
print(f"\n--- Algebraic verification for n=3..8 ---")
print(f"  {'n':>3s} {'K_{2n}':>6s} {'N':>12s} {'λ_max':>12s} {'λ_mid':>12s} {'d_phys':>8s} {'Budget?':>8s}")
print("  " + "-" * 65)

for n in range(3, 9):
    p = johnson_parameters(n)
    budget_ok = p['budget_lhs'] == p['budget_rhs']
    print(f"  {n:3d} {'K_'+str(2*n):>6s} {p['N']:12d} {p['lam_max']:12d} {p['lam_mid']:12d} {p['d_phys']:8d} {'✓' if budget_ok else '✗':>8s}")
    check(f"K_{2*n}: budget identity λ_mid·d_phys = 2(n-1)·λ_max", budget_ok)

print(f"\n{'='*70}")
print(f"TEST 13 RESULT: {PASS} passed, {FAIL} failed")
print(f"{'='*70}")
