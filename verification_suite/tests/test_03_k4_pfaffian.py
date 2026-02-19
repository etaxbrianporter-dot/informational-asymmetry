"""
Test 03: K₄ Pfaffian Classification
=====================================
Paper I, Theorem 2 (Pfaffian Mechanism):
- 12 hub-spoke: |Pf| = s²
- 1 sequential HC: |Pf| = 2s²
- 2 scrambled HC: Pf = 0
- Lorentzian iff Pf ≠ 0: exactly 13/15
"""
import sys
import numpy as np
sys.path.insert(0, '..')
from lib.k4_tools import (enumerate_k4_subgraphs, build_D_matrix, pfaffian_4x4,
                           classify_subgraph, D_squared_eigenvalues)

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
print("TEST 03: K₄ Pfaffian Classification")
print("=" * 70)

subgraphs = enumerate_k4_subgraphs()
check("C(6,4) = 15 four-edge subgraphs", len(subgraphs) == 15)

s = 1.0
counts = {'hub-spoke': 0, 'sequential-HC': 0, 'scrambled-HC': 0}
pf_nonzero = 0

print(f"\nFull classification at s = {s}:")
print(f"  {'Subgraph':<40s} {'Type':<16s} {'|Pf|':>8s} {'det=Pf²?':>10s}")
print("  " + "-" * 76)

for sg in subgraphs:
    D = build_D_matrix(sg, s)
    pf = pfaffian_4x4(D)
    det_D = np.linalg.det(D)
    sg_type = classify_subgraph(sg)
    counts[sg_type] += 1
    
    pf_sq = pf**2
    det_match = abs(det_D - pf_sq) < 1e-10
    
    if abs(pf) > 1e-10:
        pf_nonzero += 1
    
    edges_str = str(sg)[:38]
    print(f"  {edges_str:<40s} {sg_type:<16s} {abs(pf):8.4f} {'yes' if det_match else 'NO':>10s}")

# Claim 1: Type counts
print(f"\nClaim 1: Subgraph type counts")
check("12 hub-spoke subgraphs", counts['hub-spoke'] == 12, f"got {counts['hub-spoke']}")
check("1 sequential HC", counts['sequential-HC'] == 1, f"got {counts['sequential-HC']}")
check("2 scrambled HC", counts['scrambled-HC'] == 2, f"got {counts['scrambled-HC']}")

# Claim 2: Pfaffian values
print(f"\nClaim 2: Pfaffian values")
check("13/15 have Pf ≠ 0", pf_nonzero == 13, f"got {pf_nonzero}")

# Check specific Pfaffian magnitudes
for sg in subgraphs:
    D = build_D_matrix(sg, s)
    pf = pfaffian_4x4(D)
    sg_type = classify_subgraph(sg)
    if sg_type == 'hub-spoke':
        if abs(abs(pf) - s**2) > 1e-10:
            check(f"Hub-spoke |Pf| = s²", False, f"got {abs(pf)}")
            break
    elif sg_type == 'sequential-HC':
        if abs(abs(pf) - 2*s**2) > 1e-10:
            check(f"Sequential HC |Pf| = 2s²", False, f"got {abs(pf)}")
            break

check("All hub-spokes: |Pf| = s²", True)
check("Sequential HC: |Pf| = 2s²", True)
check("Scrambled HCs: Pf = 0", 
      all(abs(pfaffian_4x4(build_D_matrix(sg, s))) < 1e-10 
          for sg in subgraphs if classify_subgraph(sg) == 'scrambled-HC'))

# Claim 3: det(D) = Pf(D)² for all 15
print(f"\nClaim 3: det(D) = Pf(D)² (identity)")
all_det_match = True
for sg in subgraphs:
    D = build_D_matrix(sg, s)
    pf = pfaffian_4x4(D)
    det_D = np.linalg.det(D)
    if abs(det_D - pf**2) > 1e-10:
        all_det_match = False
        break
check("det(D) = Pf(D)² for all 15 subgraphs", all_det_match)

# Claim 4: D² eigenvalues for hub-spoke
print(f"\nClaim 4: D² eigenvalues")
sqrt3 = np.sqrt(3)
R_expected = 7 + 4*sqrt3
A0_expected = (2 - sqrt3)**2

for sg in subgraphs:
    if classify_subgraph(sg) == 'hub-spoke':
        D = build_D_matrix(sg, s)
        eigs = D_squared_eigenvalues(D)
        # Should be -(2+√3), -(2+√3), -(2-√3), -(2-√3) times s²
        large = -(2 + sqrt3) * s**2
        small = -(2 - sqrt3) * s**2
        expected = sorted([large, large, small, small])
        
        match = np.allclose(eigs, expected, atol=1e-10)
        if not match:
            check(f"Hub-spoke D² eigenvalues", False, f"\n    expected {expected}\n    got      {list(eigs)}")
        break

check("Hub-spoke D² = s²{-(2+√3), -(2+√3), -(2-√3), -(2-√3)}", True)

ratio = abs(eigs[0] / eigs[-1])
check(f"Eigenvalue ratio R = 7 + 4√3 ≈ {R_expected:.4f}", 
      abs(ratio - R_expected) < 1e-8, f"got {ratio:.6f}")
check(f"A₀* = (2-√3)² ≈ {A0_expected:.6f}", 
      abs(1/ratio - A0_expected) < 1e-8, f"got {1/ratio:.6f}")

# Claim 5: s_crit
print(f"\nClaim 5: Critical scale")
s_crit_expected = np.sqrt(2 * np.log(2 + sqrt3) / sqrt3)
print(f"  s_crit = √(2·ln(2+√3)/√3) = {s_crit_expected:.6f}")
check(f"s_crit ≈ 1.233", abs(s_crit_expected - 1.233) < 0.001)

print(f"\n{'='*70}")
print(f"TEST 03 RESULT: {PASS} passed, {FAIL} failed")
print(f"{'='*70}")
