"""
Test 06: K₄ Heat Kernel Coefficients
======================================
Paper I / Paper II:
  a₀ = Tr(1) = 20
  b₀ = Tr(γ₅) = -4
  b₀/a₀ = -1/5 (parity ratio, purely combinatorial)
"""
import sys
sys.path.insert(0, '..')
from lib.k4_tools import heat_kernel_k4_internal

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
print("TEST 06: K₄ Heat Kernel Coefficients")
print("=" * 70)

a0, b0 = heat_kernel_k4_internal()

print(f"\nK₄ internal spectral triple:")
print(f"  Eigenvalues λ_n = n, degeneracy D(n) = 2n, chirality (-1)^(n+1)")
print(f"  n=1: D=2, γ₅=+1  → contributes +2 to a₀, +2 to b₀")
print(f"  n=2: D=4, γ₅=-1  → contributes +4 to a₀, -4 to b₀")
print(f"  n=3: D=6, γ₅=+1  → contributes +6 to a₀, +6 to b₀")
print(f"  n=4: D=8, γ₅=-1  → contributes +8 to a₀, -8 to b₀")

print(f"\n  a₀ = 2+4+6+8 = {a0}")
print(f"  b₀ = 2-4+6-8 = {b0}")
print(f"  b₀/a₀ = {b0}/{a0} = {b0/a0}")

check("a₀ = 20", a0 == 20)
check("b₀ = -4", b0 == -4)
check("b₀/a₀ = -1/5", abs(b0/a0 - (-1/5)) < 1e-15)

# This ratio enters the birefringence prediction
print(f"\n  The parity ratio -1/5 is a combinatorial invariant of K₄.")
print(f"  It determines the coupling strength of parity-violating torsion.")

print(f"\n{'='*70}")
print(f"TEST 06 RESULT: {PASS} passed, {FAIL} failed")
print(f"{'='*70}")
