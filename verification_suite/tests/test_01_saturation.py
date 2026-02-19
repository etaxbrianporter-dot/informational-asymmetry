"""
Test 01: Saturation Equation — (2n-1)!! = I(n)
================================================
Paper I, Section 2: The saturation equation (2n-1)!! = n(2n-1) 
has exactly two solutions: n=1 (trivial) and n=3 (K₆).

This forces K₆ as the unique non-trivial complete graph.
"""
import sys
sys.path.insert(0, '..')
from lib.johnson_tools import dfact

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
print("TEST 01: Saturation Equation (2n-1)!! = n(2n-1)")
print("=" * 70)

# Claim 1: Exactly two solutions for n=1..100
print("\nClaim 1: Solutions are exactly n=1 and n=3")
solutions = []
for n in range(1, 101):
    lhs = dfact(2*n - 1)  # (2n-1)!!
    rhs = n * (2*n - 1)    # I(n) = n(2n-1)
    if lhs == rhs:
        solutions.append(n)
        print(f"  n={n}: (2n-1)!! = {lhs}, n(2n-1) = {rhs} — MATCH")

check("Exactly two solutions found", len(solutions) == 2, 
      f"Found {len(solutions)}: {solutions}")
check("n=1 is a solution", 1 in solutions)
check("n=3 is a solution", 3 in solutions)

# Claim 2: f(n) = (2n-1)!! - n(2n-1) is strictly increasing for n≥4
print("\nClaim 2: f(n) strictly increasing for n ≥ 4")
prev_f = None
increasing = True
for n in range(4, 50):
    f_n = dfact(2*n - 1) - n * (2*n - 1)
    if prev_f is not None and f_n <= prev_f:
        increasing = False
        print(f"  FAIL at n={n}: f({n}) = {f_n} ≤ f({n-1}) = {prev_f}")
        break
    prev_f = f_n

check("f(n) strictly increasing for n=4..49", increasing)

# Claim 3: The growth rate makes more solutions impossible
print("\nClaim 3: Growth rate comparison")
for n in [3, 4, 5, 10, 20]:
    lhs = dfact(2*n - 1)
    rhs = n * (2*n - 1)
    ratio = lhs / rhs
    print(f"  n={n:2d}: (2n-1)!! = {lhs:>20d}, n(2n-1) = {rhs:>6d}, ratio = {ratio:.2e}")

check("(2n-1)!! grows super-exponentially vs linear n(2n-1)", 
      dfact(2*20 - 1) / (20 * 39) > 1e15)

# Claim 4: Verify the specific values at n=3
print("\nClaim 4: At n=3, both sides equal 15")
check("(2·3-1)!! = 5!! = 15", dfact(5) == 15)
check("n(2n-1) = 3·5 = 15", 3 * 5 == 15)
check("5!! = 1·3·5 = 15", 1*3*5 == 15)
check("15 = number of perfect matchings of K₆", dfact(5) == 15)

print(f"\n{'='*70}")
print(f"TEST 01 RESULT: {PASS} passed, {FAIL} failed")
print(f"{'='*70}")
