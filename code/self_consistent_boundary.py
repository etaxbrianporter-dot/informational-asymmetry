#!/usr/bin/env python3
"""
Self-Consistent Boundary Problem for c₁ = 8
=============================================

The spectral action at cutoff Λ gives boundary conditions:
  1/αᵢ(Λ) = (f₂Λ²/2π²) × cᵢ  + (f₀/4π²) × aᵢ × yₜ²(Λ) + ...

where:
  cᵢ = {4, 6, 8}  (spectator mechanism: dim K₂ᵢ)
  aᵢ = {17/10, 3/2, 2}  (Yukawa Casimir sums)

Two parameters: K = f₂Λ²/2π²  and  ρ = f₀/(f₂Λ²)

Boundary conditions become:
  1/αᵢ(Λ) = K × cᵢ × (1 + ρ × (aᵢ/cᵢ) × yₜ²(Λ)/(2π²))

Key insight: the ratios aᵢ/cᵢ differ:
  a₁/c₁ = 17/80 = 0.2125
  a₂/c₂ = 3/12  = 0.2500  
  a₃/c₃ = 2/4   = 0.5000

So the f₀ correction breaks the simple c₁:c₂:c₃ = 8:6:4 ratios.
The question: does this HELP or HURT, and can we solve self-consistently?
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve, minimize

# ============================================================
# Physical constants
# ============================================================
M_Z = 91.1876       # Z mass in GeV
M_t = 173.1         # top pole mass
M_W = 80.379
M_H = 125.09
v   = 246.22        # Higgs VEV

# Observed couplings at M_Z (MS-bar)
alpha1_obs = 0.016943   # U(1)_Y with GUT normalization 5/3
alpha2_obs = 0.033801   # SU(2)
alpha3_obs = 0.118000   # SU(3)
sin2_obs   = 0.23122    # sin²θ_W

# SM 1-loop beta coefficients: dg_i/dt = b_i g_i^3 / (16π²)
# with t = ln(μ/M_Z)
b1 = 41/10   # U(1)
b2 = -19/6   # SU(2)  
b3 = -7      # SU(3)

# SM 2-loop beta coefficients: dg_i/dt = [b_i g_i^3 + Σ_j B_ij g_i^3 g_j^2 / (16π²)] / (16π²)
# Plus top Yukawa corrections
B = np.array([
    [199/50,  27/10, 44/5],
    [ 9/10,   35/6,  12  ],
    [11/10,    9/2, -26  ]
])

# Top Yukawa 2-loop contributions to gauge betas
# d(1/α_i)/dt gets +a_i y_t^2/(32π³) correction
a_yuk = np.array([17/10, 3/2, 2])  # Yukawa coefficients

# ============================================================
# RG evolution (2-loop with top Yukawa)
# ============================================================
def rg_system(t, y):
    """
    y = [1/α₁, 1/α₂, 1/α₃, y_t²]
    t = ln(μ/M_Z)
    """
    inv_a1, inv_a2, inv_a3, yt2 = y
    a1, a2, a3 = 1/inv_a1, 1/inv_a2, 1/inv_a3
    g = np.array([a1, a2, a3])
    
    # 1-loop
    b_vec = np.array([b1, b2, b3])
    
    # d(1/αᵢ)/dt = -bᵢ/(2π) - Σⱼ Bᵢⱼ αⱼ/(8π³) - aᵢ yₜ²/(32π³)
    d_inv = np.zeros(3)
    for i in range(3):
        d_inv[i] = -b_vec[i]/(2*np.pi)
        for j in range(3):
            d_inv[i] -= B[i,j] * g[j] / (8*np.pi**3)
        d_inv[i] -= a_yuk[i] * yt2 / (32*np.pi**3)
    
    # Top Yukawa running (1-loop dominant terms)
    # dy_t²/dt = y_t²/(16π²) × [9/2 y_t² - 8g₃² - 9/4 g₂² - 17/20 g₁²]
    dyt2 = yt2/(16*np.pi**2) * (
        9/2 * yt2 
        - 8 * a3 * 4*np.pi
        - 9/4 * a2 * 4*np.pi
        - 17/20 * a1 * 4*np.pi
    )
    # Correction: g² = 4πα
    dyt2 = yt2/(16*np.pi**2) * (
        9/2 * yt2 
        - 8 * (4*np.pi*a3)
        - 9/4 * (4*np.pi*a2)
        - 17/20 * (4*np.pi*a1)
    )
    
    return [d_inv[0], d_inv[1], d_inv[2], dyt2]

def run_up(log_Lambda_over_MZ):
    """Run couplings from M_Z up to Λ = M_Z × exp(t_max)"""
    yt2_MZ = 2 * M_t**2 / v**2   # ≈ 0.99
    
    y0 = [1/alpha1_obs, 1/alpha2_obs, 1/alpha3_obs, yt2_MZ]
    
    sol = solve_ivp(rg_system, [0, log_Lambda_over_MZ], y0,
                    method='RK45', rtol=1e-10, atol=1e-12,
                    dense_output=True)
    return sol

def get_couplings_at_Lambda(log_Lambda):
    """Get 1/αᵢ(Λ) and yₜ²(Λ) at log(Λ/M_Z)"""
    sol = run_up(log_Lambda)
    vals = sol.sol(log_Lambda)
    return vals[0], vals[1], vals[2], vals[3]

# ============================================================
# Step 1: Reproduce the 2.5% result
# ============================================================
print("=" * 72)
print("STEP 1: REPRODUCE KNOWN RESULT")
print("=" * 72)

# Scan for Λ where c₂/c₃ = 3/2
def c_ratios_at_Lambda(log_Lam):
    ia1, ia2, ia3, yt2 = get_couplings_at_Lambda(log_Lam)
    # Without f₀ correction: cᵢ ∝ 1/αᵢ(Λ)
    # c₂/c₃ = (1/α₂)/(1/α₃) should be 3/2
    return ia2/ia3

from scipy.optimize import brentq

# Find where c₂/c₃ = 3/2
log_Lam_candidates = np.linspace(np.log(1e6/M_Z), np.log(1e12/M_Z), 100)
ratios = [c_ratios_at_Lambda(ll) for ll in log_Lam_candidates]

print("\nScan: c₂/c₃ vs Λ")
for ll, r in zip(log_Lam_candidates[::10], ratios[::10]):
    Lam = M_Z * np.exp(ll)
    print(f"  Λ = {Lam:.2e} GeV:  c₂/c₃ = {r:.4f}")

# Find exact crossing
target = 1.5
try:
    log_Lam_match = brentq(lambda ll: c_ratios_at_Lambda(ll) - target, 
                           log_Lam_candidates[0], log_Lam_candidates[-1])
    Lam_match = M_Z * np.exp(log_Lam_match)
    
    ia1, ia2, ia3, yt2 = get_couplings_at_Lambda(log_Lam_match)
    c1_eff = ia1/ia3 * 4  # normalize so c₃ = 4
    c2_eff = ia2/ia3 * 4  # should be 6.00
    c3_eff = 4.0
    
    print(f"\n  Λ where c₂/c₃ = 3/2: {Lam_match:.3e} GeV")
    print(f"  c₁ = {c1_eff:.4f}  (target: 8, deviation: {(c1_eff-8)/8*100:.2f}%)")
    print(f"  c₂ = {c2_eff:.4f}  (target: 6, deviation: {(c2_eff-6)/6*100:.2f}%)")
    print(f"  c₃ = {c3_eff:.4f}")
    print(f"  yₜ²(Λ) = {yt2:.4f}")
    print(f"  sin²θ_W = {ia1/(ia1+ia2):.5f}  (obs: {sin2_obs})")
    # Note: sin²θ_W = α₁/(α₁+α₂) = (1/α₂)/(1/α₁+1/α₂) ... 
    # Actually sin²θ_W = g'²/(g²+g'²) = α₁/(α₁+α₂) with GUT normalization
    # = (5/3)α'/(α₂+(5/3)α') where α' is the non-GUT U(1) coupling
    
except Exception as e:
    print(f"  Error: {e}")
    log_Lam_match = np.log(1e8/M_Z)

# ============================================================
# Step 2: Include f₀ boundary correction
# ============================================================
print("\n" + "=" * 72)
print("STEP 2: SELF-CONSISTENT BOUNDARY WITH f₀ CORRECTION")
print("=" * 72)

print("""
Boundary conditions at Λ:
  1/αᵢ(Λ) = K × cᵢ × [1 + ρ × (aᵢ/cᵢ) × yₜ²(Λ)/(2π²)]

where K = f₂Λ²/(2π²), ρ = f₀/(f₂Λ²).

Define effective spectator dimensions:
  cᵢ_eff(ρ) = cᵢ × [1 + ρ × (aᵢ/cᵢ) × yₜ²(Λ)/(2π²)]

We need: c₁_eff/c₃_eff = 2, c₂_eff/c₃_eff = 3/2

This is an algebraic constraint on ρ given yₜ²(Λ).
""")

c_bare = np.array([8, 6, 4])
a_over_c = a_yuk / c_bare  # [0.2125, 0.25, 0.50]

def c_eff_ratio(rho, yt2_Lambda):
    """Compute c₁_eff/c₃_eff and c₂_eff/c₃_eff"""
    factor = rho * a_over_c * yt2_Lambda / (2*np.pi**2)
    c_eff = c_bare * (1 + factor)
    r13 = c_eff[0] / c_eff[2]
    r23 = c_eff[1] / c_eff[2]
    return r13, r23

# At the matching scale, yₜ²(Λ) is known from running
ia1, ia2, ia3, yt2_Lam = get_couplings_at_Lambda(log_Lam_match)

print(f"At Λ = {Lam_match:.2e} GeV:")
print(f"  yₜ²(Λ) = {yt2_Lam:.4f}")

# What ρ is needed to make c₁_eff/c₃_eff = 2 exactly?
# c₁(1 + ρ a₁ yt² / (c₁ 2π²)) / c₃(1 + ρ a₃ yt² / (c₃ 2π²)) = 2
# 8(1 + ρ × 0.2125 × yt²/(2π²)) / 4(1 + ρ × 0.50 × yt²/(2π²)) = 2
# 2(1 + ρ × 0.2125 × yt²/(2π²)) / (1 + ρ × 0.50 × yt²/(2π²)) = 2
# 2(1 + ρ × 0.2125 × yt²/(2π²)) = 2(1 + ρ × 0.50 × yt²/(2π²))
# 1 + ρ × 0.2125 × yt²/(2π²) = 1 + ρ × 0.50 × yt²/(2π²)
# ρ × (0.2125 - 0.50) × yt²/(2π²) = 0
# 
# WAIT. If c₁/c₃ = 8/4 = 2, then c₁_eff/c₃_eff = 2 requires:
# a₁/c₁ = a₃/c₃, which is 0.2125 ≠ 0.50.
# So c₁_eff/c₃_eff = 2 × (1 + ρ×0.2125×yt²/(2π²))/(1 + ρ×0.50×yt²/(2π²))
# This is ALWAYS less than 2 for ρ > 0 (since 0.2125 < 0.50).
# 
# But we NEED c₁_eff/c₃_eff < 2 to match observation!
# Currently c₁/c₃ runs to ~2.05 at the matching scale.
# So we need the f₀ correction to pull c₁_eff/c₃_eff DOWN from ~2.05.

print("\n--- Algebraic analysis ---")
print(f"  a₁/c₁ = {a_over_c[0]:.4f}")
print(f"  a₂/c₂ = {a_over_c[1]:.4f}")
print(f"  a₃/c₃ = {a_over_c[2]:.4f}")
print()
print("  Since a₁/c₁ < a₃/c₃, the f₀ correction (ρ > 0)")
print("  DECREASES c₁_eff/c₃_eff relative to c₁/c₃.")
print("  This is the RIGHT direction.")

# The observed ratios from RG running (without f₀):
r13_bare = ia1/ia3  # = c₁_run/c₃_run at the matching scale
r23_bare = ia2/ia3

print(f"\n  Bare RG ratios at Λ = {Lam_match:.2e}:")
print(f"    (1/α₁)/(1/α₃) = {r13_bare:.5f}  (need: 2.0000)")
print(f"    (1/α₂)/(1/α₃) = {r23_bare:.5f}  (need: 1.5000)")

# ============================================================
# Step 3: Self-consistent iteration
# ============================================================
print("\n" + "=" * 72)
print("STEP 3: SELF-CONSISTENT ITERATION")
print("=" * 72)

print("""
The problem: we can't just add ρ at one fixed Λ because:
1. The f₀ term changes the boundary conditions
2. Different boundary conditions change the matching scale
3. Different Λ changes yₜ²(Λ)
4. Different yₜ² changes the f₀ correction

Need to solve SIMULTANEOUSLY for (Λ, ρ) such that:
  c₁_eff(ρ, yₜ²(Λ))/c₃_eff(ρ, yₜ²(Λ)) matches the RG-evolved ratio at Λ
""")

# Actually, let me think about this more carefully.
# 
# The spectral action boundary conditions are:
#   1/αᵢ(Λ) = K × cᵢ × [1 + ρ × (aᵢ/cᵢ) × yₜ²(Λ)/(2π²)]
#
# RG running gives us 1/αᵢ(Λ) as functions of Λ.
# So we need:
#   1/αᵢ(Λ) / 1/αⱼ(Λ) = cᵢ × [1 + ρ × (aᵢ/cᵢ) × y/(2π²)] / 
#                          cⱼ × [1 + ρ × (aⱼ/cⱼ) × y/(2π²)]
#
# The LHS is determined by SM RG running (known).
# The RHS depends on ρ and y = yₜ²(Λ).
#
# Two ratios (1/2 and 1/3), two unknowns (Λ, ρ).

def residuals(params):
    """
    Given (log_Lambda, rho), compute:
    R₁₃ = (1/α₁)/(1/α₃) from RG vs from boundary
    R₂₃ = (1/α₂)/(1/α₃) from RG vs from boundary
    """
    log_Lam, rho = params
    
    # Get RG-evolved couplings at Λ
    try:
        ia1, ia2, ia3, yt2 = get_couplings_at_Lambda(log_Lam)
    except:
        return [100, 100]
    
    # RG ratios
    R13_rg = ia1 / ia3
    R23_rg = ia2 / ia3
    
    # Boundary ratios with f₀
    y = yt2
    eps = rho * y / (2*np.pi**2)  # small parameter
    
    # c_eff_i = c_i × (1 + eps × a_i/c_i)
    c1_eff = 8 * (1 + eps * 17/80)
    c2_eff = 6 * (1 + eps * 3/12)
    c3_eff = 4 * (1 + eps * 2/4)
    
    R13_bdy = c1_eff / c3_eff
    R23_bdy = c2_eff / c3_eff
    
    return [R13_rg - R13_bdy, R23_rg - R23_bdy]

# Initial guess: Λ ~ 10⁸, ρ ~ 0.4
from scipy.optimize import fsolve

log_Lam_init = np.log(1e8 / M_Z)
rho_init = 0.4

sol = fsolve(residuals, [log_Lam_init, rho_init], full_output=True)
x_sol, info, ier, msg = sol

if ier == 1:
    log_Lam_sol, rho_sol = x_sol
    Lam_sol = M_Z * np.exp(log_Lam_sol)
    
    ia1, ia2, ia3, yt2 = get_couplings_at_Lambda(log_Lam_sol)
    
    print(f"\nSELF-CONSISTENT SOLUTION FOUND:")
    print(f"  Λ = {Lam_sol:.3e} GeV")
    print(f"  ρ = f₀/(f₂Λ²) = {rho_sol:.4f}")
    print(f"  yₜ²(Λ) = {yt2:.4f}")
    
    # Check effective c values
    eps = rho_sol * yt2 / (2*np.pi**2)
    c1_eff = 8 * (1 + eps * 17/80)
    c2_eff = 6 * (1 + eps * 3/12)
    c3_eff = 4 * (1 + eps * 2/4)
    
    print(f"\n  Effective spectator dimensions:")
    print(f"    c₁_eff = {c1_eff:.4f}  (bare: 8)")
    print(f"    c₂_eff = {c2_eff:.4f}  (bare: 6)")
    print(f"    c₃_eff = {c3_eff:.4f}  (bare: 4)")
    print(f"    c₁_eff/c₃_eff = {c1_eff/c3_eff:.5f}")
    print(f"    c₂_eff/c₃_eff = {c2_eff/c3_eff:.5f}")
    
    # Forward predictions
    K = ia3 / c3_eff  # overall scale
    print(f"\n  Overall scale K = f₂Λ²/(2π²) = {K:.2f}")
    print(f"  f₂Λ² = {K * 2*np.pi**2:.2f}")
    
    # Predict sin²θ_W
    alpha1_pred = 1/ia1
    alpha2_pred = 1/ia2
    sin2_pred = alpha1_pred / (alpha1_pred + alpha2_pred)
    # With standard normalization: sin²θ_W = g'²/(g²+g'²) 
    # = (5/3)α₁/((5/3)α₁ + α₂) ... need to be careful
    # Actually with GUT-normalized α₁: sin²θ_W = (3/8) × c₂/c₁ at tree level
    # More precisely: sin²θ_W = α_em/α₂ = (1/α_em)⁻¹ × (1/α₂)
    # where 1/α_em = (3/5)/α₁ + 1/α₂
    
    inv_alpha_em = 3/5 * ia1 + ia2
    sin2_w = 1/(ia2 * (1/inv_alpha_em))
    # Actually: sin²θ_W = α_em/α₂ where α_em = e²/(4π)
    # 1/α_em = 5/(3α₁) + 1/α₂  (with GUT normalization of α₁)
    # Wait: with the standard GUT normalization where α₁ = (5/3)α_Y:
    # sin²θ_W = g'²/(g²+g'²) = α_Y/(α_Y + α₂)
    # = α₁/(5/3)  / (α₁/(5/3) + α₂)
    # = (3α₁)/(5α₂ + 3α₁)
    sin2_w = (3/ia1) / (5/ia2 + 3/ia1)
    
    print(f"\n  sin²θ_W = {sin2_w:.5f}  (obs: {sin2_obs})")
    print(f"    Deviation: {(sin2_w - sin2_obs)/sin2_obs * 100:.2f}%")
    
    # Check: what does ρ ≈ 0.4 mean physically?
    print(f"\n  Physical interpretation:")
    print(f"    ρ = {rho_sol:.4f} means f₀ correction is {rho_sol*100:.1f}% of f₂Λ²")
    print(f"    This is {'perturbative' if abs(rho_sol) < 0.1 else 'NON-PERTURBATIVE'}")
    if abs(rho_sol) > 0.1:
        print(f"    The Seeley-DeWitt expansion is NOT converging.")
        print(f"    But the self-consistent solution may still exist non-perturbatively.")
else:
    print(f"\n  No solution found: {msg}")
    
    # Scan ρ to see the landscape
    print("\n  Scanning ρ landscape...")
    log_Lam_test = np.log(1e8/M_Z)
    for rho_test in np.arange(-1, 2, 0.1):
        r = residuals([log_Lam_test, rho_test])
        print(f"    ρ = {rho_test:+.1f}: residuals = [{r[0]:+.4f}, {r[1]:+.4f}]")

# ============================================================
# Step 4: Beyond first-order f₀ — full non-perturbative
# ============================================================
print("\n" + "=" * 72)
print("STEP 4: NON-PERTURBATIVE BOUNDARY CONDITIONS")
print("=" * 72)

print("""
If ρ is large, we shouldn't expand to first order. Instead:

The spectral action Tr(f(D²/Λ²)) gives moments:
  f_k = ∫₀^∞ f(x) x^(k/2-1) dx

The FULL boundary condition is:
  1/αᵢ = Σ_k f_k × cᵢ^(k)(D_F) / (normalization)

where cᵢ^(k) involves higher traces of the internal Dirac operator
in each gauge sector.

For k=2 (leading): cᵢ^(2) = cᵢ = {4, 6, 8}
For k=0 (subleading): cᵢ^(0) depends on the specific spectral triple.

In the K₈ spectral triple, cᵢ^(0) involves Yukawa traces that
include the aᵢ coefficients.

The key non-perturbative question: what is the EXACT relation between
cᵢ^(0) and cᵢ^(2) in the K₂ₙ matching algebra?
""")

# Let's parametrize differently: 
# Instead of expanding in ρ, use two independent parameters α and β:
#   1/α₁(Λ) = α × 8 + β × 17/10
#   1/α₂(Λ) = α × 6 + β × 3/2
#   1/α₃(Λ) = α × 4 + β × 2

# Then α = K, β = K × f₀ yₜ²(Λ)/(2 × 2π²) = K × ρ × yₜ²(Λ)/(2 × 2π²)
# But we can solve for α, β directly from the RG couplings.

def solve_alpha_beta(log_Lam):
    ia1, ia2, ia3, yt2 = get_couplings_at_Lambda(log_Lam)
    
    # System: [8  17/10] [α]   [1/α₁]
    #         [6  3/2 ] [β] = [1/α₂]
    #         [4  2   ]       [1/α₃]
    
    # Overdetermined: 3 equations, 2 unknowns
    A = np.array([[8, 17/10], [6, 3/2], [4, 2]])
    b = np.array([ia1, ia2, ia3])
    
    # Least squares
    result = np.linalg.lstsq(A, b, rcond=None)
    alpha, beta = result[0]
    residual = A @ result[0] - b
    
    return alpha, beta, residual, yt2, ia1, ia2, ia3

print("\nScanning Λ for best fit with two-parameter boundary:")
print(f"{'Λ (GeV)':>14} {'α':>10} {'β':>10} {'ρ_eff':>10} {'|residual|':>12} {'c₁_eff':>8} {'c₂_eff':>8}")
print("-" * 80)

best_res = 1e10
best_log_Lam = None

for log_exp in np.arange(6, 12, 0.2):
    log_Lam = np.log(10**log_exp / M_Z)
    try:
        alpha, beta, res, yt2, ia1, ia2, ia3 = solve_alpha_beta(log_Lam)
        res_norm = np.linalg.norm(res)
        
        # Effective c values
        c1_eff = (alpha * 8 + beta * 17/10) / (alpha * 4 + beta * 2) * 4
        c2_eff = (alpha * 6 + beta * 3/2) / (alpha * 4 + beta * 2) * 4
        
        rho_eff = beta / (alpha * yt2 / (2*np.pi**2)) if alpha > 0 and yt2 > 0 else float('inf')
        
        if res_norm < best_res:
            best_res = res_norm
            best_log_Lam = log_Lam
            best_alpha = alpha
            best_beta = beta
            
        if abs(log_exp - round(log_exp)) < 0.01 or res_norm < 0.1:
            print(f"  {10**log_exp:12.2e}  {alpha:10.4f}  {beta:10.4f}  {rho_eff:10.4f}  {res_norm:12.6f}  {c1_eff:8.4f}  {c2_eff:8.4f}")
    except:
        pass

if best_log_Lam is not None:
    print(f"\nBest fit: Λ = {M_Z * np.exp(best_log_Lam):.3e} GeV")
    alpha, beta, res, yt2, ia1, ia2, ia3 = solve_alpha_beta(best_log_Lam)
    print(f"  α = {alpha:.6f}, β = {beta:.6f}")
    print(f"  Residuals: {res}")
    print(f"  ρ_eff = {beta/(alpha * yt2/(2*np.pi**2)):.4f}")
    
    # The residual tells us how well a TWO-parameter model fits THREE data points
    # If residual ≈ 0, the spectral action with f₀ term is self-consistent
    # If residual >> 0, need f₀₀ (third moment) or other corrections
    
    rel_res = np.abs(res / np.array([ia1, ia2, ia3]))
    print(f"\n  Relative residuals:")
    print(f"    U(1):  {rel_res[0]*100:.4f}%")
    print(f"    SU(2): {rel_res[1]*100:.4f}%")
    print(f"    SU(3): {rel_res[2]*100:.4f}%")

# ============================================================
# Step 5: The critical question — is the system exactly solvable?
# ============================================================
print("\n" + "=" * 72)
print("STEP 5: EXACT SOLVABILITY AT THE MATCHING SCALE")
print("=" * 72)

# Find Λ where the 3-equation/2-unknown system is EXACTLY consistent
# (residual = 0)
# This would mean the spectral action with TWO moments (f₂, f₀)
# exactly determines all three gauge couplings.

from scipy.optimize import minimize_scalar

def residual_norm(log_Lam):
    try:
        _, _, res, _, _, _, _ = solve_alpha_beta(log_Lam)
        return np.linalg.norm(res)
    except:
        return 1e10

# Fine scan
result = minimize_scalar(residual_norm, 
                         bounds=(np.log(1e6/M_Z), np.log(1e12/M_Z)),
                         method='bounded')

log_Lam_opt = result.x
Lam_opt = M_Z * np.exp(log_Lam_opt)
alpha, beta, res, yt2, ia1, ia2, ia3 = solve_alpha_beta(log_Lam_opt)

print(f"\nOptimal Λ (minimum residual): {Lam_opt:.4e} GeV")
print(f"  Minimum |residual| = {np.linalg.norm(res):.6e}")
print(f"  α = {alpha:.6f}, β = {beta:.6f}")
print(f"  yₜ²(Λ) = {yt2:.4f}")

rel_res = np.abs(res / np.array([ia1, ia2, ia3]))
print(f"\n  Relative residuals at optimum:")
print(f"    U(1):  {rel_res[0]*100:.4f}%")
print(f"    SU(2): {rel_res[1]*100:.4f}%")
print(f"    SU(3): {rel_res[2]*100:.4f}%")

rho_opt = beta / (alpha * yt2 / (2*np.pi**2)) if alpha > 0 and yt2 > 0 else float('inf')
print(f"\n  ρ_eff = {rho_opt:.4f}")

# Forward prediction at optimum
c1_eff_opt = (alpha * 8 + beta * 17/10) / (alpha * 4 + beta * 2) * 4
c2_eff_opt = (alpha * 6 + beta * 3/2) / (alpha * 4 + beta * 2) * 4
print(f"\n  Effective spectator dimensions at optimum:")
print(f"    c₁_eff = {c1_eff_opt:.4f}")
print(f"    c₂_eff = {c2_eff_opt:.4f}")

# sin²θ_W prediction
# With the two-parameter fit, predict all couplings
ia1_pred = alpha * 8 + beta * 17/10
ia2_pred = alpha * 6 + beta * 3/2
ia3_pred = alpha * 4 + beta * 2

# Use the fit values (which include the f₀ correction) for sin²θ_W
sin2_pred_opt = (3/ia1_pred) / (5/ia2_pred + 3/ia1_pred)

print(f"\n  sin²θ_W from self-consistent boundary:")
print(f"    Predicted: {sin2_pred_opt:.5f}")
print(f"    Observed:  {sin2_obs:.5f}")
print(f"    Deviation: {(sin2_pred_opt - sin2_obs)/sin2_obs*100:.3f}%")

# ============================================================
# Step 6: What the residual tells us
# ============================================================
print("\n" + "=" * 72)
print("STEP 6: INTERPRETATION")
print("=" * 72)

min_res = np.linalg.norm(res)
if min_res < 0.01 * min(ia1, ia2, ia3):
    print("""
  RESULT: The two-parameter spectral action (f₂, f₀) is 
  EXACTLY CONSISTENT with all three gauge couplings at a 
  specific Λ. The 2.5% residual is resolved by the f₀ term.
""")
else:
    print(f"""
  RESULT: The two-parameter spectral action (f₂, f₀) gives a 
  MINIMUM residual of {min_res:.4e} at Λ = {Lam_opt:.2e} GeV.

  Relative to the coupling values, this is:
    {max(rel_res)*100:.3f}% in the worst channel.

  If this is small (< 0.1%): the f₀ correction essentially 
  closes the gap. The 2.5% residual was due to ignoring f₀.

  If this is large (> 1%): a THIRD moment f₋₂ or threshold 
  corrections are needed.
""")

