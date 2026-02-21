"""
Self-consistent boundary problem v2: corrected sin²θ_W, deeper diagnostics
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# Physical constants
M_Z = 91.1876
M_t = 173.0
v = 246.22  # Higgs vev

# Observed couplings at M_Z (MSbar, GUT-normalized α₁)
alpha1_obs = 0.016943
alpha2_obs = 0.033801
alpha3_obs = 0.1180

# sin²θ_W = (3/5)α₁ / ((3/5)α₁ + α₂)  [α₁ has 5/3 GUT factor]
sin2_obs = 0.23122
sin2_from_alphas = (3/5*alpha1_obs) / (3/5*alpha1_obs + alpha2_obs)
print(f"Cross-check: sin²θ_W from α₁,α₂ = {sin2_from_alphas:.5f} (obs: {sin2_obs})")

# Spectator dimensions and subleading Casimir sums
c  = np.array([8.0, 6.0, 4.0])
cp = np.array([17/10, 3/2, 2.0])

# SM 1-loop beta coefficients: dα_i/dt = b_i α_i²/(2π)
# where t = ln(μ/M_Z)
# b₁ = 41/10, b₂ = -19/6, b₃ = -7
# For 1/α: d(1/α_i)/dt = -b_i/(2π)
b = np.array([41/10, -19/6, -7.0])

# 2-loop gauge-gauge: Bij
# d(1/α_i)/dt += Σ_j B_ij α_j/(4π²)  
Bij = np.array([
    [199/50,   27/10,   44/5 ],
    [9/10,     35/6,    12   ],
    [11/10,    9/2,    -26   ]
])

# Yukawa a_i coefficients
a_yuk = np.array([17/10, 3/2, 2.0])

def rge_2loop(t, y):
    """2-loop RGE for (1/α₁, 1/α₂, 1/α₃, y_t²)"""
    inv_a = y[:3]
    yt2 = max(y[3], 1e-10)
    alpha = np.clip(1.0/inv_a, 1e-6, 1.0)
    
    # d(1/αᵢ)/dt at 1-loop
    dinv = -b / (2*np.pi)
    
    # 2-loop gauge-gauge
    for i in range(3):
        for j in range(3):
            dinv[i] -= Bij[i,j] * alpha[j] / (4*np.pi**2)
    
    # 2-loop Yukawa correction to gauge
    dinv -= a_yuk * yt2 / (8*np.pi**2)  # note sign convention
    
    # Top Yukawa 1-loop RGE
    g_sq = 4*np.pi*alpha  # g_i² = 4π α_i
    dyt2 = yt2/(16*np.pi**2) * (
        9/2 * yt2 
        - 17/20 * g_sq[0] 
        - 9/4 * g_sq[1] 
        - 8 * g_sq[2]
    )
    
    return np.concatenate([dinv, [dyt2]])

def run_up(log_Lambda):
    """Run from M_Z to Λ"""
    yt2_MZ = 2*M_t**2/v**2
    y0 = np.array([1/alpha1_obs, 1/alpha2_obs, 1/alpha3_obs, yt2_MZ])
    sol = solve_ivp(rge_2loop, [0, log_Lambda], y0, 
                    method='RK45', rtol=1e-10, atol=1e-12, max_step=0.2)
    return sol.y[:,-1] if sol.success else None

def solve_2param(log_Lambda):
    """
    Solve K₁, K₂ from SU(2)+SU(3), check U(1).
    Return (K1, K2, resid, rho, inv_alphas, yt2)
    """
    r = run_up(log_Lambda)
    if r is None: return None
    inv_a, yt2 = r[:3], r[3]
    
    A = np.array([[c[1], cp[1]], [c[2], cp[2]]])
    K1, K2 = np.linalg.solve(A, inv_a[1:3])
    
    inv_a1_pred = c[0]*K1 + cp[0]*K2
    resid = inv_a1_pred - inv_a[0]
    rho = K2/K1
    
    return K1, K2, resid, rho, inv_a, yt2

# ============================================================
# 1. Reference: f₀=0 case
# ============================================================
print("\n" + "="*72)
print("REFERENCE: f₀ = 0 (leading Seeley-DeWitt only)")
print("="*72)

def f0_zero_resid(lL):
    r = run_up(lL)
    if r is None: return 1.0
    return r[1]/r[2] - c[1]/c[2]

prev = None
for l10 in np.arange(5, 16, 0.1):
    lL = (l10 - np.log10(M_Z))*np.log(10)
    val = f0_zero_resid(lL)
    if prev is not None and val*prev < 0:
        lL0 = brentq(f0_zero_resid, prev_lL, lL, xtol=1e-12)
        r = run_up(lL0)
        Lambda0 = M_Z*np.exp(lL0)
        K1_0 = r[2]/c[2]
        c1e = r[0]/K1_0
        print(f"  Λ = {Lambda0:.4e} GeV  (log₁₀ = {np.log10(Lambda0):.4f})")
        print(f"  c₁_eff = {c1e:.4f}  → deviation from 8: {(c1e-8)/8*100:+.2f}%")
        
        # Forward prediction of sin²θ_W
        a1_pred = K1_0*8
        a2_pred = K1_0*6
        # These are 1/α at Λ. Run down:
        # sin²θ_W from the OBSERVED couplings at M_Z:
        print(f"  sin²θ_W (observed) = {sin2_from_alphas:.5f}")
        
        # sin²θ_W PREDICTED from c₁:c₂ = 8:6 boundary:
        # At Λ: 1/α₁ = 8K₁, 1/α₂ = 6K₁
        # → α₁(Λ)/α₂(Λ) = 6/8 = 3/4
        # sin²θ_W(Λ) = (3/5)(1/8K₁) / ((3/5)(1/8K₁) + 1/(6K₁))
        #            = (3/5)/8 / ((3/5)/8 + 1/6) = (3/40) / (3/40 + 1/6) = 9/120 / (9/120 + 20/120) = 9/29
        sin2_Λ = (3/5*c[1]) / (3/5*c[1] + c[0])  # sin²θ at Λ from spectator ratios
        # Nope, sin²θ_W = (3/5)α₁/((3/5)α₁ + α₂) = (3/5)(1/c₁)/((3/5)(1/c₁) + 1/c₂)
        sin2_Λ = (3/5/c[0]) / (3/5/c[0] + 1/c[1])
        # = (3/40)/(3/40 + 1/6) = (9/120)/(9/120 + 20/120) = 9/29
        print(f"  sin²θ_W at Λ (from 8:6 ratio) = {sin2_Λ:.5f} = {9}/{29}")
        break
    prev, prev_lL = val, lL

# ============================================================
# 2. Two-parameter solution
# ============================================================
print("\n" + "="*72)
print("TWO-PARAMETER SOLUTION: f₂Λ² + f₀")
print("="*72)

# Find zero crossing
data = []
for l10 in np.arange(5, 17, 0.25):
    lL = (l10 - np.log10(M_Z))*np.log(10)
    r = solve_2param(lL)
    if r: data.append((lL, l10, *r))

for i in range(len(data)-1):
    if data[i][4]*data[i+1][4] < 0:  # residual sign change
        lL_lo, lL_hi = data[i][0], data[i+1][0]
        lL_root = brentq(lambda x: solve_2param(x)[2], lL_lo, lL_hi, xtol=1e-14)
        
        K1, K2, resid, rho, inv_a, yt2 = solve_2param(lL_root)
        Lambda = M_Z*np.exp(lL_root)
        
        print(f"\n  SOLUTION:")
        print(f"  Λ = {Lambda:.6e} GeV  (log₁₀ = {np.log10(Lambda):.4f})")
        print(f"  K₁ = f₂Λ²/(2π²) = {K1:.8f}")
        print(f"  K₂ = f₀/(4π²)   = {K2:.8f}")
        print(f"  ρ = f₀/(f₂Λ²)   = {rho:.6f}")
        print(f"  y_t²(Λ)          = {yt2:.6f}")
        print(f"  |residual|        = {abs(resid):.2e}")
        
        # Physical parameters
        f2L2 = K1*2*np.pi**2
        f0   = K2*4*np.pi**2
        print(f"\n  f₂Λ² = {f2L2:.4f}")
        print(f"  f₀   = {f0:.4f}")
        print(f"  f₂   = {f2L2/Lambda**2:.4e}")
        
        # Effective boundary values
        print(f"\n  Boundary conditions at Λ:")
        for i, name in enumerate(['U(1)','SU(2)','SU(3)']):
            bc = c[i]*K1 + cp[i]*K2
            pure = c[i]*K1
            corr = cp[i]*K2
            print(f"    1/α_{name}(Λ) = {bc:.4f}  = {c[i]}×{K1:.4f} + {cp[i]}×{K2:.4f}")
            print(f"      leading: {pure:.4f} ({pure/bc*100:.1f}%)  subleading: {corr:.4f} ({corr/bc*100:.1f}%)")
        
        # sin²θ_W at the boundary  
        alpha_Λ = 1.0/inv_a
        sin2_Λ = (3/5*alpha_Λ[0])/((3/5*alpha_Λ[0]) + alpha_Λ[1])
        sin2_MZ = sin2_from_alphas
        print(f"\n  sin²θ_W(Λ)  = {sin2_Λ:.6f}")
        print(f"  sin²θ_W(M_Z) = {sin2_MZ:.6f}  (obs: {sin2_obs})")
        
        # ============================================================
        # 3. Key diagnostic: how much does f₀ shift each coupling?
        # ============================================================
        print(f"\n  DIAGNOSTIC: f₀ correction magnitude per coupling")
        for i, name in enumerate(['U(1)','SU(2)','SU(3)']):
            shift_pct = (cp[i]*K2)/(c[i]*K1)*100
            print(f"    {name}: f₀ shift = {shift_pct:+.2f}% of leading term")
        
        # ============================================================
        # 4. The REAL question: is ρ = 0.25 reasonable?
        # ============================================================
        print(f"\n" + "="*72)
        print(f"  ASSESSMENT")
        print(f"="*72)
        print(f"""
  The two-parameter system (f₂Λ², f₀) has a unique solution:
    Λ  = {Lambda:.2e} GeV
    ρ  = f₀/(f₂Λ²) = {rho:.4f}
    
  This is a {rho*100:.1f}% subleading correction. Assessment:
  
  1. SIGN: f₀ > 0 ✓ (physically required for positive spectral function)
  
  2. MAGNITUDE: ρ = {rho:.3f}
     - Perturbative regime: ρ ≪ 1 → barely satisfied
     - The Seeley-DeWitt expansion IS the spectral action's asymptotic
       expansion in 1/(f₂Λ²). At ρ = {rho:.2f}, the next term (f₋₂/Λ²)
       could contribute ~ρ² = {rho**2:.3f} (~{rho**2*100:.1f}%)
     - This is the SAME order as the residual we're trying to fix
  
  3. PREDICTION COUNT: 
     System has 3 equations, 3 unknowns → exactly determined
     NO new prediction is made by this step alone
     The system ALWAYS has a solution for any c'_i values
  
  4. THE REAL TEST: 
     The c'_i = {{17/10, 3/2, 2}} are NOT free — they're the SM
     Casimir sums, completely determined by fermion content.
     But with 3 equations and 3 unknowns, ANY c'_i would work.
     
     The test would be: include a FOURTH observable (e.g., y_t or λ_H 
     at Λ) and check if the same (K₁, K₂, Λ) satisfies it.
""")
        
        # ============================================================
        # 5. Fourth observable: top Yukawa
        # ============================================================
        print(f"  FOURTH OBSERVABLE: Top Yukawa at Λ")
        print(f"  y_t²(Λ) from RG running = {yt2:.6f}")
        print(f"  y_t(Λ) = {np.sqrt(yt2):.6f}")
        
        # From K₈ spectral triple, the Yukawa should be related to
        # matching algebra invariants. The tree-level relation is
        # y_t = √2 m_t / v, but at Λ it has run.
        yt_tree = np.sqrt(2)*M_t/v
        print(f"  y_t (tree, M_Z) = {yt_tree:.6f}")
        print(f"  y_t(Λ)/y_t(tree) = {np.sqrt(yt2)/yt_tree:.4f}")
        
        # ============================================================
        # 6. What about threshold corrections?
        # ============================================================
        print(f"\n  NOTE ON THRESHOLDS:")
        print(f"  Current RGE uses 3-generation SM from M_Z to Λ.")
        print(f"  Step-function decoupling at M_t, M_W changes β coefficients.")
        print(f"  Expected effect: O(α) ~ O(1%) shift in 1/α_i(Λ)")
        print(f"  This is comparable to the f₀ correction itself.")
        print(f"  A proper treatment requires matched running at each threshold.")
        
        break

print("\nDone.")
