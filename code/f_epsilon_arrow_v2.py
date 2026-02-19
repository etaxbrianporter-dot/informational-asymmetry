#!/usr/bin/env python3
"""
f(Îµ) with Arrow of Time: Positive Energy Only
================================================

The particle-hole symmetric computation (Â±E) gives f ~ 1/Îµâ´
because âŸ¨E|kâŸ© = 0 for all k â€” no information in the mean.

The arrow of time (â„¤â‚ƒ flux â†’ C = âˆ’2 â†’ D â‰  D*) breaks this
symmetry. A physical observer sees positive-energy excitations
propagating forward in time.

For E > 0 only: âŸ¨E|kâŸ© = Eâ‚(k) varies from 0 to 1.73.
Var_k[Eâ‚] = 0.1735 (nonzero!).
Leading term: f(Îµ) ~ Var_k[Eâ‚] / (2ÎµÂ² Â· H(E)) ~ 1/ÎµÂ²

This script computes all three cases:
  (A) Single channel, E > 0 only (physical: arrow of time)
  (B) Single channel, Â±E (particle-hole symmetric)
  (C) Analytical prediction for coefficient

The arrow of time doesn't just give time a direction.
It gives dark energy a magnitude.

Brian Porter â€” February 2026
"""

import numpy as np
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

omega = np.exp(2j * np.pi / 3)
t_dem = 1.0 / np.sqrt(3)

def cone_energy(k1, k2):
    """Single gapless channel energy Eâ‚(k) = |dâ‚(k)| â‰¥ 0."""
    e1 = np.exp(1j * k1)
    e2 = np.exp(1j * k2)
    d1 = t_dem * (1 + omega * e1 + omega**2 * e2)
    return np.abs(d1)

# =============================================================================
# PART 1: VARIANCE STRUCTURE
# =============================================================================

print("="*72)
print("f(Îµ) with Arrow of Time: Positive Energy Only")
print("="*72)

Nk_fine = 400
k1s = np.linspace(-np.pi, np.pi, Nk_fine, endpoint=False)
k2s = np.linspace(-np.pi, np.pi, Nk_fine, endpoint=False)

E1_grid = np.zeros((Nk_fine, Nk_fine))
for i, k1 in enumerate(k1s):
    for j, k2 in enumerate(k2s):
        E1_grid[i,j] = cone_energy(k1, k2)

E1_flat = E1_grid.flatten()

print(f"\n--- Cone energy statistics ---")
print(f"  âŸ¨Eâ‚âŸ©_BZ       = {np.mean(E1_flat):.8f}")
print(f"  âŸ¨Eâ‚Â²âŸ©_BZ      = {np.mean(E1_flat**2):.8f}")
print(f"  Var_k[Eâ‚]      = {np.var(E1_flat):.8f}  â† THIS drives 1/ÎµÂ²")
print(f"  Var_k[Eâ‚Â²]     = {np.var(E1_flat**2):.8f}  â† THIS drove 1/Îµâ´")
print(f"")
print(f"  Â±E case:  âŸ¨EâŸ©_k = 0 for all k â†’ leading term ~ Var[EÂ²]/Îµâ´")
print(f"  +E case:  âŸ¨EâŸ©_k = Eâ‚(k) varies â†’ leading term ~ Var[Eâ‚]/ÎµÂ²")
print(f"")
print(f"  Arrow of time promotes scaling: 1/Îµâ´ â†’ 1/ÎµÂ²")

# =============================================================================
# PART 2: COMPUTE f(Îµ) â€” POSITIVE ENERGY ONLY
# =============================================================================

def compute_f_positive_E(epsilon, Nk=100, NE=300, E_range_sigmas=6):
    """
    f(Îµ) for POSITIVE ENERGY ONLY (arrow of time).
    A(k,E) = G_Îµ(E - Eâ‚(k)) for E â‰¥ 0.
    """
    k1_grid = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
    k2_grid = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
    dk = (2*np.pi / Nk)**2
    
    E_min = -3 * epsilon  # buffer for Gaussian tails
    E_max = np.max(E1_flat) + E_range_sigmas * epsilon
    E_grid = np.linspace(E_min, E_max, NE)
    dE = E_grid[1] - E_grid[0]
    
    H_Ek_sum = 0.0
    p_E_marginal = np.zeros(NE)
    
    for k1 in k1_grid:
        for k2 in k2_grid:
            E_cone = cone_energy(k1, k2)
            
            # Single peak at E = Eâ‚(k), positive energy
            A_k = np.exp(-0.5 * ((E_grid - E_cone) / epsilon)**2)
            
            # Enforce E â‰¥ 0
            A_k[E_grid < 0] = 0.0
            
            A_k_sum = np.sum(A_k) * dE
            if A_k_sum > 1e-300:
                p_Ek = A_k / A_k_sum
                mask = p_Ek > 1e-300
                H_k = -np.sum(p_Ek[mask] * np.log(p_Ek[mask]) * dE)
                H_Ek_sum += H_k * dk
            
            p_E_marginal += A_k * dk
    
    p_E_sum = np.sum(p_E_marginal) * dE
    if p_E_sum > 1e-300:
        p_E_marginal /= p_E_sum
    
    mask = p_E_marginal > 1e-300
    H_E = -np.sum(p_E_marginal[mask] * np.log(p_E_marginal[mask]) * dE)
    
    BZ_area = (2*np.pi)**2
    H_Ek = H_Ek_sum / BZ_area
    
    if H_E > 1e-300:
        d_eff = 2.0 + H_Ek / H_E
        f_eps = 3.0 - d_eff
    else:
        d_eff = 2.0
        f_eps = 1.0
    
    return d_eff, f_eps, H_E, H_Ek


def compute_f_pm_E(epsilon, Nk=100, NE=300, E_range_sigmas=6):
    """
    f(Îµ) for Â±E (particle-hole symmetric, the 1/Îµâ´ case).
    """
    k1_grid = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
    k2_grid = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
    dk = (2*np.pi / Nk)**2
    
    E_range = np.max(E1_flat) + E_range_sigmas * epsilon
    E_grid = np.linspace(-E_range, E_range, NE)
    dE = E_grid[1] - E_grid[0]
    
    H_Ek_sum = 0.0
    p_E_marginal = np.zeros(NE)
    
    for k1 in k1_grid:
        for k2 in k2_grid:
            E_cone = cone_energy(k1, k2)
            A_k = (np.exp(-0.5 * ((E_grid - E_cone) / epsilon)**2) +
                   np.exp(-0.5 * ((E_grid + E_cone) / epsilon)**2))
            
            A_k_sum = np.sum(A_k) * dE
            if A_k_sum > 1e-300:
                p_Ek = A_k / A_k_sum
                mask = p_Ek > 1e-300
                H_k = -np.sum(p_Ek[mask] * np.log(p_Ek[mask]) * dE)
                H_Ek_sum += H_k * dk
            
            p_E_marginal += A_k * dk
    
    p_E_sum = np.sum(p_E_marginal) * dE
    if p_E_sum > 1e-300:
        p_E_marginal /= p_E_sum
    
    mask = p_E_marginal > 1e-300
    H_E = -np.sum(p_E_marginal[mask] * np.log(p_E_marginal[mask]) * dE)
    
    BZ_area = (2*np.pi)**2
    H_Ek = H_Ek_sum / BZ_area
    
    if H_E > 1e-300:
        d_eff = 2.0 + H_Ek / H_E
        f_eps = 3.0 - d_eff
    else:
        d_eff = 2.0
        f_eps = 1.0
    
    return d_eff, f_eps, H_E, H_Ek

# =============================================================================
# PART 3: SCAN â€” POSITIVE ENERGY (ARROW OF TIME)
# =============================================================================

epsilon_values = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7,
                  1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0,
                  30.0, 50.0, 70.0, 100.0]

print("\n" + "="*72)
print("CASE A: Positive Energy Only (Arrow of Time)")
print("="*72)
print(f"\n{'Îµ':>8s} {'d_eff':>12s} {'f(Îµ)':>14s} {'ÎµÂ²Â·f':>12s}")
print("-" * 50)

results_pos = []
for eps in epsilon_values:
    if eps < 0.05:
        Nk, NE = 80, 500
    elif eps < 0.5:
        Nk, NE = 100, 400
    elif eps < 5.0:
        Nk, NE = 120, 300
    else:
        Nk, NE = 150, 400
    
    d_eff, f_eps, H_E, H_Ek = compute_f_positive_E(eps, Nk=Nk, NE=NE)
    results_pos.append((eps, d_eff, f_eps, H_E, H_Ek))
    
    if f_eps > 1e-15:
        eps2f = eps**2 * f_eps
        print(f"{eps:8.3f} {d_eff:12.8f} {f_eps:14.8e} {eps2f:12.8e}")
    else:
        print(f"{eps:8.3f} {d_eff:12.8f} {'<noise':>14s} {'---':>12s}")

# =============================================================================
# PART 4: HEAD-TO-HEAD COMPARISON
# =============================================================================

print("\n" + "="*72)
print("CASE B: Head-to-Head â€” Arrow of Time vs Particle-Hole Symmetric")
print("="*72)
print(f"\n{'Îµ':>8s} {'f(+E)':>14s} {'f(Â±E)':>14s} {'ÎµÂ²Â·f(+E)':>12s} {'Îµâ´Â·f(Â±E)':>12s}")
print("-" * 66)

eps_compare = [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 50.0]

for eps in eps_compare:
    if eps < 5.0:
        Nk, NE = 120, 300
    else:
        Nk, NE = 150, 400
    
    _, f_pos, _, _ = compute_f_positive_E(eps, Nk=Nk, NE=NE)
    _, f_pm, _, _ = compute_f_pm_E(eps, Nk=Nk, NE=NE)
    
    e2fp = eps**2 * f_pos if f_pos > 1e-15 else 0
    e4fm = eps**4 * f_pm if f_pm > 1e-15 else 0
    
    fp_str = f"{f_pos:.6e}" if f_pos > 1e-15 else "<noise"
    fm_str = f"{f_pm:.6e}" if f_pm > 1e-15 else "<noise"
    e2_str = f"{e2fp:.6e}" if e2fp > 0 else "---"
    e4_str = f"{e4fm:.6e}" if e4fm > 0 else "---"
    
    print(f"{eps:8.3f} {fp_str:>14s} {fm_str:>14s} {e2_str:>12s} {e4_str:>12s}")

# =============================================================================
# PART 5: SCALING ANALYSIS â€” POSITIVE ENERGY
# =============================================================================

print("\n" + "="*72)
print("PART 5: Scaling Analysis (Positive Energy)")
print("="*72)

results_arr = np.array(results_pos)
eps_arr = results_arr[:, 0]
f_arr = results_arr[:, 2]

# Local slopes
print("\n  Local power law exponent: Î±(Îµ) = -d(ln f)/d(ln Îµ)")
print(f"{'Îµ_mid':>10s} {'Î±_local':>10s}")
print("-" * 25)

slopes = []
eps_mids = []
for i in range(len(results_arr)-1):
    eps1, _, f1, _, _ = results_arr[i]
    eps2, _, f2, _, _ = results_arr[i+1]
    
    if f1 > 1e-14 and f2 > 1e-14 and eps1 >= 0.3:
        alpha_local = -(np.log(f2) - np.log(f1)) / (np.log(eps2) - np.log(eps1))
        eps_mid = np.sqrt(eps1 * eps2)
        slopes.append(alpha_local)
        eps_mids.append(eps_mid)
        print(f"{eps_mid:10.3f} {alpha_local:10.4f}")

# Tail fits (Îµ â‰¥ 2)
mask_tail = (eps_arr >= 2.0) & (f_arr > 1e-14)
eps_tail = eps_arr[mask_tail]
f_tail = f_arr[mask_tail]

if len(eps_tail) >= 3:
    log_eps = np.log(eps_tail)
    log_f = np.log(np.maximum(f_tail, 1e-300))
    
    if np.all(np.isfinite(log_f)):
        # Free power law
        coeffs = np.polyfit(log_eps, log_f, 1)
        alpha_fit = -coeffs[0]
        A_fit = np.exp(coeffs[1])
        f_pred = A_fit / eps_tail**alpha_fit
        resid = np.mean(np.abs(f_tail - f_pred) / f_tail)
        print(f"\n  Free power law: f(Îµ) = {A_fit:.6e} / Îµ^{alpha_fit:.4f}")
        print(f"  Mean residual: {resid:.4f}")
        
        # Pure 1/ÎµÂ²
        C2 = np.mean(f_tail * eps_tail**2)
        f_pred2 = C2 / eps_tail**2
        resid2 = np.mean(np.abs(f_tail - f_pred2) / f_tail)
        print(f"\n  Pure 1/ÎµÂ²: f(Îµ) = {C2:.6e} / ÎµÂ²")
        print(f"  Mean residual: {resid2:.4f}")
        
        # 1/(ÎµÂ² log Îµ) â€” expected from analytical expansion
        try:
            mask_log = eps_tail > 2.5
            if np.sum(mask_log) >= 3:
                eps_log = eps_tail[mask_log]
                f_log = f_tail[mask_log]
                C_log = np.mean(f_log * eps_log**2 * np.log(eps_log))
                f_pred_log = C_log / (eps_tail**2 * np.log(eps_tail))
                resid_log = np.mean(np.abs(f_tail - f_pred_log) / f_tail)
                print(f"\n  1/(ÎµÂ² ln Îµ): f(Îµ) = {C_log:.6e} / (ÎµÂ² ln Îµ)")
                print(f"  Mean residual: {resid_log:.4f}")
        except:
            pass
        
        # 1/ÎµÂ² + 1/Îµâ´
        try:
            def pw24(x, a, b):
                return a / x**2 + b / x**4
            popt24, _ = curve_fit(pw24, eps_tail, f_tail, p0=[C2, 0.01])
            f_pred24 = pw24(eps_tail, *popt24)
            resid24 = np.mean(np.abs(f_tail - f_pred24) / f_tail)
            print(f"\n  1/ÎµÂ² + 1/Îµâ´: f(Îµ) = {popt24[0]:.6e}/ÎµÂ² + {popt24[1]:.6e}/Îµâ´")
            print(f"  Mean residual: {resid24:.4f}")
        except:
            pass
        
        # Exponential (for comparison â€” should be bad)
        try:
            def exp_model(x, B, beta):
                return B * np.exp(-beta * x)
            popt_e, _ = curve_fit(exp_model, eps_tail, f_tail, p0=[0.1, 0.5], maxfev=10000)
            f_pred_e = exp_model(eps_tail, *popt_e)
            resid_e = np.mean(np.abs(f_tail - f_pred_e) / f_tail)
            print(f"\n  Exponential: f(Îµ) = {popt_e[0]:.6e} * exp(-{popt_e[1]:.4f}Îµ)")
            print(f"  Mean residual: {resid_e:.4f}")
        except:
            resid_e = 999

# =============================================================================
# PART 6: THE DEFINITIVE COLUMN â€” ÎµÂ²Â·fÂ·ln(Îµ)
# =============================================================================

print("\n" + "="*72)
print("PART 6: Stabilization Tests")
print("="*72)

print(f"""
  If f ~ C/ÎµÂ²:           ÎµÂ²Â·f â†’ C (constant)
  If f ~ C/(ÎµÂ² ln Îµ):    ÎµÂ²Â·fÂ·ln(Îµ) â†’ C (constant)  â† analytical prediction
  If f ~ C/Îµâ´:           Îµâ´Â·f â†’ C (constant)
""")

print(f"{'Îµ':>8s} {'ÎµÂ²Â·f':>14s} {'ÎµÂ²Â·fÂ·ln(Îµ)':>14s} {'Îµâ´Â·f':>14s}")
print("-" * 55)

for eps, d_eff, f_eps, H_E, H_Ek in results_pos:
    if f_eps > 1e-15 and eps >= 0.5:
        e2f = eps**2 * f_eps
        e2f_log = e2f * np.log(eps) if eps > 1 else 0
        e4f = eps**4 * f_eps
        log_str = f"{e2f_log:14.8e}" if eps > 1 else f"{'---':>14s}"
        print(f"{eps:8.3f} {e2f:14.8e} {log_str} {e4f:14.8e}")

# =============================================================================
# PART 7: COSMOLOGICAL CONSTANT
# =============================================================================

Var_E1 = np.var(E1_flat)

print("\n" + "="*72)
print("PART 7: Cosmological Constant Estimate")
print("="*72)

print(f"""
  If f(Îµ) â‰ˆ Var_k[Eâ‚] / (2ÎµÂ² Â· H(E)) and H(E) ~ ln(Îµ) at large Îµ:
  
    f(Îµ) ~ {Var_E1:.6f} / (2 ÎµÂ² ln Îµ)
  
  For Hubble-scale observer: Îµ_H = L_H/L_Pl â‰ˆ 10â¶Â¹
  
    f(Îµ_H) â‰ˆ {Var_E1:.4f} / (2 Ã— 10Â¹Â²Â² Ã— 140.5)
           â‰ˆ {Var_E1/(2 * 140.5):.4e} Ã— 10â»Â¹Â²Â²
           â‰ˆ {Var_E1/(2 * 140.5) * 1e122:.2f} Ã— 10â»Â¹Â²Â²
  
  Observed: Î›_CC/Î›_Pl â‰ˆ 2.9 Ã— 10â»Â¹Â²Â²
  
  The 10â»Â¹Â²Â² comes from 1/ÎµÂ² = (L_Pl/L_H)Â². 
  Zero free parameters. Var_k[Eâ‚] = {Var_E1:.6f} is a Kâ‚„ invariant.
""")

# =============================================================================
# PART 8: VERDICT
# =============================================================================

print("="*72)
print("VERDICT")  
print("="*72)

if len(slopes) >= 3:
    tail_slopes = [s for s, e in zip(slopes, eps_mids) if e > 3.0]
    if len(tail_slopes) >= 2:
        mean_slope = np.mean(tail_slopes)
        std_slope = np.std(tail_slopes)
        
        print(f"\n  Local slope (Îµ > 3): Î± = {mean_slope:.4f} Â± {std_slope:.4f}")
        
        if abs(mean_slope - 2.0) < 0.5:
            print("""
  â˜…â˜…â˜… 1/ÎµÂ² CONFIRMED â˜…â˜…â˜…
  
  The arrow of time promotes f(Îµ) from 1/Îµâ´ to 1/ÎµÂ².
  The four-consequence theorem STANDS.
  Î›_CC ~ (L_Pl/L_H)Â² ~ 10â»Â¹Â²Â² from first principles.
  
  Axiom 3 gives dark energy. Axiom 4 gives its magnitude.
""")
        elif abs(mean_slope - 2.0) < 1.0:
            print(f"""
  â˜… NEAR 1/ÎµÂ² (Î± = {mean_slope:.2f})
  
  Likely 1/(ÎµÂ² Â· (ln Îµ)^p) with logarithmic correction.
  This is STILL power-law. The 10â»Â¹Â²Â² survives.
  The log correction modifies the prefactor, not the exponent.
  
  Check: does ÎµÂ²Â·fÂ·(ln Îµ)^p stabilize for p â‰ˆ {mean_slope - 2:.1f}?
""")
        elif mean_slope < 3.5:
            print(f"""
  â˜… POWER LAW but Î± â‰ˆ {mean_slope:.1f}, not 2.0.
  
  This could be:
  (a) Genuine Î± â‰ˆ {mean_slope:.1f} â†’ Î›_CC ~ 10^{{-{61*mean_slope:.0f}}} (wrong)
  (b) 1/ÎµÂ² with large log correction not yet in asymptotic regime
  (c) Crossover effect from finite bandwidth
  
  Need: larger Îµ range (Îµ > 100) with high-precision arithmetic,
  or analytical control of the entropy expansion to all orders.
""")
        else:
            print(f"""
  Î± â‰ˆ {mean_slope:.1f} â€” too fast for dark energy.
  Three consequences survive. Dark energy not explained.
""")

print("\n" + "="*72)
print("THE CHAIN")
print("="*72)
print("""
  Â±E (no arrow):   âŸ¨E|kâŸ© = 0      â†’ Var = 0    â†’ f ~ 1/Îµâ´ â†’ no Î›_CC
  +E (arrow):      âŸ¨E|kâŸ© = Eâ‚(k)  â†’ Var > 0    â†’ f ~ 1/ÎµÂ² â†’ Î›_CC âœ“
  
  What breaks the symmetry: Axiom 3 (D â‰  D*)
  What sets the scale:      Axiom 4 (Îµ = L_obs/L_Pl)
  What fixes the coefficient: Kâ‚„ band structure (Var_k[Eâ‚])
  
  Dark energy = the arrow of time Ã— observer finitude Ã— Kâ‚„ geometry.
""")
