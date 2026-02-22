#!/usr/bin/env python3
"""
Geometric Dark Matter Probe v2 — CORRECTED
===========================================
Uses the CORRECT K₄ band structure (ℤ₃ flux!) and d_eff method 
(mutual information, not return probability) from d_eff_campaign_v3.

The K₄ Dirac cone with ℤ₃ flux:
  d₁(k) = t₁ + t₂·ω·exp(ik₁) + t₃·ω²·exp(ik₂)
  ω = exp(2πi/3), democratic point: t₁ = t₂ = t₃ = 1/√3

d_eff = 2 + H(E|k)/H(E) via Gaussian spectral function at resolution ε.

Question: how does f(ε) = 3 - d_eff respond to perturbations δt_i 
away from the democratic point? 
"""

import numpy as np
import time
import sys

# ============================================================
# K₄ BAND STRUCTURE WITH ℤ₃ FLUX
# ============================================================

omega = np.exp(2j * np.pi / 3)
t_dem = 1.0 / np.sqrt(3)


def cone_energy_grid(Nk, t1=None, t2=None, t3=None):
    """
    E₁(k) on BZ grid with ℤ₃ flux.
    d₁(k) = t₁ + t₂·ω·e^{ik₁} + t₃·ω²·e^{ik₂}
    """
    if t1 is None: t1 = t_dem
    if t2 is None: t2 = t_dem
    if t3 is None: t3 = t_dem
    
    k1 = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
    k2 = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
    K1, K2 = np.meshgrid(k1, k2, indexing='ij')
    
    d1 = t1 + t2 * omega * np.exp(1j * K1) + t3 * omega**2 * np.exp(1j * K2)
    E1 = np.abs(d1)
    
    return E1, K1, K2


# ============================================================
# d_eff COMPUTATION (MUTUAL INFORMATION METHOD)
# ============================================================

def compute_f_positive_E(epsilon, E1_flat, Nk, NE=None):
    """
    f(ε) for E > 0 (arrow of time) via mutual information.
    d_eff = 2 + H(E|k)/H(E)
    f(ε) = 3 - d_eff = 1 - H(E|k)/H(E)
    
    Adapted from d_eff_campaign_v3.py.
    """
    if NE is None:
        NE = max(400, min(2000, int(20 * epsilon + 200)))
    
    N_k = len(E1_flat)
    dk = (2 * np.pi)**2 / N_k
    
    E_min = -3.0 * epsilon
    E_max = np.max(E1_flat) + 6 * epsilon
    E_grid = np.linspace(E_min, E_max, NE)
    dE = E_grid[1] - E_grid[0]
    
    chunk_size = min(N_k, 10000)
    p_E_marginal = np.zeros(NE)
    H_Ek_sum = 0.0
    
    for start in range(0, N_k, chunk_size):
        end = min(start + chunk_size, N_k)
        E1_chunk = E1_flat[start:end]
        
        diff = E_grid[np.newaxis, :] - E1_chunk[:, np.newaxis]
        A = np.exp(-0.5 * (diff / epsilon)**2)
        A[:, E_grid < 0] = 0.0
        
        row_sums = np.sum(A, axis=1) * dE
        valid = row_sums > 1e-300
        
        if np.any(valid):
            p_Ek = np.zeros_like(A)
            p_Ek[valid] = A[valid] / row_sums[valid, np.newaxis]
            
            log_p = np.zeros_like(p_Ek)
            mask = p_Ek > 1e-300
            log_p[mask] = np.log(p_Ek[mask])
            
            H_k = -np.sum(p_Ek * log_p * dE, axis=1)
            H_Ek_sum += np.sum(H_k[valid]) * dk
            
            p_E_marginal += np.sum(A[valid] * dk, axis=0)
    
    p_E_sum = np.sum(p_E_marginal) * dE
    if p_E_sum > 1e-300:
        p_E_marginal /= p_E_sum
    
    mask = p_E_marginal > 1e-300
    H_E = -np.sum(p_E_marginal[mask] * np.log(p_E_marginal[mask]) * dE)
    
    BZ_area = (2 * np.pi)**2
    H_Ek = H_Ek_sum / BZ_area
    
    if H_E > 1e-300:
        d_eff = 2.0 + H_Ek / H_E
        f_eps = 3.0 - d_eff
    else:
        d_eff = 2.0
        f_eps = 1.0
    
    return d_eff, f_eps


# ============================================================
# PART 1: VALIDATION
# ============================================================

def validate(Nk=400):
    print("=" * 72)
    print("  PART 1: VALIDATION AT DEMOCRATIC POINT")
    print("=" * 72)
    
    E1, _, _ = cone_energy_grid(Nk)
    E1_flat = E1.flatten()
    
    print(f"\n  K₄ band structure (Nk={Nk}):")
    print(f"    ⟨E₁⟩   = {np.mean(E1_flat):.6f}")
    print(f"    Var[E₁] = {np.var(E1_flat):.6f}  (expected 0.1736)")
    print(f"    E_min   = {np.min(E1_flat):.6f}  (should be ~0)")
    print(f"    E_max   = {np.max(E1_flat):.6f}")
    
    H0 = 0.5 * (1 + np.log(np.pi/2))
    
    print(f"\n  {'ε':>8s} {'d_eff':>12s} {'f(ε)':>14s} {'ε²·f·(lnε+H₀)':>18s}")
    print("  " + "-" * 56)
    
    for eps in [5, 10, 20, 50, 100, 200, 500]:
        d_eff, f_eps = compute_f_positive_E(eps, E1_flat, Nk)
        if f_eps > 0:
            stab = eps**2 * f_eps * (np.log(eps) + H0)
            print(f"  {eps:8.0f} {d_eff:12.8f} {f_eps:14.8e} {stab:18.6f}")
        else:
            print(f"  {eps:8.0f} {d_eff:12.8f} {f_eps:14.8e} {'---':>18s}")
    
    print(f"\n  Expected stabilization A ≈ 0.0867")
    return E1_flat


# ============================================================
# PART 2: PERTURBATION RESPONSE (ANISOTROPIC)
# ============================================================

def perturbation_response(Nk=400):
    """
    Compute f(ε, δ) for anisotropic perturbation breaking ℤ₃.
    t₁ = t_dem(1+δ), t₂ = t_dem(1-δ/2), t₃ = t_dem(1-δ/2)
    
    This preserves the trace (total bandwidth) while breaking symmetry.
    Physically: tidal gravitational field.
    """
    print(f"\n{'='*72}")
    print("  PART 2: ANISOTROPIC PERTURBATION RESPONSE")
    print(f"{'='*72}")
    print(f"  t_i = t_dem × (1 + δ_i), traceless: δ₁ = δ, δ₂ = δ₃ = -δ/2")
    
    eps_values = [10, 50, 100, 200]
    delta_values = [0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.8]
    
    results = {}
    
    for eps in eps_values:
        print(f"\n  ε = {eps}:")
        print(f"  {'δ':>8s} {'f(ε,δ)':>14s} {'Δf/f₀':>12s} {'c₃=Δf/(f₀δ²)':>14s} {'E_min':>8s}")
        print("  " + "-" * 60)
        
        f0 = None
        f_list = []
        
        for delta in delta_values:
            t1 = t_dem * (1 + delta)
            t2 = t_dem * (1 - delta/2)
            t3 = t_dem * (1 - delta/2)
            
            E1, _, _ = cone_energy_grid(Nk, t1, t2, t3)
            E1_flat = E1.flatten()
            E_min = np.min(E1_flat)
            
            d_eff, f_eps = compute_f_positive_E(eps, E1_flat, Nk)
            f_list.append(f_eps)
            
            if delta == 0:
                f0 = f_eps
                print(f"  {delta:8.4f} {f_eps:14.8e} {'---':>12s} {'---':>14s} {E_min:8.5f}")
            elif f0 is not None and f0 > 1e-20 and f_eps > 1e-20:
                df_rel = (f_eps - f0) / f0
                c3 = df_rel / delta**2 if delta > 0 else 0
                print(f"  {delta:8.4f} {f_eps:14.8e} {df_rel:12.4e} {c3:14.4e} {E_min:8.5f}")
            elif f0 is not None:
                print(f"  {delta:8.4f} {f_eps:14.8e} {'---':>12s} {'---':>14s} {E_min:8.5f}")
        
        results[eps] = {'deltas': delta_values, 'f_values': f_list, 'f0': f0}
    
    return results


# ============================================================
# PART 3: GAP OPENING (ℤ₃-BREAKING MASS TERM)
# ============================================================

def gap_scan(Nk=400):
    """
    The democratic point is gapless (Dirac cone touches E=0).
    Breaking the triangle inequality opens a gap.
    
    If t₁ > t₂ + t₃, there's a gap Δ = t₁ - t₂ - t₃.
    Parameterize: t₁ = t_dem + Δ, t₂ = t₃ = t_dem.
    Gap opens when Δ > 0 (barely, since triangle condition is 
    t₁ < t₂ + t₃ → t_dem + Δ < 2·t_dem → Δ < t_dem ≈ 0.577)
    
    More precisely, the gap at K-point for general t_i with ℤ₃ phases is:
    E_gap = min_k |d₁(k)| which depends on detailed geometry.
    """
    print(f"\n{'='*72}")
    print("  PART 3: GAP OPENING (ℤ₃-BREAKING)")
    print(f"{'='*72}")
    print(f"  t₁ = t_dem + Δ/√3, t₂ = t₃ = t_dem (gap opens at Δ > 0)")
    
    eps_values = [10, 50, 100]
    # delta in units of t_dem, so the physical gap is delta * t_dem
    delta_values = [0, 0.01, 0.03, 0.1, 0.2, 0.3, 0.5, 0.8]
    
    for eps in eps_values:
        print(f"\n  ε = {eps}:")
        print(f"  {'Δ':>8s} {'E_gap':>10s} {'f(ε,Δ)':>14s} {'Δf/f₀':>12s}")
        print("  " + "-" * 48)
        
        f0 = None
        
        for delta in delta_values:
            t1 = t_dem + delta / np.sqrt(3)
            t2 = t_dem
            t3 = t_dem
            
            E1, _, _ = cone_energy_grid(Nk, t1, t2, t3)
            E1_flat = E1.flatten()
            E_gap = np.min(E1_flat)
            
            d_eff, f_eps = compute_f_positive_E(eps, E1_flat, Nk)
            
            if delta == 0:
                f0 = f_eps
                print(f"  {delta:8.4f} {E_gap:10.6f} {f_eps:14.8e} {'---':>12s}")
            elif f0 is not None and f0 > 1e-20 and f_eps > 1e-20:
                df_rel = (f_eps - f0) / f0
                print(f"  {delta:8.4f} {E_gap:10.6f} {f_eps:14.8e} {df_rel:12.4e}")
            elif f0 is not None:
                print(f"  {delta:8.4f} {E_gap:10.6f} {f_eps:14.8e} {'---':>12s}")


# ============================================================
# PART 4: SCALE DEPENDENCE OF RESPONSE
# ============================================================

def scale_dependence(Nk=400):
    """
    How does c₃(ε) depend on ε? If it grows, the mechanism might
    amplify at cosmological scales.
    """
    print(f"\n{'='*72}")
    print("  PART 4: c₃(ε) SCALE DEPENDENCE")
    print(f"{'='*72}")
    
    delta_test = 0.01  # Small traceless perturbation
    
    eps_range = [5, 10, 20, 50, 100, 200, 500]
    
    E1_dem, _, _ = cone_energy_grid(Nk)
    E1_dem_flat = E1_dem.flatten()
    
    t1_p = t_dem * (1 + delta_test)
    t2_p = t_dem * (1 - delta_test/2)
    t3_p = t_dem * (1 - delta_test/2)
    E1_pert, _, _ = cone_energy_grid(Nk, t1_p, t2_p, t3_p)
    E1_pert_flat = E1_pert.flatten()
    
    print(f"\n  δ = {delta_test} (traceless anisotropic)")
    print(f"\n  {'ε':>8s} {'f₀':>14s} {'f_pert':>14s} {'Δf/f₀':>12s} {'c₃':>14s}")
    print("  " + "-" * 66)
    
    c3_list = []
    eps_list = []
    
    for eps in eps_range:
        _, f0 = compute_f_positive_E(eps, E1_dem_flat, Nk)
        _, f_pert = compute_f_positive_E(eps, E1_pert_flat, Nk)
        
        if f0 > 1e-20 and f_pert > 1e-20:
            df_rel = (f_pert - f0) / f0
            c3 = df_rel / delta_test**2
            c3_list.append(c3)
            eps_list.append(eps)
            print(f"  {eps:8.0f} {f0:14.8e} {f_pert:14.8e} {df_rel:12.4e} {c3:14.4e}")
        elif f0 > 1e-20:
            print(f"  {eps:8.0f} {f0:14.8e} {f_pert:14.8e} {'(noise)':>12s}")
        else:
            print(f"  {eps:8.0f} {'(noise)':>14s}")
    
    # Fit power law
    if len(c3_list) >= 3:
        le = np.log(np.array(eps_list))
        lc = np.log(np.abs(np.array(c3_list)))
        coeffs = np.polyfit(le, lc, 1)
        alpha = coeffs[0]
        A = np.exp(coeffs[1])
        
        print(f"\n  Power-law fit: c₃(ε) ≈ {A:.4e} × ε^{alpha:.3f}")
        
        # Extrapolation
        c3_gal = A * (1e57)**alpha
        c3_hub = A * (1e61)**alpha
        print(f"  Extrapolated to ε=10⁵⁷ (galactic): c₃ ≈ {c3_gal:.2e}")
        print(f"  Extrapolated to ε=10⁶¹ (Hubble):   c₃ ≈ {c3_hub:.2e}")
        
        # DM test: need c₃ × (Φ/c²)² ~ O(1)
        # Φ/c² ~ 10⁻⁶ for galaxies, so (Φ/c²)² ~ 10⁻¹²
        dm_effect_gal = c3_gal * 1e-12
        print(f"\n  DM effect at galaxy scale: c₃(10⁵⁷) × (10⁻⁶)² = {dm_effect_gal:.2e}")
        print(f"  Needed for DM: ~1")
        
        return c3_list, eps_list, alpha
    
    return c3_list, eps_list, 0


# ============================================================
# PART 5: DIRECTION-DEPENDENT RESPONSE
# ============================================================

def direction_dependence(Nk=400, eps=50):
    """
    The anisotropic perturbation has a direction. Does f depend on
    WHICH hopping we perturb, or only on the magnitude?
    
    By ℤ₃ symmetry of the democratic point, the three perturbation
    directions should give the same f. Verify this.
    """
    print(f"\n{'='*72}")
    print("  PART 5: DIRECTION DEPENDENCE (ℤ₃ SYMMETRY CHECK)")
    print(f"{'='*72}")
    
    delta = 0.1
    
    configs = [
        ("δ₁ = +δ, δ₂ = δ₃ = -δ/2", t_dem*(1+delta), t_dem*(1-delta/2), t_dem*(1-delta/2)),
        ("δ₂ = +δ, δ₁ = δ₃ = -δ/2", t_dem*(1-delta/2), t_dem*(1+delta), t_dem*(1-delta/2)),
        ("δ₃ = +δ, δ₁ = δ₂ = -δ/2", t_dem*(1-delta/2), t_dem*(1-delta/2), t_dem*(1+delta)),
        ("δ₁ = +δ, δ₂ = -δ, δ₃ = 0", t_dem*(1+delta), t_dem*(1-delta), t_dem),
        ("δ₁ = +δ, δ₂ = 0, δ₃ = -δ", t_dem*(1+delta), t_dem, t_dem*(1-delta)),
    ]
    
    print(f"\n  ε = {eps}, δ = {delta}")
    
    # Democratic baseline
    E1_dem, _, _ = cone_energy_grid(Nk)
    _, f0 = compute_f_positive_E(eps, E1_dem.flatten(), Nk)
    print(f"  Democratic: f₀ = {f0:.8e}")
    
    print(f"\n  {'Config':<38s} {'f(ε,δ)':>14s} {'Δf/f₀':>12s} {'E_min':>10s}")
    print("  " + "-" * 76)
    
    for label, t1, t2, t3 in configs:
        E1, _, _ = cone_energy_grid(Nk, t1, t2, t3)
        E1_flat = E1.flatten()
        _, f_eps = compute_f_positive_E(eps, E1_flat, Nk)
        E_min = np.min(E1_flat)
        
        if f0 > 1e-20 and f_eps > 1e-20:
            df_rel = (f_eps - f0) / f0
            print(f"  {label:<38s} {f_eps:14.8e} {df_rel:12.4e} {E_min:10.6f}")
        else:
            print(f"  {label:<38s} {f_eps:14.8e} {'---':>12s} {E_min:10.6f}")


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 72)
    print("  GEOMETRIC DARK MATTER PROBE v2 (CORRECTED)")
    print("  d_eff response to boundary perturbations")
    print("  K₄ Dirac cone with ℤ₃ flux, mutual information method")
    print("=" * 72)
    
    Nk = 400
    if "--hires" in sys.argv: Nk = 600
    if "--quick" in sys.argv: Nk = 200
    
    t0 = time.time()
    
    # Part 1: Validate
    E1_flat = validate(Nk)
    
    # Part 2: Anisotropic response
    results = perturbation_response(Nk)
    
    # Part 3: Gap opening
    gap_scan(Nk)
    
    # Part 4: Scale dependence
    c3_list, eps_list, alpha = scale_dependence(Nk)
    
    # Part 5: Direction check
    direction_dependence(Nk)
    
    # ================================================================
    # VERDICT
    # ================================================================
    dt = time.time() - t0
    
    print(f"\n{'='*72}")
    print(f"  VERDICT (computed in {dt:.0f}s)")
    print(f"{'='*72}")
    
    if len(c3_list) >= 3 and alpha != 0:
        A_fit = np.exp(np.polyfit(np.log(eps_list), np.log(np.abs(c3_list)), 1)[1])
        c3_gal = A_fit * (1e57)**alpha
        dm_effect = c3_gal * 1e-12  # times (Φ/c²)²
        
        if abs(dm_effect) > 0.01:
            print(f"\n  INTERESTING: c₃ scaling (α = {alpha:.2f}) produces non-negligible")
            print(f"  effect at galactic scales: Δf/f₀ ~ {dm_effect:.2e}")
            print(f"  STATUS: ALIVE — further investigation needed")
        elif abs(dm_effect) > 1e-10:
            print(f"\n  WEAK: c₃ scaling gives {dm_effect:.2e} at galactic scales.")
            print(f"  Many orders too small for dark matter.")
            print(f"  STATUS: SOFT KILL — mechanism present but too weak")
        else:
            print(f"\n  NEGLIGIBLE: c₃ scaling gives {dm_effect:.2e} at galactic scales.")
            print(f"  STATUS: KILL — geometric DM via boundary perturbation doesn't work")
    else:
        print(f"\n  INSUFFICIENT DATA for verdict.")
        print(f"  c₃ values: {c3_list}")
    
    print(f"\n{'='*72}")


if __name__ == "__main__":
    main()
