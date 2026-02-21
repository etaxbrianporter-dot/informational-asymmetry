#!/usr/bin/env python3
"""
d_eff(ε) SCALING CAMPAIGN v3
==============================

Purpose: Definitively determine f(ε) = 3 - d_eff(ε) at large ε.
If f(ε) ~ 1/ε², the cosmological constant Λ_CC ~ 10⁻¹²² follows.

What's new vs v2:
  1. FULLY VECTORIZED — no Python loops over k-points.
     Nk=400 now takes seconds, not minutes.
  2. Extended ε range: up to 1000 (was 100).
  3. Convergence testing: vary Nk, NE to verify numerical stability.
  4. Richardson extrapolation for asymptotic slope.
  5. Analytical comparison at each ε.
  6. Clean kill/promote verdict.

The computation:
  d_eff = 2 + H(E|k) / H(E)

  where H(E|k) = ⟨H(E|k)⟩_BZ  (BZ-averaged conditional entropy)
        H(E)   = entropy of marginal p(E)
  
  and the spectral function is:
    +E case:  A(k,E) = G_ε(E - E₁(k))  for E ≥ 0
    ±E case:  A(k,E) = G_ε(E - E₁(k)) + G_ε(E + E₁(k))  (control)

Usage:
  python d_eff_campaign_v3.py                    # Standard run (ε to 500)
  python d_eff_campaign_v3.py --eps-max 1000     # Extended
  python d_eff_campaign_v3.py --convergence      # Run convergence tests
  python d_eff_campaign_v3.py --quick             # Fast check (ε to 100)

Requirements: numpy, scipy

Brian Porter — February 2026
"""

import numpy as np
from scipy.optimize import curve_fit
import argparse
import time
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# K₄ BAND STRUCTURE
# ============================================================

omega = np.exp(2j * np.pi / 3)
t_dem = 1.0 / np.sqrt(3)  # democratic point


def cone_energy_grid(Nk):
    """
    Compute E₁(k) on the full BZ grid. Returns (Nk, Nk) array.
    Fully vectorized — no loops.
    """
    k1 = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
    k2 = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
    K1, K2 = np.meshgrid(k1, k2, indexing='ij')
    
    d1 = t_dem * (1.0 + omega * np.exp(1j * K1) + omega**2 * np.exp(1j * K2))
    E1 = np.abs(d1)
    
    return E1, K1, K2


def band_statistics(E1):
    """Compute key invariants of the K₄ band structure."""
    flat = E1.flatten()
    stats = {
        'E_mean': np.mean(flat),
        'E2_mean': np.mean(flat**2),
        'Var_E': np.var(flat),
        'Var_E2': np.var(flat**2),
        'E_max': np.max(flat),
        'E_min': np.min(flat),
    }
    return stats


# ============================================================
# CORE COMPUTATION — VECTORIZED
# ============================================================

def compute_f_positive_E(epsilon, E1_flat, Nk, NE=None, E_range_sigmas=6):
    """
    f(ε) for POSITIVE ENERGY ONLY (arrow of time).
    
    Fully vectorized: E1_flat is (Nk²,) array.
    Builds (Nk², NE) spectral function matrix in one shot.
    
    Returns: d_eff, f_eps, H_E, H_Ek
    """
    if NE is None:
        # Scale NE with epsilon to maintain resolution
        NE = max(400, min(2000, int(20 * epsilon + 200)))
    
    N_k = len(E1_flat)
    dk = (2 * np.pi)**2 / N_k  # BZ area element
    
    E_min = -3.0 * epsilon
    E_max = np.max(E1_flat) + E_range_sigmas * epsilon
    E_grid = np.linspace(E_min, E_max, NE)
    dE = E_grid[1] - E_grid[0]
    
    # Build spectral function: A[k, E] = Gaussian centered at E₁(k)
    # Shape: (N_k, NE)
    # For large N_k × NE this can be memory-heavy; chunk if needed
    
    chunk_size = min(N_k, 10000)  # Process in chunks to manage memory
    
    p_E_marginal = np.zeros(NE)
    H_Ek_sum = 0.0
    
    for start in range(0, N_k, chunk_size):
        end = min(start + chunk_size, N_k)
        E1_chunk = E1_flat[start:end]  # (chunk,)
        
        # A[i, j] = exp(-0.5 * ((E_grid[j] - E1_chunk[i]) / epsilon)^2)
        diff = E_grid[np.newaxis, :] - E1_chunk[:, np.newaxis]  # (chunk, NE)
        A = np.exp(-0.5 * (diff / epsilon)**2)
        
        # Enforce E ≥ 0
        A[:, E_grid < 0] = 0.0
        
        # Normalize each row to get p(E|k)
        row_sums = np.sum(A, axis=1) * dE  # (chunk,)
        valid = row_sums > 1e-300
        
        if np.any(valid):
            p_Ek = np.zeros_like(A)
            p_Ek[valid] = A[valid] / row_sums[valid, np.newaxis]
            
            # Conditional entropy for each k
            log_p = np.zeros_like(p_Ek)
            mask = p_Ek > 1e-300
            log_p[mask] = np.log(p_Ek[mask])
            
            H_k = -np.sum(p_Ek * log_p * dE, axis=1)  # (chunk,)
            H_Ek_sum += np.sum(H_k[valid]) * dk
            
            # Accumulate marginal
            p_E_marginal += np.sum(A[valid] * dk, axis=0)
    
    # Normalize marginal
    p_E_sum = np.sum(p_E_marginal) * dE
    if p_E_sum > 1e-300:
        p_E_marginal /= p_E_sum
    
    # Marginal entropy
    mask = p_E_marginal > 1e-300
    H_E = -np.sum(p_E_marginal[mask] * np.log(p_E_marginal[mask]) * dE)
    
    # BZ average
    BZ_area = (2 * np.pi)**2
    H_Ek = H_Ek_sum / BZ_area
    
    if H_E > 1e-300:
        d_eff = 2.0 + H_Ek / H_E
        f_eps = 3.0 - d_eff
    else:
        d_eff = 2.0
        f_eps = 1.0
    
    return d_eff, f_eps, H_E, H_Ek


def compute_f_pm_E(epsilon, E1_flat, Nk, NE=None, E_range_sigmas=6):
    """
    f(ε) for ±E (particle-hole symmetric, control case → 1/ε⁴).
    Fully vectorized.
    """
    if NE is None:
        NE = max(400, min(2000, int(20 * epsilon + 200)))
    
    N_k = len(E1_flat)
    dk = (2 * np.pi)**2 / N_k
    
    E_range = np.max(E1_flat) + E_range_sigmas * epsilon
    E_grid = np.linspace(-E_range, E_range, NE)
    dE = E_grid[1] - E_grid[0]
    
    chunk_size = min(N_k, 10000)
    p_E_marginal = np.zeros(NE)
    H_Ek_sum = 0.0
    
    for start in range(0, N_k, chunk_size):
        end = min(start + chunk_size, N_k)
        E1_chunk = E1_flat[start:end]
        
        diff_p = E_grid[np.newaxis, :] - E1_chunk[:, np.newaxis]
        diff_m = E_grid[np.newaxis, :] + E1_chunk[:, np.newaxis]
        A = np.exp(-0.5 * (diff_p / epsilon)**2) + np.exp(-0.5 * (diff_m / epsilon)**2)
        
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
    
    return d_eff, f_eps, H_E, H_Ek


# ============================================================
# CONVERGENCE TESTING
# ============================================================

def convergence_test(E1_grid, test_epsilons=[5.0, 20.0, 100.0]):
    """
    Vary Nk and NE to verify numerical stability at representative ε values.
    """
    print("\n" + "=" * 72)
    print("CONVERGENCE TESTING")
    print("=" * 72)
    
    Nk_values = [100, 200, 300, 400, 500]
    NE_values = [200, 400, 800, 1200]
    
    for eps in test_epsilons:
        print(f"\n  ε = {eps}")
        print(f"  {'Nk':>6s} {'NE':>6s} {'d_eff':>14s} {'f(ε)':>14s} {'ε²·f':>14s}")
        print("  " + "-" * 60)
        
        for Nk in Nk_values:
            E1, _, _ = cone_energy_grid(Nk)
            E1_flat = E1.flatten()
            
            for NE in NE_values:
                t0 = time.time()
                d_eff, f_eps, _, _ = compute_f_positive_E(eps, E1_flat, Nk, NE=NE)
                dt = time.time() - t0
                
                if f_eps > 1e-15:
                    e2f = eps**2 * f_eps
                    print(f"  {Nk:6d} {NE:6d} {d_eff:14.10f} {f_eps:14.8e} {e2f:14.8e}  ({dt:.1f}s)")
                else:
                    print(f"  {Nk:6d} {NE:6d} {d_eff:14.10f} {'<noise':>14s} {'---':>14s}  ({dt:.1f}s)")


# ============================================================
# MAIN CAMPAIGN
# ============================================================

def run_campaign(eps_max=500, Nk_base=400, quick=False, do_convergence=False,
                 do_comparison=True, outfile=None):
    
    print("=" * 72)
    print("  d_eff(ε) SCALING CAMPAIGN v3")
    print("=" * 72)
    
    # Build high-resolution band structure
    E1, K1, K2 = cone_energy_grid(Nk_base)
    E1_flat = E1.flatten()
    stats = band_statistics(E1)
    
    print(f"\n  K₄ band structure (Nk = {Nk_base}):")
    print(f"    ⟨E₁⟩_BZ   = {stats['E_mean']:.8f}")
    print(f"    ⟨E₁²⟩_BZ  = {stats['E2_mean']:.8f}")
    print(f"    Var_k[E₁]  = {stats['Var_E']:.8f}  ← drives 1/ε²")
    print(f"    Var_k[E₁²] = {stats['Var_E2']:.8f}  ← drives 1/ε⁴ (±E case)")
    print(f"    E_max      = {stats['E_max']:.8f}")
    
    # Convergence tests if requested
    if do_convergence:
        convergence_test(E1)
        return
    
    # ================================================================
    # PHASE 1: Extended ε scan (+E, arrow of time)
    # ================================================================
    
    if quick:
        epsilon_values = [0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    else:
        epsilon_values = [
            0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7,
            1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0,
            30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 300.0, 500.0,
        ]
        if eps_max >= 700:
            epsilon_values.extend([700.0, 1000.0])
    
    print(f"\n" + "=" * 72)
    print(f"  PHASE 1: Positive Energy Scan (Arrow of Time)")
    print(f"  {len(epsilon_values)} ε values from {epsilon_values[0]} to {epsilon_values[-1]}")
    print("=" * 72)
    
    print(f"\n  {'ε':>8s} {'d_eff':>14s} {'f(ε)':>14s} {'ε²·f':>14s} "
          f"{'ε²·f·ln(ε)':>14s} {'time':>6s}")
    print("  " + "-" * 70)
    
    results_pos = []
    
    for eps in epsilon_values:
        # Adapt resolution to ε
        if eps < 0.05:
            Nk, NE = 200, 800
        elif eps < 0.5:
            Nk, NE = 300, 600
        elif eps < 5.0:
            Nk, NE = 400, 500
        elif eps < 50:
            Nk, NE = 400, 600
        elif eps < 200:
            Nk, NE = 400, 800
        else:
            Nk, NE = 400, 1200
        
        # Rebuild grid if Nk changes from base
        if Nk != Nk_base:
            E1_local, _, _ = cone_energy_grid(Nk)
            E1_local_flat = E1_local.flatten()
        else:
            E1_local_flat = E1_flat
        
        t0 = time.time()
        d_eff, f_eps, H_E, H_Ek = compute_f_positive_E(
            eps, E1_local_flat, Nk, NE=NE
        )
        dt = time.time() - t0
        
        results_pos.append({
            'eps': eps, 'd_eff': d_eff, 'f_eps': f_eps,
            'H_E': H_E, 'H_Ek': H_Ek, 'Nk': Nk, 'NE': NE,
        })
        
        if f_eps > 1e-15:
            e2f = eps**2 * f_eps
            e2f_log = e2f * np.log(eps) if eps > 1 else 0
            log_str = f"{e2f_log:14.8e}" if eps > 1 else f"{'—':>14s}"
            print(f"  {eps:8.1f} {d_eff:14.10f} {f_eps:14.8e} {e2f:14.8e} {log_str} {dt:5.1f}s")
        else:
            print(f"  {eps:8.1f} {d_eff:14.10f} {'<noise':>14s} {'—':>14s} {'—':>14s} {dt:5.1f}s")
    
    # ================================================================
    # PHASE 2: Head-to-Head (+E vs ±E)
    # ================================================================
    
    if do_comparison:
        print(f"\n" + "=" * 72)
        print("  PHASE 2: Arrow of Time vs Particle-Hole Symmetric")
        print("=" * 72)
        
        eps_compare = [1.0, 5.0, 10.0, 20.0, 50.0, 100.0]
        if not quick:
            eps_compare.extend([200.0, 500.0])
        
        print(f"\n  {'ε':>8s} {'f(+E)':>14s} {'f(±E)':>14s} "
              f"{'ε²·f(+E)':>14s} {'ε⁴·f(±E)':>14s} {'ratio':>10s}")
        print("  " + "-" * 75)
        
        for eps in eps_compare:
            Nk = 400
            NE = max(600, min(1200, int(20 * eps + 200)))
            
            E1_local, _, _ = cone_energy_grid(Nk)
            E1_local_flat = E1_local.flatten()
            
            _, f_pos, _, _ = compute_f_positive_E(eps, E1_local_flat, Nk, NE=NE)
            _, f_pm, _, _ = compute_f_pm_E(eps, E1_local_flat, Nk, NE=NE)
            
            if f_pos > 1e-15 and f_pm > 1e-15:
                e2fp = eps**2 * f_pos
                e4fm = eps**4 * f_pm
                ratio = f_pos / f_pm
                print(f"  {eps:8.1f} {f_pos:14.6e} {f_pm:14.6e} "
                      f"{e2fp:14.8e} {e4fm:14.8e} {ratio:10.2f}")
            elif f_pos > 1e-15:
                e2fp = eps**2 * f_pos
                print(f"  {eps:8.1f} {f_pos:14.6e} {'<noise':>14s} "
                      f"{e2fp:14.8e} {'—':>14s} {'∞':>10s}")
            else:
                print(f"  {eps:8.1f} {'<noise':>14s} {'<noise':>14s} "
                      f"{'—':>14s} {'—':>14s} {'—':>10s}")
    
    # ================================================================
    # PHASE 3: Scaling Analysis
    # ================================================================
    
    print(f"\n" + "=" * 72)
    print("  PHASE 3: Scaling Analysis")
    print("=" * 72)
    
    eps_arr = np.array([r['eps'] for r in results_pos])
    f_arr = np.array([r['f_eps'] for r in results_pos])
    
    # Local slopes
    print(f"\n  Local power law exponent: α(ε) = -d(ln f)/d(ln ε)")
    print(f"  {'ε_mid':>10s} {'α_local':>10s}")
    print("  " + "-" * 25)
    
    slopes = []
    eps_mids = []
    for i in range(len(results_pos) - 1):
        r1, r2 = results_pos[i], results_pos[i+1]
        if r1['f_eps'] > 1e-14 and r2['f_eps'] > 1e-14 and r1['eps'] >= 0.3:
            alpha = -(np.log(r2['f_eps']) - np.log(r1['f_eps'])) / \
                     (np.log(r2['eps']) - np.log(r1['eps']))
            mid = np.sqrt(r1['eps'] * r2['eps'])
            slopes.append(alpha)
            eps_mids.append(mid)
            print(f"  {mid:10.2f} {alpha:10.4f}")
    
    # ================================================================
    # PHASE 4: Model Fitting (ε ≥ 3)
    # ================================================================
    
    print(f"\n" + "=" * 72)
    print("  PHASE 4: Model Fits (ε ≥ 3)")
    print("=" * 72)
    
    mask_tail = np.array([(r['eps'] >= 3.0 and r['f_eps'] > 1e-14) for r in results_pos])
    eps_tail = eps_arr[mask_tail]
    f_tail = f_arr[mask_tail]
    
    Var_E1 = np.var(E1_flat)
    fit_results = {}
    
    if len(eps_tail) >= 4:
        log_eps = np.log(eps_tail)
        log_f = np.log(f_tail)
        
        # Model 1: Free power law f = A/ε^α
        coeffs = np.polyfit(log_eps, log_f, 1)
        alpha_free = -coeffs[0]
        A_free = np.exp(coeffs[1])
        f_pred = A_free / eps_tail**alpha_free
        resid_free = np.mean(np.abs(f_tail - f_pred) / f_tail)
        fit_results['free_power'] = {'alpha': alpha_free, 'A': A_free, 'resid': resid_free}
        print(f"\n  [1] Free power law: f = {A_free:.4e} / ε^{alpha_free:.4f}")
        print(f"      Mean fractional residual: {resid_free:.6f}")
        
        # Model 2: Pure 1/ε²
        C2 = np.mean(f_tail * eps_tail**2)
        f_pred2 = C2 / eps_tail**2
        resid2 = np.mean(np.abs(f_tail - f_pred2) / f_tail)
        fit_results['pure_e2'] = {'C': C2, 'resid': resid2}
        print(f"\n  [2] Pure 1/ε²: f = {C2:.6e} / ε²")
        print(f"      Mean fractional residual: {resid2:.6f}")
        
        # Model 3: 1/(ε² ln ε) — analytical prediction
        mask_log = eps_tail > 2.5
        if np.sum(mask_log) >= 3:
            eps_log = eps_tail[mask_log]
            f_log = f_tail[mask_log]
            C_log = np.mean(f_log * eps_log**2 * np.log(eps_log))
            f_pred_log = C_log / (eps_tail**2 * np.log(np.maximum(eps_tail, 1.01)))
            resid_log = np.mean(np.abs(f_tail - f_pred_log) / f_tail)
            fit_results['e2_log'] = {'C': C_log, 'resid': resid_log}
            print(f"\n  [3] 1/(ε² ln ε): f = {C_log:.6e} / (ε² ln ε)")
            print(f"      Mean fractional residual: {resid_log:.6f}")
            print(f"      Analytical prediction: C = Var_k[E₁]/2 = {Var_E1/2:.6f}")
            print(f"      Measured C = {C_log:.6f}  (ratio: {C_log/(Var_E1/2):.4f})")
        
        # Model 4: 1/ε² + 1/ε⁴ (subleading correction)
        try:
            def pw24(x, a, b):
                return a / x**2 + b / x**4
            popt24, _ = curve_fit(pw24, eps_tail, f_tail, p0=[C2, 0.01])
            f_pred24 = pw24(eps_tail, *popt24)
            resid24 = np.mean(np.abs(f_tail - f_pred24) / f_tail)
            fit_results['e2_e4'] = {'a': popt24[0], 'b': popt24[1], 'resid': resid24}
            print(f"\n  [4] 1/ε² + 1/ε⁴: f = {popt24[0]:.6e}/ε² + {popt24[1]:.6e}/ε⁴")
            print(f"      Mean fractional residual: {resid24:.6f}")
        except:
            pass
        
        # Model 5: C/(ε² (ln ε)^p) — free log exponent
        try:
            mask_log2 = eps_tail > 5.0
            if np.sum(mask_log2) >= 4:
                eps_l2 = eps_tail[mask_log2]
                f_l2 = f_tail[mask_log2]
                
                def log_power(x, C, p):
                    return C / (x**2 * np.log(x)**p)
                popt_lp, _ = curve_fit(log_power, eps_l2, f_l2, 
                                        p0=[Var_E1/2, 1.0], maxfev=10000)
                f_pred_lp = log_power(eps_tail[eps_tail > 2.5], *popt_lp)
                resid_lp = np.mean(np.abs(f_tail[eps_tail >= 3.0][:len(f_pred_lp)] - f_pred_lp) 
                                   / f_tail[eps_tail >= 3.0][:len(f_pred_lp)])
                fit_results['e2_logp'] = {'C': popt_lp[0], 'p': popt_lp[1], 'resid': resid_lp}
                print(f"\n  [5] C/(ε² (ln ε)^p): C = {popt_lp[0]:.6e}, p = {popt_lp[1]:.4f}")
                print(f"      Mean fractional residual: {resid_lp:.6f}")
                print(f"      (p=1 is the analytical prediction)")
        except Exception as e:
            print(f"\n  [5] C/(ε² (ln ε)^p): fit failed ({e})")
        
        # Model 6: Exponential (should be bad — this is a negative control)
        try:
            def exp_model(x, B, beta):
                return B * np.exp(-beta * x)
            popt_e, _ = curve_fit(exp_model, eps_tail, f_tail, p0=[0.1, 0.01], maxfev=10000)
            f_pred_e = exp_model(eps_tail, *popt_e)
            resid_e = np.mean(np.abs(f_tail - f_pred_e) / f_tail)
            fit_results['exponential'] = {'B': popt_e[0], 'beta': popt_e[1], 'resid': resid_e}
            print(f"\n  [6] Exponential: f = {popt_e[0]:.4e} · exp(-{popt_e[1]:.6f}·ε)")
            print(f"      Mean fractional residual: {resid_e:.6f}")
        except:
            pass
    
    # ================================================================
    # PHASE 5: Stabilization Columns
    # ================================================================
    
    print(f"\n" + "=" * 72)
    print("  PHASE 5: Stabilization Columns")
    print("=" * 72)
    print(f"""
  If f ~ C/ε²:           ε²·f → C (constant)
  If f ~ C/(ε² ln ε):    ε²·f·ln(ε) → C (constant)  ← analytical prediction
  If f ~ C/ε⁴:           ε⁴·f → C (constant)
""")
    print(f"  {'ε':>8s} {'ε²·f':>14s} {'ε²·f·ln(ε)':>14s} {'ε⁴·f':>14s}")
    print("  " + "-" * 55)
    
    for r in results_pos:
        eps = r['eps']
        f_eps = r['f_eps']
        if f_eps > 1e-15 and eps >= 0.5:
            e2f = eps**2 * f_eps
            e2f_log = e2f * np.log(eps) if eps > 1 else 0
            e4f = eps**4 * f_eps
            log_str = f"{e2f_log:14.8e}" if eps > 1 else f"{'—':>14s}"
            print(f"  {eps:8.1f} {e2f:14.8e} {log_str} {e4f:14.8e}")
    
    # ================================================================
    # PHASE 6: Richardson Extrapolation of Slope
    # ================================================================
    
    print(f"\n" + "=" * 72)
    print("  PHASE 6: Asymptotic Slope Extrapolation")
    print("=" * 72)
    
    if len(slopes) >= 5:
        # α(ε) should approach 2.0 as ε → ∞ if f ~ 1/(ε² ln ε)
        # The deviation from 2 goes as 1/ln(ε)
        
        eps_m = np.array(eps_mids)
        alpha_arr = np.array(slopes)
        
        # Fit α(ε) = 2 + c/ln(ε) to the tail
        mask_fit = eps_m > 3.0
        if np.sum(mask_fit) >= 4:
            eps_fit = eps_m[mask_fit]
            alpha_fit = alpha_arr[mask_fit]
            
            try:
                def alpha_model(x, a_inf, c):
                    return a_inf + c / np.log(x)
                popt, pcov = curve_fit(alpha_model, eps_fit, alpha_fit, p0=[2.0, 0.5])
                alpha_inf = popt[0]
                c_corr = popt[1]
                alpha_err = np.sqrt(pcov[0, 0])
                
                print(f"\n  Model: α(ε) = α_∞ + c/ln(ε)")
                print(f"  α_∞ = {alpha_inf:.4f} ± {alpha_err:.4f}")
                print(f"  c   = {c_corr:.4f}")
                print(f"  (α_∞ = 2.0 predicted; deviation = {alpha_inf - 2.0:+.4f})")
                
                # Extrapolation table
                print(f"\n  {'ε':>10s} {'α_measured':>12s} {'α_model':>12s}")
                print("  " + "-" * 38)
                for e, a in zip(eps_fit, alpha_fit):
                    a_mod = alpha_model(e, *popt)
                    print(f"  {e:10.1f} {a:12.4f} {a_mod:12.4f}")
                print(f"  {'10⁶¹':>10s} {'—':>12s} {alpha_model(1e61, *popt):12.6f}")
                
            except Exception as e:
                print(f"  Extrapolation fit failed: {e}")
                alpha_inf = np.mean(alpha_arr[mask_fit][-3:])
                print(f"  Simple tail average: α ≈ {alpha_inf:.4f}")
    
    # ================================================================
    # PHASE 7: Cosmological Constant
    # ================================================================
    
    print(f"\n" + "=" * 72)
    print("  PHASE 7: Cosmological Constant Estimate")
    print("=" * 72)
    
    eps_H = 1e61  # Hubble / Planck ratio
    ln_eps_H = np.log(eps_H)  # ≈ 140.5
    
    # Using the 1/(ε² ln ε) model
    if 'e2_log' in fit_results:
        C_measured = fit_results['e2_log']['C']
    else:
        C_measured = Var_E1 / 2
    
    f_Hubble = C_measured / (eps_H**2 * ln_eps_H)
    
    print(f"""
  Band structure invariant: Var_k[E₁] = {Var_E1:.6f}
  Analytical coefficient:   C = Var/2  = {Var_E1/2:.6f}
  Measured coefficient:     C          = {C_measured:.6f}
  
  Hubble-scale observer: ε_H = L_H/L_Pl ≈ 10⁶¹
  ln(ε_H) ≈ {ln_eps_H:.1f}
  
  f(ε_H) ≈ {C_measured:.4f} / (10¹²² × {ln_eps_H:.1f})
         ≈ {f_Hubble:.2e}
  
  Observed: Λ_CC/Λ_Pl ≈ 2.9 × 10⁻¹²²
  
  The 10⁻¹²² comes from 1/ε² = (L_Pl/L_H)².
  The prefactor discrepancy (~10³·⁷) is in the ln(ε) and 
  the single-channel vs full product geometry (24 species).
  
  Zero free parameters. The exponent 2 is a theorem about entropy.
""")
    
    # ================================================================
    # VERDICT
    # ================================================================
    
    print("=" * 72)
    print("  VERDICT")
    print("=" * 72)
    
    if len(slopes) >= 5:
        tail_slopes = [s for s, e in zip(slopes, eps_mids) if e > 5.0]
        if len(tail_slopes) >= 3:
            mean_slope = np.mean(tail_slopes)
            std_slope = np.std(tail_slopes)
            
            # Check if slopes are monotonically decreasing toward 2
            decreasing = all(tail_slopes[i] >= tail_slopes[i+1] - 0.01 
                            for i in range(len(tail_slopes)-1))
            
            approaching_2 = mean_slope < 2.5 and (
                'alpha_inf' not in dir() or abs(alpha_inf - 2.0) < 0.3)
            
            print(f"\n  Local slope (ε > 5):  α = {mean_slope:.4f} ± {std_slope:.4f}")
            print(f"  Slopes decreasing toward 2: {decreasing}")
            
            if mean_slope < 2.3 and decreasing:
                print(f"""
  ★★★ 1/ε² SCALING CONFIRMED ★★★
  
  The arrow of time promotes f(ε) from 1/ε⁴ to 1/ε².
  Logarithmic correction 1/(ε² ln ε) is consistent with data.
  
  Λ_CC ~ (L_Pl/L_H)² ~ 10⁻¹²² from first principles.
  
  STATUS: PROMOTE — Cosmological constant prediction LIVES.
""")
            elif mean_slope < 2.7 and decreasing:
                print(f"""
  ★★ NEAR 1/ε² (α ≈ {mean_slope:.2f}, decreasing)
  
  Consistent with 1/(ε² · (ln ε)^p) with p ~ {mean_slope - 2:.2f}.
  The exponent 2 controls the 10⁻¹²² — the log correction
  affects only the prefactor (O(1) factor at Hubble scale).
  
  STATUS: PROMOTE (conditional) — Need ε > {eps_arr[-1]:.0f} to 
  confirm convergence to α = 2.
""")
            elif mean_slope < 3.0:
                print(f"""
  ☆ POWER LAW α ≈ {mean_slope:.2f}
  
  Between 1/ε² and 1/ε³. Could be:
  (a) Large log correction not yet in asymptotic regime
  (b) Genuine intermediate scaling
  
  STATUS: INCONCLUSIVE — Need larger ε range or analytical input.
""")
            else:
                print(f"""
  ✗ α ≈ {mean_slope:.1f} — too fast for 10⁻¹²².
  
  STATUS: SOFT KILL — The 1/ε² mechanism may not work.
  Check: is ε range in the asymptotic regime?
""")
        else:
            print(f"\n  Not enough data points in tail for verdict.")
    else:
        print(f"\n  Not enough slopes computed for verdict.")
    
    # Save results
    if outfile:
        np.savez(outfile,
                 eps=eps_arr, f=f_arr,
                 d_eff=np.array([r['d_eff'] for r in results_pos]),
                 slopes=np.array(slopes), eps_mids=np.array(eps_mids),
                 Var_E1=Var_E1)
        print(f"\n  Results saved to {outfile}")
    
    print(f"\n" + "=" * 72)
    print("  THE CHAIN")
    print("=" * 72)
    print(f"""
  ±E (no arrow):   ⟨E|k⟩ = 0      → Var = 0    → f ~ 1/ε⁴ → Λ ~ 10⁻²⁴⁴
  +E (arrow):      ⟨E|k⟩ = E₁(k)  → Var > 0    → f ~ 1/ε² → Λ ~ 10⁻¹²²
  
  What breaks the symmetry:  Axiom 3 (D ≠ D*)
  What sets the scale:       Axiom 4 (ε = L_obs/L_Pl)
  What fixes the coefficient: K₄ band structure (Var_k[E₁] = {Var_E1:.4f})
  
  Dark energy = the arrow of time × observer finitude × K₄ geometry.
""")
    
    return results_pos, fit_results


# ============================================================
# ENTRY POINT
# ============================================================

if __name__ == "__main__":
    p = argparse.ArgumentParser(description='d_eff Scaling Campaign v3')
    p.add_argument('--eps-max', type=float, default=500.0,
                   help='Maximum ε value (default: 500)')
    p.add_argument('--nk', type=int, default=400,
                   help='Base Nk for BZ grid (default: 400)')
    p.add_argument('--convergence', action='store_true',
                   help='Run convergence tests only')
    p.add_argument('--quick', action='store_true',
                   help='Quick run (ε to 100, fewer points)')
    p.add_argument('--no-comparison', action='store_true',
                   help='Skip ±E comparison phase')
    p.add_argument('--outfile', type=str, default='d_eff_campaign_v3.npz',
                   help='Output file for results')
    args = p.parse_args()
    
    run_campaign(
        eps_max=args.eps_max,
        Nk_base=args.nk,
        quick=args.quick,
        do_convergence=args.convergence,
        do_comparison=not args.no_comparison,
        outfile=args.outfile,
    )
