"""
johnson_tools.py — Johnson association scheme for K_{2n} matching overlaps
=========================================================================
All eigenvalues, dimensions, and budget identities in closed form.
"""
import numpy as np


def dfact(k):
    """Double factorial: k!! = k(k-2)(k-4)...1 for odd k, with (-1)!!=1, 1!!=1"""
    if k <= 1:
        return 1
    return k * dfact(k - 2)


def johnson_parameters(n):
    """
    Compute all Johnson scheme parameters for K_{2n}.
    n = number of edges per matching (= half the vertices).
    
    Returns dict with all closed-form quantities.
    """
    N = dfact(2*n - 1)          # (2n-1)!! = total matchings
    E = n * (2*n - 1)            # total edges of K_{2n}
    r = n                        # edges per matching (row sum of A)
    c = dfact(2*n - 3)           # matchings per edge (column sum of A) = N_hub
    
    # Overlap matrix eigenvalues
    lam_max = 2 * n * dfact(2*n - 3)
    lam_mid = 4 * (n - 1) * dfact(2*n - 5) if n >= 3 else (4 if n == 2 else 0)
    
    # Eigenspace dimensions
    d_trivial = 1
    d_kernel_johnson = 2*n - 1   # V₁ dimension in J(2n,2)
    d_phys = n * (2*n - 3)       # V₂ dimension = hook length for (2n-2, 2)
    d_kernel_total = N - d_phys - d_trivial
    
    # Budget identity
    budget_lhs = lam_mid * d_phys
    budget_rhs = 2 * (n - 1) * lam_max
    
    # Hub memory fractions
    hub_trivial = 1.0 / (2*n - 1)
    hub_physical = (2*n - 2.0) / (2*n - 1)
    
    # Physical weight fraction
    phys_weight = (2*n - 2.0) / (2*n - 1)
    
    return {
        'n': n,
        'N': N,
        'E': E,
        'r': r,
        'c': c,
        'lam_max': lam_max,
        'lam_mid': lam_mid,
        'd_trivial': d_trivial,
        'd_phys': d_phys,
        'd_kernel_total': d_kernel_total,
        'd_kernel_johnson': d_kernel_johnson,
        'budget_lhs': budget_lhs,
        'budget_rhs': budget_rhs,
        'hub_trivial_frac': hub_trivial,
        'hub_physical_frac': hub_physical,
        'phys_weight_frac': phys_weight,
        'trace_O': 2 * n * N,
    }


def physical_projector_from_overlap(O, n):
    """
    Compute the exact polynomial projector onto V₂ (physical sector).
    
    Π_phys = O(O - λ_max I) / [λ_mid(λ_mid - λ_max)]
    """
    p = johnson_parameters(n)
    lam_max = p['lam_max']
    lam_mid = p['lam_mid']
    N = O.shape[0]
    
    numerator = O @ (O - lam_max * np.eye(N))
    denominator = lam_mid * (lam_mid - lam_max)
    
    return numerator / denominator


def verify_overlap_spectrum(O, n):
    """
    Check that the overlap matrix O has exactly the Johnson scheme eigenvalues.
    Returns dict with results.
    """
    p = johnson_parameters(n)
    evals = np.sort(np.linalg.eigvalsh(O))[::-1]
    
    # Expected: lam_max (mult 1), lam_mid (mult d_phys), 0 (rest)
    tol = 1e-8
    
    # Count eigenvalues near each expected value
    n_max = np.sum(np.abs(evals - p['lam_max']) < tol)
    n_mid = np.sum(np.abs(evals - p['lam_mid']) < tol)
    n_zero = np.sum(np.abs(evals) < tol)
    
    return {
        'lam_max_expected': p['lam_max'],
        'lam_max_found': evals[0],
        'lam_max_mult_expected': 1,
        'lam_max_mult_found': n_max,
        'lam_mid_expected': p['lam_mid'],
        'lam_mid_found': evals[1],
        'lam_mid_mult_expected': p['d_phys'],
        'lam_mid_mult_found': n_mid,
        'n_zero_found': n_zero,
        'n_zero_expected': p['d_kernel_total'],
        'all_correct': (n_max == 1 and n_mid == p['d_phys'] 
                        and n_zero == p['d_kernel_total']),
    }
