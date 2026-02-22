#!/usr/bin/env python3
"""
Dark Sector Structural Probe
=============================
Three computations:
  1. Cross-level coupling g²(internal/aperture/cross) for K₆→K₈ through K₁₄→K₁₆
  2. Mass scale at each level from edge Gram matrix spectral invariants
  3. Multi-level threshold model for d_eff at intermediate scales

Key technique: Edge Gram matrix G = E·Eᵀ (n_edges × n_edges) has same
nonzero eigenvalues as overlap matrix O = Eᵀ·E (n_matchings × n_matchings).
Accumulated incrementally — never stores all matchings.

Cross-level coupling via Frobenius inner product:
  ||O_cross||²_F = ⟨G_A, G_B⟩_F
  Decomposed by edge type (internal vs aperture).
"""

import numpy as np
import time
import sys

# ============================================================
# MATCHING GENERATION
# ============================================================

def generate_matchings(n_vertices):
    """Generate all perfect matchings of K_{n_vertices} as lists of sorted edge tuples."""
    def _gen(remaining):
        if len(remaining) == 0:
            yield ()
            return
        if len(remaining) == 2:
            yield ((remaining[0], remaining[1]),)
            return
        first = remaining[0]
        rest = remaining[1:]
        for i in range(len(rest)):
            partner = rest[i]
            new_remaining = rest[:i] + rest[i+1:]
            for sub in _gen(new_remaining):
                yield ((first, partner),) + sub
    
    vertices = list(range(n_vertices))
    yield from _gen(vertices)


def double_factorial(n):
    """(2n-1)!! = number of perfect matchings of K_{2n}."""
    result = 1
    for k in range(1, 2*n, 2):
        result *= k
    return result


# ============================================================
# EDGE GRAM MATRIX (INCREMENTAL)
# ============================================================

def edge_index_map(n_vertices):
    """Create (i,j) → index mapping for all edges of K_{n_vertices}."""
    idx = {}
    count = 0
    for i in range(n_vertices):
        for j in range(i+1, n_vertices):
            idx[(i,j)] = count
            count += 1
    return idx, count


def compute_gram_matrix(n_vertices, edge_idx, n_edges, extra_edge=None, 
                        progress_interval=50000, label=""):
    """
    Compute edge Gram matrix G[e,f] = #{matchings containing both e and f}.
    If extra_edge given, append it to each matching before recording.
    Uses incremental accumulation — O(1) memory in matchings.
    
    Returns: G (n_edges × n_edges), count
    """
    G = np.zeros((n_edges, n_edges), dtype=np.float64)
    count = 0
    t0 = time.time()
    
    # For efficiency, convert edge tuples to index arrays
    for matching in generate_matchings(n_vertices):
        edges = list(matching)
        if extra_edge is not None:
            edges.append(extra_edge)
        
        # Get indices
        eidxs = []
        for a, b in edges:
            key = (min(a,b), max(a,b))
            if key in edge_idx:
                eidxs.append(edge_idx[key])
        
        # Rank-1 update to G
        for i in eidxs:
            for j in eidxs:
                G[i, j] += 1.0
        
        count += 1
        if progress_interval and count % progress_interval == 0:
            elapsed = time.time() - t0
            print(f"    {label} {count:>10,d} matchings  ({elapsed:.1f}s)")
    
    elapsed = time.time() - t0
    if count > 10000:
        print(f"    {label} {count:>10,d} matchings DONE ({elapsed:.1f}s)")
    
    return G, count


# ============================================================
# PART 1: CROSS-LEVEL COUPLING DECOMPOSITION
# ============================================================

def cross_level_coupling(n_lower, verbose=True):
    """
    Compute cross-level coupling g² decomposition for K_{2n} → K_{2n+2}.
    
    n_lower: half the number of vertices of the lower level (e.g., 3 for K₆)
    
    Method:
    - Embed K_{2n} matchings into K_{2n+2} by adding edge (2n, 2n+1)
    - Compute Gram matrices G_A (embedded) and G_B (full K_{2n+2})
    - Frobenius inner product ⟨G_A, G_B⟩ decomposed by edge type
    """
    nv_low = 2 * n_lower       # vertices of lower level
    nv_high = 2 * (n_lower + 1) # vertices of higher level
    
    n_match_low = double_factorial(n_lower)
    n_match_high = double_factorial(n_lower + 1)
    
    if verbose:
        print(f"\n  K_{nv_low} → K_{nv_high}:  {n_match_low:,d} embedded × {n_match_high:,d} full")
    
    # Edge indexing in K_{nv_high} space
    edge_idx, n_edges = edge_index_map(nv_high)
    
    # New vertices
    new_v = (nv_low, nv_low + 1)
    extra_edge = (new_v[0], new_v[1])
    
    # Classify edges
    internal_idx = []
    aperture_idx = []
    for (a, b), idx in edge_idx.items():
        if a >= nv_low or b >= nv_low:
            aperture_idx.append(idx)
        else:
            internal_idx.append(idx)
    
    internal_idx = np.array(internal_idx)
    aperture_idx = np.array(aperture_idx)
    
    if verbose:
        print(f"    Edges: {n_edges} total, {len(internal_idx)} internal, {len(aperture_idx)} aperture")
    
    # G_A: Gram of embedded K_{2n} matchings (generated in K_{2n} then add extra edge)
    if verbose:
        print(f"    Computing G_A (embedded K_{nv_low})...")
    G_A, cnt_A = compute_gram_matrix(
        nv_low, edge_idx, n_edges, extra_edge=extra_edge,
        progress_interval=50000, label=f"G_A(K{nv_low}→K{nv_high})"
    )
    
    # G_B: Gram of all K_{nv_high} matchings
    if verbose:
        print(f"    Computing G_B (full K_{nv_high})...")
    G_B, cnt_B = compute_gram_matrix(
        nv_high, edge_idx, n_edges, extra_edge=None,
        progress_interval=50000, label=f"G_B(K{nv_high})"
    )
    
    # Frobenius inner product and decomposition
    total = np.sum(G_A * G_B)
    
    int_int = np.sum(G_A[np.ix_(internal_idx, internal_idx)] * 
                     G_B[np.ix_(internal_idx, internal_idx)])
    
    apt_apt = np.sum(G_A[np.ix_(aperture_idx, aperture_idx)] * 
                     G_B[np.ix_(aperture_idx, aperture_idx)])
    
    cross = total - int_int - apt_apt
    
    result = {
        'transition': f'K{nv_low}→K{nv_high}',
        'n_low': n_lower, 'n_high': n_lower + 1,
        'nv_low': nv_low, 'nv_high': nv_high,
        'n_match_low': cnt_A, 'n_match_high': cnt_B,
        'total': total,
        'internal': int_int,
        'aperture': apt_apt,
        'cross_term': cross,
        'g2_total': 1.0,  # by construction
        'g2_internal': int_int / total,
        'g2_aperture': apt_apt / total,
        'g2_cross': cross / total,
        'G_A': G_A, 'G_B': G_B,
    }
    
    if verbose:
        print(f"\n    COUPLING DECOMPOSITION:")
        print(f"      g²_total    = {result['g2_total']:.6f}")
        print(f"      g²_internal = {result['g2_internal']:.6f}  ({result['g2_internal']*100:.1f}%)")
        print(f"      g²_aperture = {result['g2_aperture']:.6f}  ({result['g2_aperture']*100:.1f}%)")
        print(f"      g²_cross    = {result['g2_cross']:.6f}  ({result['g2_cross']*100:.1f}%)")
    
    return result


# ============================================================
# PART 2: MASS SCALES FROM SPECTRAL INVARIANTS
# ============================================================

def spectral_invariants(n, verbose=True):
    """
    Compute spectral invariants of K_{2n} matching algebra from edge Gram matrix.
    
    a₂ ~ Tr(G²)/nm² = sum of squared eigenvalues of G, normalized
    a₄ ~ Tr(G⁴)/nm⁴ = sum of fourth-power eigenvalues, normalized  
    R = a₄ / a₂²
    """
    nv = 2 * n
    n_match = double_factorial(n)
    edge_idx, n_edges = edge_index_map(nv)
    
    if verbose:
        print(f"\n  K_{nv}: {n_match:,d} matchings, {n_edges} edges")
    
    G, cnt = compute_gram_matrix(
        nv, edge_idx, n_edges, extra_edge=None,
        progress_interval=50000, label=f"G(K{nv})"
    )
    
    # Eigenvalues of G
    evals = np.linalg.eigvalsh(G)
    evals = np.sort(evals)[::-1]  # descending
    
    # Spectral invariants (normalized by matching count)
    nm = cnt
    
    # Tr(G) = sum of all entries of G = sum over all matchings of (n edges)² 
    # Actually Tr(G) = sum_e #{matchings containing e}
    trG = np.trace(G)
    
    # a₂ analog: Tr(G) / n_matchings 
    # This is the average edge multiplicity
    a2_raw = trG / nm
    
    # Frobenius norm squared: Tr(G²) = sum of squared eigenvalues
    trG2 = np.sum(evals**2)
    a2_spectral = trG2 / nm**2
    
    # Fourth moment
    trG4 = np.sum(evals**4)
    a4_spectral = trG4 / nm**4
    
    # R ratio
    R = a4_spectral / a2_spectral**2 if a2_spectral > 0 else 0
    
    # Maximum eigenvalue (dominant mode)
    lambda_max = evals[0]
    
    # Effective rank
    eff_rank = (np.sum(evals[evals > 1e-10]))**2 / np.sum(evals[evals > 1e-10]**2)
    
    result = {
        'level': f'K_{nv}',
        'n': n, 'nv': nv,
        'n_matchings': nm,
        'n_edges': n_edges,
        'eigenvalues': evals,
        'lambda_max': lambda_max,
        'a2_raw': a2_raw,
        'a2_spectral': a2_spectral,
        'a4_spectral': a4_spectral,
        'R': R,
        'eff_rank': eff_rank,
        'G': G,
    }
    
    if verbose:
        print(f"    λ_max = {lambda_max:.4f}")
        print(f"    a₂ (raw) = {a2_raw:.4f}")
        print(f"    a₂ (spectral) = {a2_spectral:.6f}")
        print(f"    a₄ (spectral) = {a4_spectral:.8f}")
        print(f"    R = a₄/a₂² = {R:.8f}")
        print(f"    eff_rank = {eff_rank:.2f}")
        top5 = evals[:min(5, len(evals))]
        print(f"    Top eigenvalues: {', '.join(f'{v:.2f}' for v in top5)}")
    
    return result


# ============================================================
# PART 3: MULTI-LEVEL d_eff MODEL  
# ============================================================

def deff_single_level(eps, Var_E1=0.17355, H0=0.5*(1 + np.log(np.pi/2))):
    """
    f(ε) = A / (ε²(ln ε + H₀))  for single K₄ channel with arrow of time.
    A = Var_k[E₁] / 2 (analytical) or 0.0867 (measured).
    """
    A = 0.0867  # measured value
    if eps <= 1:
        return 0.0
    ln_eps = np.log(eps)
    return A / (eps**2 * (ln_eps + H0))


def multi_level_deff(eps, mass_scales, coupling_fractions, 
                     m_H=126.1, Lambda_cutoff=1e19):
    """
    Model: d_eff(ε) = 3 - f_K4(ε) - Σ_n δf_n(ε)
    
    Each level n with mass m_n contributes a threshold correction:
    δf_n(ε) = g²_aperture(n) × f_K4(ε) × (m_n/Λ)^{-4} × S(μ/m_n)
    
    where S is a step function (level active above threshold)
    and the suppression (m_n/Λ)^{-4} comes from spectral action.
    
    The question: does Σ δf_n vary enough between galactic and cosmological ε
    to mimic dark matter effects?
    """
    f_base = deff_single_level(eps)
    
    corrections = []
    for i, (mass, g2_apt) in enumerate(zip(mass_scales, coupling_fractions)):
        # The correction from this frozen level
        # Suppression: (m/Λ)^4 from spectral action freeze-out  
        # Enhancement: g²_aperture from cross-level coupling
        suppression = (mass / Lambda_cutoff)**4
        delta_f = g2_apt * f_base * suppression
        corrections.append({
            'level': i,
            'mass': mass,
            'g2_aperture': g2_apt,
            'suppression': suppression,
            'delta_f': delta_f,
        })
    
    total_correction = sum(c['delta_f'] for c in corrections)
    
    return {
        'eps': eps,
        'f_base': f_base,
        'f_total': f_base + total_correction,
        'corrections': corrections,
        'total_correction': total_correction,
        'fractional_correction': total_correction / f_base if f_base > 0 else 0,
    }


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 72)
    print("  DARK SECTOR STRUCTURAL PROBE")
    print("  Cross-level coupling | Mass scales | Multi-level d_eff")
    print("=" * 72)
    
    max_n = 7  # Default: through K₁₄
    if "--k16" in sys.argv:
        max_n = 8
    if "--quick" in sys.argv:
        max_n = 5  # Through K₁₀ only
    
    # ================================================================
    # PART 1: CROSS-LEVEL COUPLING
    # ================================================================
    
    print("\n" + "=" * 72)
    print("  PART 1: CROSS-LEVEL COUPLING DECOMPOSITION")
    print("=" * 72)
    
    coupling_results = []
    for n in range(3, max_n):  # K₆→K₈ through K_{2(max_n-1)}→K_{2·max_n}
        t0 = time.time()
        result = cross_level_coupling(n, verbose=True)
        dt = time.time() - t0
        result['time'] = dt
        coupling_results.append(result)
    
    # Summary table
    print(f"\n{'='*72}")
    print("  COUPLING SUMMARY")
    print(f"{'='*72}")
    print(f"\n  {'Transition':<12s} {'g²_int':>8s} {'g²_apt':>8s} {'g²_cross':>8s} {'int_frac':>9s} {'time':>8s}")
    print("  " + "-" * 58)
    for r in coupling_results:
        print(f"  {r['transition']:<12s} {r['g2_internal']:8.4f} {r['g2_aperture']:8.4f} "
              f"{r['g2_cross']:8.4f} {r['g2_internal']*100:8.1f}% {r['time']:7.1f}s")
    
    # Asymptotic trend
    if len(coupling_results) >= 3:
        ints = [r['g2_internal'] for r in coupling_results]
        apts = [r['g2_aperture'] for r in coupling_results]
        print(f"\n  TREND ANALYSIS:")
        print(f"    g²_internal sequence: {', '.join(f'{v:.4f}' for v in ints)}")
        print(f"    g²_aperture sequence: {', '.join(f'{v:.4f}' for v in apts)}")
        
        # Fit exponential approach to limit: g²_int = L - C·r^n
        if len(ints) >= 4:
            # Simple: check ratio of increments
            deltas = [ints[i+1] - ints[i] for i in range(len(ints)-1)]
            if len(deltas) >= 2 and all(d > 0 for d in deltas):
                ratios = [deltas[i+1]/deltas[i] for i in range(len(deltas)-1)]
                avg_ratio = np.mean(ratios)
                # Extrapolated limit: L = last + delta_last / (1 - ratio)
                if avg_ratio < 1:
                    L = ints[-1] + deltas[-1] * avg_ratio / (1 - avg_ratio)
                    print(f"    Increment ratios: {', '.join(f'{r:.3f}' for r in ratios)}")
                    print(f"    Extrapolated g²_int limit: {L:.4f}")
                    print(f"    Extrapolated g²_apt limit: {1 - L:.4f}" if L < 1 else "")
    
    # ================================================================
    # PART 2: MASS SCALES
    # ================================================================
    
    print(f"\n{'='*72}")
    print("  PART 2: SPECTRAL INVARIANTS & MASS SCALES")
    print(f"{'='*72}")
    
    spectral_results = []
    for n in range(3, max_n + 1):
        t0 = time.time()
        result = spectral_invariants(n, verbose=True)
        result['time'] = time.time() - t0
        spectral_results.append(result)
    
    # Mass calibration: K₆ a₂_spectral → 126.1 GeV
    a2_K6 = spectral_results[0]['a2_spectral']
    lmax_K6 = spectral_results[0]['lambda_max']
    
    print(f"\n{'='*72}")
    print("  MASS SCALE CALIBRATION (K₆ = 126.1 GeV)")
    print(f"{'='*72}")
    print(f"\n  {'Level':<6s} {'n_match':>10s} {'λ_max':>10s} {'a₂':>12s} {'R':>10s} {'m (GeV)':>10s} {'m/m_H':>8s}")
    print("  " + "-" * 70)
    
    mass_scales = []
    for r in spectral_results:
        # Mass ~ sqrt(a₂) calibrated to K₆
        ratio = r['a2_spectral'] / a2_K6
        mass = 126.1 * np.sqrt(ratio)
        mass_scales.append(mass)
        print(f"  {r['level']:<6s} {r['n_matchings']:>10,d} {r['lambda_max']:>10.2f} "
              f"{r['a2_spectral']:>12.6f} {r['R']:>10.6f} {mass:>10.1f} {mass/126.1:>8.2f}x")
    
    # Also calibrate by λ_max
    print(f"\n  Alternative calibration (λ_max):")
    print(f"  {'Level':<6s} {'λ_max':>10s} {'m_alt (GeV)':>12s}")
    print("  " + "-" * 32)
    for r in spectral_results:
        m_alt = 126.1 * np.sqrt(r['lambda_max'] / lmax_K6)
        print(f"  {r['level']:<6s} {r['lambda_max']:>10.2f} {m_alt:>12.1f}")
    
    # ================================================================
    # PART 3: MULTI-LEVEL d_eff MODEL
    # ================================================================
    
    print(f"\n{'='*72}")
    print("  PART 3: MULTI-LEVEL d_eff THRESHOLD MODEL")
    print(f"{'='*72}")
    
    # Use aperture fractions from coupling results
    g2_apts = [r['g2_aperture'] for r in coupling_results]
    # Mass scales for the HIGHER level in each transition
    transition_masses = mass_scales[1:]  # K₈, K₁₀, K₁₂, ...
    
    # Truncate to matching lengths
    n_levels = min(len(g2_apts), len(transition_masses))
    g2_apts = g2_apts[:n_levels]
    transition_masses = transition_masses[:n_levels]
    
    print(f"\n  Level contributions to d_eff:")
    print(f"  {'Level':>6s} {'mass (GeV)':>12s} {'g²_apt':>8s} {'(m/Λ)⁴':>14s} {'δf/f':>14s}")
    print("  " + "-" * 56)
    
    Lambda = 1e19  # Planck scale in GeV
    eps_values = {
        'galactic (1 kpc)':   1e57,
        'cluster (1 Mpc)':    1e60,
        'Hubble (14 Gly)':    1e61,
    }
    
    for i in range(n_levels):
        supp = (transition_masses[i] / Lambda)**4
        delta_rel = g2_apts[i] * supp
        print(f"  K_{2*(4+i):<4d} {transition_masses[i]:>12.1f} {g2_apts[i]:>8.4f} "
              f"{supp:>14.2e} {delta_rel:>14.2e}")
    
    print(f"\n  d_eff at different scales:")
    print(f"  {'Scale':<22s} {'ε':>10s} {'f_base':>14s} {'f_total':>14s} {'δf/f':>14s}")
    print("  " + "-" * 76)
    
    for name, eps in eps_values.items():
        result = multi_level_deff(eps, transition_masses, g2_apts, Lambda_cutoff=Lambda)
        print(f"  {name:<22s} {eps:>10.0e} {result['f_base']:>14.6e} "
              f"{result['f_total']:>14.6e} {result['fractional_correction']:>14.6e}")
    
    # ================================================================
    # VERDICT
    # ================================================================
    
    print(f"\n{'='*72}")
    print("  VERDICT")
    print(f"{'='*72}")
    
    # Check if tower corrections are relevant
    result_gal = multi_level_deff(1e57, transition_masses, g2_apts, Lambda_cutoff=Lambda)
    result_hub = multi_level_deff(1e61, transition_masses, g2_apts, Lambda_cutoff=Lambda)
    
    # Dark matter needs ~10% variation in effective gravity between galaxy core and outskirts
    # That requires δf to vary by O(f) between these scales
    frac_gal = result_gal['fractional_correction']
    frac_hub = result_hub['fractional_correction']
    
    print(f"""
  The multi-level tower corrections to d_eff are:
    At galactic scale: δf/f = {frac_gal:.2e}
    At Hubble scale:   δf/f = {frac_hub:.2e}
    
  Dark matter requires ~10% variation in effective gravity.
  The tower provides {abs(frac_gal):.2e} fractional correction.
""")
    
    if abs(frac_gal) < 1e-50:
        print("  KILL: Tower corrections are negligible (< 10⁻⁵⁰).")
        print("  The (m/Λ)⁴ suppression kills all particle-like contributions.")
        print("  If geometric dark matter exists in this framework, it CANNOT come")
        print("  from the K₂ₙ tower thresholds.")
        print()
        print("  SURVIVING PATH: Local d_eff variation from boundary entanglement")
        print("  geometry — the d_eff correction must be environment-dependent,")
        print("  not level-dependent.")
    else:
        print("  INTERESTING: Tower corrections are non-negligible.")
        print("  Further investigation warranted.")
    
    print(f"\n{'='*72}")
    
    return coupling_results, spectral_results


if __name__ == "__main__":
    main()
