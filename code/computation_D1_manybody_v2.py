#!/usr/bin/env python3
"""
COMPUTATION D1-MB v2: Many-Body K₄ Model — CORRECTED
======================================================

CRITICAL FIX (Feb 2026): build_triangular_lattice() now deduplicates
undirected bonds. The v1 code listed every bond in both directions with
the same ℤ₃ phase, causing complex phases to cancel (Im/Re = 0.000000
on the 2×2 torus). This was a code artifact, not physics.

See D1MB_bug_analysis.md for full diagnosis.

Bond counts after fix:
  3 sites (triangle):  3 bonds  (unchanged — was hand-written)
  4 sites (2×2 torus): 6 bonds  (was 12 — all involutive)
  6 sites (3×2 torus): 9 bonds  (was 18 — direction 1 involutive)

Purpose: Determine whether the interacting K₄ boundary theory has
anomalous thermalization properties that could support cosmological cycling.

What this computes:
  1. Many-body energy spectrum E_n(U) for the K₄ Hubbard model
  2. Level statistics r-ratio (GUE vs Poisson → thermalization vs localization)
  3. Spectral asymmetry in the many-body spectrum (T-breaking signature)
  4. Entanglement entropy of half-system at various energies
  5. Comparison WITH vs WITHOUT ℤ₃ flux

Usage:
  python computation_D1_manybody_v2.py --nsites 4 --U 0 1 2 3 4 5 6 7 8 --full --noflux
  python computation_D1_manybody_v2.py --nsites 6 --U 0 2 4 6 8 --nev 500

Flags:
  --nsites N    Number of lattice sites (3, 4, or 6)
  --U val ...   List of U/t values to scan
  --nev N       Number of eigenvalues for sparse solver (default: 200)
  --full        Use full (dense) diagonalization (only for nsites ≤ 4)
  --noflux      Also run without ℤ₃ flux for comparison
  --outdir DIR  Output directory (default: ./d1mb_results_v2/)
  --no-ee       Skip entanglement entropy (faster)
  --skip-selftest  Skip the single-particle self-test

KEY CAVEATS for interpreting results:
  - At U=0 the system is integrable → Poisson expected regardless of topology
  - Small systems have symmetry sectors that contaminate level statistics.
    The 3-site results have large error bars. 4-site is the minimum for
    meaningful level statistics. 6-site is where things get reliable.
  - Must eventually resolve V₄ × translation symmetry sectors for clean ETH test.
    This code does NOT do that yet — it's the full Hilbert space.
  - With the corrected code, ⟨r⟩ should be compared to GUE (0.5996), not
    GOE (0.5307), since T-symmetry is genuinely broken by ℤ₃ flux.

Requirements: numpy, scipy
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.special import comb
import argparse
import os
import time
from itertools import combinations

# ============================================================
# MODEL DEFINITION
# ============================================================

omega = np.exp(2j * np.pi / 3)

# K₄ matching matrices (4×4)
M1 = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]], dtype=complex)
M2 = np.array([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]], dtype=complex)
M3 = np.array([[0,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,0]], dtype=complex)

a1 = np.array([1.0, 0.0])
a2 = np.array([-0.5, np.sqrt(3)/2])


def build_triangular_lattice(nsites):
    """
    Build triangular lattice with PBC. Returns positions, neighbor_list.
    
    CORRECTED: Tracks seen undirected bonds as (min(i,j), max(i,j)) and
    skips duplicates. This prevents the involution bug on small tori where
    stepping forward from i reaches j AND stepping forward from j reaches i
    via the same direction, causing double-listing that kills complex phases.
    """
    if nsites == 3:
        positions = [np.array([0,0]), a1, a2]
        neighbors = [
            (0, 1, 0),  # δ₁ direction → M₁, phase 1
            (0, 2, 1),  # δ₂ direction → M₂, phase ω
            (1, 2, 2),  # δ₃ direction → M₃, phase ω²
        ]
        return positions, neighbors
    
    elif nsites == 4:
        Lx, Ly = 2, 2
    elif nsites == 6:
        Lx, Ly = 3, 2
    else:
        raise ValueError(f"nsites={nsites} not supported. Use 3, 4, or 6.")
    
    positions = []
    site_index = {}
    idx = 0
    for iy in range(Ly):
        for ix in range(Lx):
            positions.append(ix * a1 + iy * a2)
            site_index[(ix, iy)] = idx
            idx += 1
    
    # CRITICAL FIX: deduplicate undirected bonds
    neighbors = []
    seen_bonds = set()
    for iy in range(Ly):
        for ix in range(Lx):
            i = site_index[(ix, iy)]
            for j, d in [(site_index[((ix+1)%Lx, iy)], 0),
                         (site_index[(ix, (iy+1)%Ly)], 1),
                         (site_index[((ix-1)%Lx, (iy-1)%Ly)], 2)]:
                bond = (min(i, j), max(i, j))
                if bond not in seen_bonds:
                    seen_bonds.add(bond)
                    neighbors.append((i, j, d))
    
    return positions, neighbors


# ============================================================
# SELF-TEST: Verify the fix before any computation
# ============================================================

def verify_fix():
    """
    Build the single-particle hopping matrix for 4 sites and verify:
    1. Exactly 6 bonds (not 12)
    2. Im/Re > 0 (complex Hermitian, not real)
    3. Hermitian
    4. 4 distinct eigenvalues (not 3)
    
    Hard-fails if any check fails — we do NOT proceed with corrupted physics.
    """
    print("=" * 70)
    print("  SELF-TEST: Verifying bond-deduplication fix")
    print("=" * 70)
    
    nsites = 4
    n_orb = 4
    n_sp = nsites * n_orb  # 16
    
    positions, neighbors = build_triangular_lattice(nsites)
    n_bonds = len(neighbors)
    print(f"  Bond count: {n_bonds} (expected: 6)")
    assert n_bonds == 6, f"FAIL: Expected 6 bonds, got {n_bonds}. Fix not applied!"
    
    # Build single-particle hopping matrix
    M_list = [M1, M2, M3]
    phases = [1.0, omega, omega**2]
    
    H_sp = np.zeros((n_sp, n_sp), dtype=complex)
    for (site_i, site_j, direction) in neighbors:
        M = M_list[direction]
        phase = phases[direction]
        for alpha in range(n_orb):
            for beta in range(n_orb):
                t_hop = phase * M[alpha, beta]
                if abs(t_hop) < 1e-12:
                    continue
                sp_i = site_i * n_orb + alpha
                sp_j = site_j * n_orb + beta
                H_sp[sp_i, sp_j] += t_hop
                H_sp[sp_j, sp_i] += np.conj(t_hop)
    
    # Check Im/Re
    im_norm = np.linalg.norm(np.imag(H_sp))
    re_norm = np.linalg.norm(np.real(H_sp))
    im_re = im_norm / max(re_norm, 1e-15)
    print(f"  ||Im(H_sp)||/||Re(H_sp)|| = {im_re:.6f} (expected: ~1.0)")
    assert im_re > 0.5, f"FAIL: Im/Re = {im_re:.6f}, H is still (nearly) real!"
    
    # Check Hermitian
    herm_err = np.linalg.norm(H_sp - H_sp.conj().T)
    print(f"  ||H - H†|| = {herm_err:.2e} (expected: ~0)")
    assert herm_err < 1e-10, f"FAIL: H not Hermitian, error = {herm_err:.2e}"
    
    # Check eigenvalue count
    evals = np.linalg.eigvalsh(H_sp)
    unique_evals = len(set(np.round(evals, 6)))
    print(f"  Distinct eigenvalues: {unique_evals} (expected: ≥4)")
    assert unique_evals >= 4, f"FAIL: Only {unique_evals} distinct eigenvalues"
    
    print("  ✓ All self-tests PASSED. Fix verified.")
    print("=" * 70)
    print()


# ============================================================
# MANY-BODY HILBERT SPACE
# ============================================================

def build_basis(n_sp, n_particles):
    """Build Fock-space basis as sorted list of occupation-number integers."""
    basis = []
    for combo in combinations(range(n_sp), n_particles):
        state = 0
        for orb in combo:
            state |= (1 << orb)
        basis.append(state)
    basis.sort()
    return np.array(basis, dtype=np.int64)


def state_to_index(basis):
    """Fock state → basis index."""
    return {int(state): idx for idx, state in enumerate(basis)}


def sp_index(site, orbital, n_orb=4):
    """(site, orbital) → single-particle index."""
    return site * n_orb + orbital


# ============================================================
# HAMILTONIAN CONSTRUCTION
# ============================================================

def build_hamiltonian(nsites, U, use_flux=True, basis=None, state_map=None):
    """
    Build the many-body K₄ Hubbard Hamiltonian as a sparse matrix.
    
    H = H₀ + U·H_int
    H₀ = Σ_bonds t_{αβ} c†_{iα} c_{jβ} + h.c.
    H_int = U Σ_i Σ_{α<β} n_{iα} n_{iβ}
    """
    n_orb = 4
    n_sp = nsites * n_orb
    n_particles = nsites * 2  # half-filling
    
    if basis is None:
        print(f"  Building basis: {nsites} sites, {n_sp} orbitals, "
              f"{n_particles} particles...")
        basis = build_basis(n_sp, n_particles)
        state_map = state_to_index(basis)
    
    dim = len(basis)
    print(f"  Hilbert space dimension: {dim:,}")
    
    M_list = [M1, M2, M3]
    phases = [1.0, omega, omega**2] if use_flux else [1.0, 1.0, 1.0]
    _, neighbors = build_triangular_lattice(nsites)
    
    rows, cols, vals = [], [], []
    
    # --- Hopping term ---
    for (site_i, site_j, direction) in neighbors:
        M = M_list[direction]
        phase = phases[direction]
        
        for alpha in range(n_orb):
            for beta in range(n_orb):
                t_hop = phase * M[alpha, beta]
                if abs(t_hop) < 1e-12:
                    continue
                
                sp_i = sp_index(site_i, alpha)
                sp_j = sp_index(site_j, beta)
                
                for idx in range(dim):
                    state = basis[idx]
                    
                    # c†_{iα} c_{jβ}: need sp_j occupied, sp_i empty
                    if (state >> sp_j) & 1 and not (state >> sp_i) & 1:
                        new_state = state ^ (1 << sp_j) ^ (1 << sp_i)
                        
                        # Fermionic sign
                        lo, hi = min(sp_i, sp_j), max(sp_i, sp_j)
                        mask = 0
                        for b in range(lo+1, hi):
                            mask |= (1 << b)
                        sign = (-1) ** bin(state & mask).count('1')
                        
                        new_idx = state_map.get(int(new_state))
                        if new_idx is not None:
                            rows.append(new_idx)
                            cols.append(idx)
                            vals.append(t_hop * sign)
                    
                    # h.c.: c†_{jβ} c_{iα}: need sp_i occupied, sp_j empty
                    if (state >> sp_i) & 1 and not (state >> sp_j) & 1:
                        new_state = state ^ (1 << sp_i) ^ (1 << sp_j)
                        
                        lo, hi = min(sp_i, sp_j), max(sp_i, sp_j)
                        mask = 0
                        for b in range(lo+1, hi):
                            mask |= (1 << b)
                        sign = (-1) ** bin(state & mask).count('1')
                        
                        new_idx = state_map.get(int(new_state))
                        if new_idx is not None:
                            rows.append(new_idx)
                            cols.append(idx)
                            vals.append(np.conj(t_hop) * sign)
    
    # --- Interaction term (diagonal) ---
    if abs(U) > 1e-12:
        for idx in range(dim):
            state = basis[idx]
            diag_val = 0.0
            for site in range(nsites):
                n_occ = 0
                for alpha in range(n_orb):
                    if (state >> sp_index(site, alpha)) & 1:
                        n_occ += 1
                diag_val += U * n_occ * (n_occ - 1) / 2.0
            rows.append(idx)
            cols.append(idx)
            vals.append(diag_val)
    
    H = sparse.csr_matrix((vals, (rows, cols)), shape=(dim, dim), dtype=complex)
    H = (H + H.conj().T) / 2.0  # enforce hermiticity
    
    return H, basis, state_map


# ============================================================
# ANALYSIS FUNCTIONS
# ============================================================

def level_statistics(energies, n_skip=0):
    """
    Mean r-ratio for level spacing statistics.
    GUE (T-broken, thermalizing):   ⟨r⟩ ≈ 0.5996
    GOE (T-symmetric, thermalizing): ⟨r⟩ ≈ 0.5307  
    Poisson (non-thermalizing):      ⟨r⟩ ≈ 0.3863
    """
    E = np.sort(np.real(energies))
    if n_skip > 0:
        E = E[n_skip:-n_skip] if len(E) > 2*n_skip else E
    
    spacings = np.diff(E)
    spacings = spacings[spacings > 1e-12]
    
    if len(spacings) < 3:
        return 0, 0, np.array([])
    
    r_values = np.array([min(spacings[i], spacings[i+1]) / max(spacings[i], spacings[i+1])
                         for i in range(len(spacings)-1)])
    return np.mean(r_values), np.std(r_values)/np.sqrt(len(r_values)), r_values


def entanglement_entropy(psi, basis, nsites, n_orb=4):
    """
    Von Neumann entanglement entropy for bipartition into 
    first nsites//2 sites vs rest.
    """
    n_A = (nsites // 2) * n_orb
    
    from collections import defaultdict
    groups = defaultdict(list)
    for idx, state in enumerate(basis):
        A_config = int(state) & ((1 << n_A) - 1)
        B_config = int(state) >> n_A
        c = psi[idx]
        if abs(c) > 1e-15:
            groups[A_config].append((B_config, c))
    
    A_configs = sorted(groups.keys())
    A_map = {cfg: i for i, cfg in enumerate(A_configs)}
    dim_A = len(A_configs)
    
    rho_A = np.zeros((dim_A, dim_A), dtype=complex)
    for a1 in A_configs:
        i = A_map[a1]
        b_dict = {b: np.conj(c) for b, c in groups[a1]}
        for a2 in A_configs:
            j = A_map[a2]
            for b, c in groups[a2]:
                if b in b_dict:
                    rho_A[i, j] += b_dict[b] * c
    
    evals = np.real(np.linalg.eigvalsh(rho_A))
    evals = evals[evals > 1e-15]
    return -np.sum(evals * np.log(evals))


# ============================================================
# MAIN COMPUTATION
# ============================================================

def run_analysis(nsites, U_values, nev=200, full_diag=False, 
                 run_noflux=False, compute_ee=True, outdir='d1mb_results_v2'):
    
    os.makedirs(outdir, exist_ok=True)
    
    n_orb = 4
    n_sp = nsites * n_orb
    n_particles = nsites * 2
    dim = int(comb(n_sp, n_particles, exact=True))
    
    print("=" * 70)
    print(f"  K₄ MANY-BODY COMPUTATION D1-MB v2 (CORRECTED)")
    print(f"  Sites: {nsites}  |  Orbitals: {n_sp}  |  Particles: {n_particles}")
    print(f"  Hilbert space: {dim:,}")
    print(f"  U values: {U_values}")
    print(f"  Method: {'full diag' if full_diag else f'Lanczos (nev={nev})'}")
    print("=" * 70)
    
    if dim > 50000 and full_diag:
        print(f"  WARNING: dim={dim:,} too large for dense diag. Switching to Lanczos.")
        full_diag = False
    
    # Build basis once
    t0 = time.time()
    basis = build_basis(n_sp, n_particles)
    state_map = state_to_index(basis)
    print(f"  Basis built in {time.time()-t0:.1f}s")
    
    r_GUE, r_GOE, r_Poi = 0.5996, 0.5307, 0.3863
    results = []
    
    for U in U_values:
        for use_flux in ([True, False] if run_noflux else [True]):
            flux_label = "WITH ℤ₃ flux" if use_flux else "WITHOUT flux (control)"
            print(f"\n{'='*70}")
            print(f"  U/t = {U:.2f}  ({flux_label})")
            print(f"{'='*70}")
            
            t0 = time.time()
            H, _, _ = build_hamiltonian(nsites, U, use_flux=use_flux,
                                         basis=basis, state_map=state_map)
            build_time = time.time() - t0
            print(f"  H built in {build_time:.1f}s  (nnz={H.nnz:,})")
            
            # T-breaking measure
            H_d = H.toarray() if dim < 20000 else None
            if H_d is not None:
                im_re = np.linalg.norm(np.imag(H_d)) / max(np.linalg.norm(np.real(H_d)), 1e-15)
            else:
                # Sample random entries
                sample = H.toarray()[:1000, :1000] if dim > 1000 else H.toarray()
                im_re = np.linalg.norm(np.imag(sample)) / max(np.linalg.norm(np.real(sample)), 1e-15)
            print(f"  ||Im(H)||/||Re(H)|| = {im_re:.6f}  "
                  f"({'T-BROKEN' if im_re > 0.01 else 'T-symmetric'})")
            
            # RUNTIME ASSERTION: If flux case at U=0 has Im/Re ~ 0, the fix failed
            if use_flux and abs(U) < 1e-12 and im_re < 0.01:
                raise RuntimeError(
                    f"CRITICAL: Im/Re = {im_re:.6f} at U=0 with ℤ₃ flux. "
                    f"Bond deduplication fix not working! Aborting."
                )
            
            # Diagonalize
            t0 = time.time()
            if full_diag and dim < 50000:
                if H_d is None:
                    H_d = H.toarray()
                energies, eigvecs = np.linalg.eigh(H_d)
            else:
                eigvecs = None
                try:
                    energies, eigvecs = eigsh(H, k=min(nev, dim-2), 
                                              sigma=0, which='LM')
                    idx_s = np.argsort(energies)
                    energies = energies[idx_s]
                    eigvecs = eigvecs[:, idx_s]
                except Exception as e:
                    print(f"  Shift-invert failed: {e}")
                    print(f"  Trying extremal eigenvalues...")
                    energies, eigvecs = eigsh(H, k=min(nev, dim-2), which='SA')
                    idx_s = np.argsort(energies)
                    energies = energies[idx_s]
                    eigvecs = eigvecs[:, idx_s]
            
            n_ev = len(energies)
            print(f"  Diagonalized in {time.time()-t0:.1f}s  ({n_ev} eigenvalues)")
            print(f"  E_min={energies[0]:.4f}  E_max={energies[-1]:.4f}")
            
            # Level statistics
            n_skip = max(n_ev // 20, 2)
            mean_r, std_r, r_vals = level_statistics(energies, n_skip=n_skip)
            
            dists = {'GUE': abs(mean_r-r_GUE), 'GOE': abs(mean_r-r_GOE), 
                     'Poisson': abs(mean_r-r_Poi)}
            verdict = min(dists, key=dists.get)
            
            print(f"\n  LEVEL STATISTICS ({n_ev - 2*n_skip} levels):")
            print(f"    ⟨r⟩ = {mean_r:.4f} ± {std_r:.4f}")
            print(f"    GUE={r_GUE:.4f}  GOE={r_GOE:.4f}  Poisson={r_Poi:.4f}")
            print(f"    → {verdict}")
            
            # Entanglement entropy
            S_gs = S_mid = None
            if compute_ee and eigvecs is not None and dim < 20000:
                S_gs = entanglement_entropy(eigvecs[:, 0], basis, nsites)
                S_max = (nsites // 2) * n_orb * np.log(2)
                mid = n_ev // 2
                S_mid = entanglement_entropy(eigvecs[:, mid], basis, nsites)
                print(f"\n  ENTANGLEMENT ENTROPY:")
                print(f"    S_gs = {S_gs:.4f}  ({S_gs/S_max:.3f} of max)")
                print(f"    S_mid = {S_mid:.4f}  ({S_mid/S_max:.3f} of max) at E={energies[mid]:.4f}")
            
            results.append({
                'U': U, 'flux': use_flux, 'mean_r': mean_r, 'std_r': std_r,
                'im_re': im_re, 'verdict': verdict, 'S_gs': S_gs, 'S_mid': S_mid,
                'E_min': energies[0], 'E_max': energies[-1],
            })
            
            np.savez(os.path.join(outdir, 
                     f'spectrum_N{nsites}_U{U:.1f}_{"flux" if use_flux else "noflux"}.npz'),
                     energies=energies, U=U, nsites=nsites, mean_r=mean_r)
    
    # ============================================================
    # SUMMARY
    # ============================================================
    
    print(f"\n{'='*70}")
    print(f"  SUMMARY — D1-MB v2 (CORRECTED)")
    print(f"{'='*70}")
    print(f"  {'U/t':>5s}  {'Flux':>4s}  {'⟨r⟩':>7s}  {'±':>6s}  "
          f"{'Im/Re':>7s}  {'Verdict':>10s}  {'S_gs':>6s}  {'S_mid':>6s}")
    print(f"  {'-'*5}  {'-'*4}  {'-'*7}  {'-'*6}  {'-'*7}  {'-'*10}  {'-'*6}  {'-'*6}")
    
    for r in results:
        f_str = "ℤ₃" if r['flux'] else "—"
        s1 = f"{r['S_gs']:.3f}" if r['S_gs'] is not None else "—"
        s2 = f"{r['S_mid']:.3f}" if r['S_mid'] is not None else "—"
        print(f"  {r['U']:5.1f}  {f_str:>4s}  {r['mean_r']:7.4f}  {r['std_r']:6.4f}  "
              f"{r['im_re']:7.4f}  {r['verdict']:>10s}  {s1:>6s}  {s2:>6s}")
    
    # ============================================================
    # FLUX vs NO-FLUX COMPARISON TABLE
    # ============================================================
    
    if run_noflux:
        print(f"\n{'='*70}")
        print(f"  FLUX vs NO-FLUX COMPARISON")
        print(f"{'='*70}")
        print(f"  {'U/t':>5s}  {'⟨r⟩_flux':>9s}  {'⟨r⟩_noflux':>11s}  "
              f"{'Δ⟨r⟩':>7s}  {'Im/Re_flux':>11s}  {'Im/Re_noflux':>13s}")
        print(f"  {'-'*5}  {'-'*9}  {'-'*11}  {'-'*7}  {'-'*11}  {'-'*13}")
        
        flux_by_U = {r['U']: r for r in results if r['flux']}
        noflux_by_U = {r['U']: r for r in results if not r['flux']}
        
        for U in sorted(flux_by_U.keys()):
            if U in noflux_by_U:
                rf = flux_by_U[U]
                rn = noflux_by_U[U]
                delta_r = rf['mean_r'] - rn['mean_r']
                print(f"  {U:5.1f}  {rf['mean_r']:9.4f}  {rn['mean_r']:11.4f}  "
                      f"{delta_r:+7.4f}  {rf['im_re']:11.6f}  {rn['im_re']:13.6f}")
    
    # ============================================================
    # DECISION-POINT ANALYSIS
    # ============================================================
    
    flux_results = [r for r in results if r['flux']]
    if len(flux_results) > 1:
        print(f"\n{'='*70}")
        print(f"  DECISION-POINT ANALYSIS")
        print(f"{'='*70}")
        
        r_vals = [r['mean_r'] for r in flux_results]
        im_vals = [r['im_re'] for r in flux_results]
        U_vals = [r['U'] for r in flux_results]
        
        print(f"  Δ⟨r⟩ across U scan: {max(r_vals)-min(r_vals):.4f}")
        
        # Does Im/Re decrease with U?
        if im_vals[0] > im_vals[-1]:
            print(f"  Im/Re ratio DECREASES with U: {im_vals[0]:.4f} → {im_vals[-1]:.4f}")
            print(f"  → Interaction DILUTES the ℤ₃ T-breaking effect")
        else:
            print(f"  Im/Re ratio INCREASES with U: {im_vals[0]:.4f} → {im_vals[-1]:.4f}")
            print(f"  → Interaction ENHANCES or preserves T-breaking")
        
        # Is there a peak in ⟨r⟩ suggesting a transition?
        max_r_idx = np.argmax(r_vals)
        if 0 < max_r_idx < len(r_vals)-1:
            print(f"  Peak ⟨r⟩ = {r_vals[max_r_idx]:.4f} at U = {U_vals[max_r_idx]:.1f}")
            print(f"  → Possible enhanced thermalization at this coupling")
        
        # Interacting non-integrability check
        interacting = [r for r in flux_results if r['U'] > 0.5]
        if interacting:
            max_r_inter = max(r['mean_r'] for r in interacting)
            mean_r_inter = np.mean([r['mean_r'] for r in interacting])
            print(f"\n  INTERACTING SECTOR (U > 0):")
            print(f"    max ⟨r⟩ = {max_r_inter:.4f}  mean ⟨r⟩ = {mean_r_inter:.4f}")
            
            if mean_r_inter > 0.55:
                print(f"    → GUE-like: T-broken AND thermalizing")
                print(f"    → SUPPORTS cycling (approach to equilibrium with broken T)")
            elif mean_r_inter > 0.50:
                print(f"    → Between GOE and GUE: thermalizing, T-breaking unclear")
                print(f"    → Need 6-site to resolve.")
            elif mean_r_inter > 0.45:
                print(f"    → GOE-like: thermalizing but T-breaking may be washed out")
            else:
                print(f"    → Poisson-like: possible integrable structure or sector mixing")
                print(f"    → Need symmetry-resolved analysis before concluding")
    
    # Save text summary
    with open(os.path.join(outdir, f'summary_N{nsites}.txt'), 'w') as f:
        f.write(f"K₄ D1-MB v2 (CORRECTED): N={nsites}, dim={dim}\n")
        f.write(f"Bond deduplication fix applied.\n\n")
        for r in results:
            f.write(f"U={r['U']:.1f} flux={r['flux']} <r>={r['mean_r']:.4f} "
                    f"Im/Re={r['im_re']:.4f} -> {r['verdict']}\n")
    
    print(f"\n  Saved to {outdir}/")
    return results


# ============================================================
# ENTRY POINT
# ============================================================

if __name__ == "__main__":
    p = argparse.ArgumentParser(description='K₄ Many-Body D1 Computation v2 (CORRECTED)')
    p.add_argument('--nsites', type=int, default=4)
    p.add_argument('--U', type=float, nargs='+', default=[0, 1, 2, 3, 4, 5, 6, 7, 8])
    p.add_argument('--nev', type=int, default=200)
    p.add_argument('--full', action='store_true')
    p.add_argument('--noflux', action='store_true')
    p.add_argument('--no-ee', action='store_true', help='Skip entanglement entropy')
    p.add_argument('--outdir', type=str, default='d1mb_results_v2')
    p.add_argument('--skip-selftest', action='store_true', help='Skip single-particle self-test')
    args = p.parse_args()
    
    # Always run self-test first (unless explicitly skipped)
    if not args.skip_selftest:
        verify_fix()
    
    run_analysis(
        nsites=args.nsites, U_values=args.U, nev=args.nev,
        full_diag=args.full, run_noflux=args.noflux,
        compute_ee=not args.no_ee, outdir=args.outdir,
    )
