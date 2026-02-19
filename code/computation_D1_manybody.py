#!/usr/bin/env python3
"""
COMPUTATION D1-MB: Many-Body Kâ‚„ Model â€” Interacting Case
=========================================================

Purpose: Determine whether the interacting Kâ‚„ boundary theory has 
anomalous thermalization properties that could support cosmological cycling.

What this computes:
  1. Many-body energy spectrum E_n(U) for the Kâ‚„ Hubbard model
  2. Level statistics r-ratio (GUE vs Poisson â†’ thermalization vs localization)
  3. Spectral asymmetry in the many-body spectrum (T-breaking signature)
  4. Entanglement entropy of half-system at various energies
  5. Comparison WITH vs WITHOUT â„¤â‚ƒ flux

System sizes:
  3 sites:  C(12,6)  =      924 states  [seconds]
  4 sites:  C(16,8)  =   12,870 states  [seconds-minutes]
  6 sites:  C(24,12) = 2,704,156 states [hours, sparse only]

Usage:
  python computation_D1_manybody.py --nsites 3 --U 0 2 4 6 8 --full --noflux
  python computation_D1_manybody.py --nsites 4 --U 0 2 4 6 8 --full --noflux
  python computation_D1_manybody.py --nsites 6 --U 0 4 8 --nev 500

Flags:
  --nsites N    Number of lattice sites (3, 4, or 6)
  --U val ...   List of U/t values to scan
  --nev N       Number of eigenvalues for sparse solver (default: 200)
  --full        Use full (dense) diagonalization (only for nsites â‰¤ 4)
  --noflux      Also run without â„¤â‚ƒ flux for comparison
  --outdir DIR  Output directory (default: ./d1mb_results/)
  --no-ee       Skip entanglement entropy (faster)

KEY CAVEATS for interpreting results:
  - At U=0 the system is integrable â†’ Poisson expected regardless of topology
  - Small systems have symmetry sectors that contaminate level statistics.
    The 3-site results have large error bars. 4-site is the minimum for
    meaningful level statistics. 6-site is where things get reliable.
  - Must eventually resolve Vâ‚„ Ã— translation symmetry sectors for clean ETH test.
    This code does NOT do that yet â€” it's the full Hilbert space.

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

# Kâ‚„ matching matrices (4Ã—4)
M1 = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]], dtype=complex)
M2 = np.array([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]], dtype=complex)
M3 = np.array([[0,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,0]], dtype=complex)

a1 = np.array([1.0, 0.0])
a2 = np.array([-0.5, np.sqrt(3)/2])


def build_triangular_lattice(nsites):
    """Build triangular lattice with PBC. Returns positions, neighbor_list."""
    if nsites == 3:
        positions = [np.array([0,0]), a1, a2]
        neighbors = [
            (0, 1, 0),  # Î´â‚ direction â†’ Mâ‚, phase 1
            (0, 2, 1),  # Î´â‚‚ direction â†’ Mâ‚‚, phase Ï‰
            (1, 2, 2),  # Î´â‚ƒ direction â†’ Mâ‚ƒ, phase Ï‰Â²
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
    
    neighbors = []
    for iy in range(Ly):
        for ix in range(Lx):
            i = site_index[(ix, iy)]
            # Direction 0 (Î´â‚): Mâ‚, phase 1
            j = site_index[((ix+1)%Lx, iy)]
            neighbors.append((i, j, 0))
            # Direction 1 (Î´â‚‚): Mâ‚‚, phase Ï‰
            j = site_index[(ix, (iy+1)%Ly)]
            neighbors.append((i, j, 1))
            # Direction 2 (Î´â‚ƒ): Mâ‚ƒ, phase Ï‰Â²
            j = site_index[((ix-1)%Lx, (iy-1)%Ly)]
            neighbors.append((i, j, 2))
    
    return positions, neighbors


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
    """Fock state â†’ basis index."""
    return {int(state): idx for idx, state in enumerate(basis)}


def sp_index(site, orbital, n_orb=4):
    """(site, orbital) â†’ single-particle index."""
    return site * n_orb + orbital


# ============================================================
# HAMILTONIAN CONSTRUCTION
# ============================================================

def build_hamiltonian(nsites, U, use_flux=True, basis=None, state_map=None):
    """
    Build the many-body Kâ‚„ Hubbard Hamiltonian as a sparse matrix.
    
    H = Hâ‚€ + UÂ·H_int
    Hâ‚€ = Î£_bonds t_{Î±Î²} câ€ _{iÎ±} c_{jÎ²} + h.c.
    H_int = U Î£_i Î£_{Î±<Î²} n_{iÎ±} n_{iÎ²}
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
                    
                    # câ€ _{iÎ±} c_{jÎ²}: need sp_j occupied, sp_i empty
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
                    
                    # h.c.: câ€ _{jÎ²} c_{iÎ±}: need sp_i occupied, sp_j empty
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
    GUE (T-broken, thermalizing):   âŸ¨râŸ© â‰ˆ 0.5996
    GOE (T-symmetric, thermalizing): âŸ¨râŸ© â‰ˆ 0.5307  
    Poisson (non-thermalizing):      âŸ¨râŸ© â‰ˆ 0.3863
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
                 run_noflux=False, compute_ee=True, outdir='d1mb_results'):
    
    os.makedirs(outdir, exist_ok=True)
    
    n_orb = 4
    n_sp = nsites * n_orb
    n_particles = nsites * 2
    dim = int(comb(n_sp, n_particles, exact=True))
    
    print("=" * 70)
    print(f"  Kâ‚„ MANY-BODY COMPUTATION D1-MB")
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
            flux_label = "WITH â„¤â‚ƒ flux" if use_flux else "WITHOUT flux (control)"
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
            print(f"    âŸ¨râŸ© = {mean_r:.4f} Â± {std_r:.4f}")
            print(f"    GUE={r_GUE:.4f}  GOE={r_GOE:.4f}  Poisson={r_Poi:.4f}")
            print(f"    â†’ {verdict}")
            
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
    print(f"  SUMMARY")
    print(f"{'='*70}")
    print(f"  {'U/t':>5s}  {'Flux':>4s}  {'âŸ¨râŸ©':>7s}  {'Â±':>6s}  "
          f"{'Im/Re':>7s}  {'Verdict':>10s}  {'S_gs':>6s}  {'S_mid':>6s}")
    print(f"  {'-'*5}  {'-'*4}  {'-'*7}  {'-'*6}  {'-'*7}  {'-'*10}  {'-'*6}  {'-'*6}")
    
    for r in results:
        f_str = "â„¤â‚ƒ" if r['flux'] else "â€”"
        s1 = f"{r['S_gs']:.3f}" if r['S_gs'] is not None else "â€”"
        s2 = f"{r['S_mid']:.3f}" if r['S_mid'] is not None else "â€”"
        print(f"  {r['U']:5.1f}  {f_str:>4s}  {r['mean_r']:7.4f}  {r['std_r']:6.4f}  "
              f"{r['im_re']:7.4f}  {r['verdict']:>10s}  {s1:>6s}  {s2:>6s}")
    
    # Key physics
    flux_results = [r for r in results if r['flux']]
    if len(flux_results) > 1:
        r_vals = [r['mean_r'] for r in flux_results]
        print(f"\n  Î”âŸ¨râŸ© across U scan: {max(r_vals)-min(r_vals):.4f}")
        
        # Does Im/Re decrease with U?
        im_vals = [r['im_re'] for r in flux_results]
        if im_vals[0] > im_vals[-1]:
            print(f"  Im/Re ratio DECREASES with U: {im_vals[0]:.4f} â†’ {im_vals[-1]:.4f}")
            print(f"  â†’ Interaction DILUTES the â„¤â‚ƒ T-breaking effect")
        
        # Is there a peak in âŸ¨râŸ© suggesting a transition?
        max_r_idx = np.argmax(r_vals)
        if 0 < max_r_idx < len(r_vals)-1:
            print(f"  Peak âŸ¨râŸ© = {r_vals[max_r_idx]:.4f} at U = {flux_results[max_r_idx]['U']:.1f}")
            print(f"  â†’ Possible enhanced thermalization at this coupling")
    
    # Save text summary
    with open(os.path.join(outdir, f'summary_N{nsites}.txt'), 'w') as f:
        f.write(f"Kâ‚„ D1-MB: N={nsites}, dim={dim}\n\n")
        for r in results:
            f.write(f"U={r['U']:.1f} flux={r['flux']} <r>={r['mean_r']:.4f} "
                    f"Im/Re={r['im_re']:.4f} â†’ {r['verdict']}\n")
    
    print(f"\n  Saved to {outdir}/")
    return results


# ============================================================
# ENTRY POINT
# ============================================================

if __name__ == "__main__":
    p = argparse.ArgumentParser(description='Kâ‚„ Many-Body D1 Computation')
    p.add_argument('--nsites', type=int, default=4)
    p.add_argument('--U', type=float, nargs='+', default=[0, 1, 2, 4, 6, 8])
    p.add_argument('--nev', type=int, default=200)
    p.add_argument('--full', action='store_true')
    p.add_argument('--noflux', action='store_true')
    p.add_argument('--no-ee', action='store_true', help='Skip entanglement entropy')
    p.add_argument('--outdir', type=str, default='d1mb_results')
    args = p.parse_args()
    
    run_analysis(
        nsites=args.nsites, U_values=args.U, nev=args.nev,
        full_diag=args.full, run_noflux=args.noflux,
        compute_ee=not args.no_ee, outdir=args.outdir,
    )
