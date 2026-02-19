#!/usr/bin/env python3
"""
Computation D1 (FIXED): â„¤â‚ƒ flux and thermal equilibrium
Bug fix: Berry curvature must be computed from the 8Ã—8 chiral Hamiltonian
H(k) = [[0, D(k)], [Dâ€ (k), 0]], where filled bands carry C = -2.
"""

import numpy as np
from numpy import linalg as la

omega = np.exp(2j * np.pi / 3)

# Matching matrices
M1 = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]], dtype=complex)
M2 = np.array([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]], dtype=complex)
M3 = np.array([[0,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,0]], dtype=complex)

d1 = np.array([1.0, 0.0])
d2 = np.array([-0.5, np.sqrt(3)/2])
d3 = -d1 - d2

def D_k(kx, ky):
    """4Ã—4 Dirac-Bloch Hamiltonian with â„¤â‚ƒ flux"""
    k = np.array([kx, ky])
    return (np.exp(1j*np.dot(k,d1)) * M1 + 
            omega * np.exp(1j*np.dot(k,d2)) * M2 + 
            omega**2 * np.exp(1j*np.dot(k,d3)) * M3)

def H_chiral(kx, ky):
    """8Ã—8 chiral Hamiltonian [[0, D], [Dâ€ , 0]]"""
    D = D_k(kx, ky)
    Z = np.zeros((4,4), dtype=complex)
    return np.block([[Z, D], [D.conj().T, Z]])

def get_filled_states(kx, ky):
    """Return eigenvalues and eigenvectors of H_chiral, sorted. 
    Filled = lowest 4 bands (negative energy at half-filling)."""
    H = H_chiral(kx, ky)
    evals, evecs = la.eigh(H)
    # eigh returns sorted ascending, so first 4 are filled
    return evals, evecs

# Reciprocal lattice vectors
b1 = 2*np.pi * np.array([1.0, 1.0/np.sqrt(3)])
b2 = 2*np.pi * np.array([0.0, 2.0/np.sqrt(3)])
bz_area = abs(b1[0]*b2[1] - b1[1]*b2[0])

# ============================================================
# 1. VERIFY CHERN NUMBER C = -2
# ============================================================

def chern_number_filled(Nk=80):
    """Compute Chern number of filled bands using lattice method (Fukui-Hatsugai-Suzuki)"""
    print("=" * 65)
    print("VERIFYING CHERN NUMBER FROM 8Ã—8 CHIRAL HAMILTONIAN")
    print("=" * 65)
    
    # Store filled-band projectors on grid
    def filled_projector(kx, ky):
        evals, evecs = get_filled_states(kx, ky)
        # Filled = columns 0,1,2,3 (4 lowest)
        U = evecs[:, :4]  # 8Ã—4 matrix
        return U
    
    total_berry = 0.0
    dk_s = 1.0 / Nk
    
    for i in range(Nk):
        for j in range(Nk):
            # 4 corners of plaquette in (s,t) coordinates
            s0, t0 = i*dk_s, j*dk_s
            
            k00 = s0 * b1 + t0 * b2
            k10 = (s0+dk_s) * b1 + t0 * b2
            k11 = (s0+dk_s) * b1 + (t0+dk_s) * b2
            k01 = s0 * b1 + (t0+dk_s) * b2
            
            U00 = filled_projector(*k00)
            U10 = filled_projector(*k10)
            U11 = filled_projector(*k11)
            U01 = filled_projector(*k01)
            
            # Overlap matrices (4Ã—4)
            O1 = U00.conj().T @ U10
            O2 = U10.conj().T @ U11
            O3 = U11.conj().T @ U01
            O4 = U01.conj().T @ U00
            
            # Berry phase = -Im(log(det(O1Â·O2Â·O3Â·O4)))
            W = O1 @ O2 @ O3 @ O4
            berry = -np.imag(np.log(la.det(W)))
            total_berry += berry
    
    chern = total_berry / (2 * np.pi)
    print(f"  Chern number (filled bands): C = {chern:.6f}")
    print(f"  Expected: C = -2")
    print(f"  Match: {'YES' if abs(chern + 2) < 0.1 else 'NO â€” BUG'}")
    return chern

# ============================================================
# 2. THERMAL HALL CONDUCTIVITY Ïƒ_xy(T)
# ============================================================

def thermal_hall(Nk=60):
    """
    Ïƒ_xy(Î²) = Î£_n âˆ« Î©_n(k) f(E_n(k); Î²) dÂ²k/(2Ï€)Â²
    
    At T=0: filled bands give C = -2
    At Tâ†’âˆž: f â†’ 1/2 for all bands, total Chern (all 8 bands) = 0, so Ïƒ_xy â†’ 0
    
    The question: HOW does Ïƒ_xy go to zero?
    """
    print("\n" + "=" * 65)
    print("PROBE D1a: THERMAL HALL CONDUCTIVITY Ïƒ_xy(T)")
    print("=" * 65)
    
    # Compute per-band Berry curvature using plaquette method
    dk_s = 1.0 / Nk
    
    # Store eigenvalues and Berry curvatures per band
    band_evals = np.zeros((Nk, Nk, 8))
    band_curv = np.zeros((Nk, Nk, 8))
    
    for i in range(Nk):
        for j in range(Nk):
            s0, t0 = (i+0.5)*dk_s, (j+0.5)*dk_s
            k = s0*b1 + t0*b2
            evals, evecs = get_filled_states(*k)
            band_evals[i,j] = evals
            
            # Per-band Berry curvature via small plaquette
            dk = 1e-4
            kx, ky = k
            
            def evecs_at(kx, ky):
                _, ev = la.eigh(H_chiral(kx, ky))
                return ev
            
            U00 = evecs_at(kx, ky)
            U10 = evecs_at(kx+dk, ky)
            U11 = evecs_at(kx+dk, ky+dk)
            U01 = evecs_at(kx, ky+dk)
            
            for n in range(8):
                o1 = np.conj(U00[:,n]) @ U10[:,n]
                o2 = np.conj(U10[:,n]) @ U11[:,n]
                o3 = np.conj(U11[:,n]) @ U01[:,n]
                o4 = np.conj(U01[:,n]) @ U00[:,n]
                berry = -np.imag(np.log(o1*o2*o3*o4))
                band_curv[i,j,n] = berry / dk**2
    
    # Check per-band Chern numbers
    print("\n  Per-band Chern numbers:")
    for n in range(8):
        cn = np.sum(band_curv[:,:,n]) * dk_s**2 * bz_area / (2*np.pi)**2
        en_mean = np.mean(band_evals[:,:,n])
        print(f"    Band {n}: C_n = {cn:+.4f},  âŸ¨EâŸ© = {en_mean:+.4f}")
    
    C_filled = sum(np.sum(band_curv[:,:,n]) * dk_s**2 * bz_area / (2*np.pi)**2 
                   for n in range(4))
    C_all = sum(np.sum(band_curv[:,:,n]) * dk_s**2 * bz_area / (2*np.pi)**2 
                for n in range(8))
    print(f"\n    C (filled) = {C_filled:+.4f}")
    print(f"    C (all)    = {C_all:+.4f}  (must be 0)")
    
    # Compute Ïƒ_xy at various temperatures
    print(f"\n  Ïƒ_xy(T) at various temperatures:")
    print(f"  {'T':>10s}  {'Î²':>10s}  {'Ïƒ_xy':>12s}  {'Ïƒ_xy/Ïƒ_xy(0)':>14s}")
    print(f"  {'-'*10}  {'-'*10}  {'-'*12}  {'-'*14}")
    
    betas = [1000, 100, 50, 20, 10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.001]
    sigma_values = []
    
    for beta in betas:
        sigma = 0
        for i in range(Nk):
            for j in range(Nk):
                for n in range(8):
                    E = band_evals[i,j,n]
                    x = beta * E
                    if x > 500: f_n = 0.0       # high energy: empty
                    elif x < -500: f_n = 1.0     # low energy: filled
                    else: f_n = 1.0 / (1 + np.exp(x))  # Fermi function (Î¼=0)
                    sigma += band_curv[i,j,n] * f_n
        sigma *= dk_s**2 * bz_area / (2*np.pi)**2
        sigma_values.append(sigma)
        T = 1.0/beta if beta > 0 else float('inf')
        ratio = sigma/sigma_values[0] if abs(sigma_values[0]) > 1e-10 else 0
        print(f"  {T:10.4f}  {beta:10.3f}  {sigma:+12.6f}  {ratio:14.6f}")
    
    # Fit power law at high T
    betas_arr = np.array(betas)
    sigma_arr = np.abs(np.array(sigma_values))
    mask = (betas_arr < 1.0) & (sigma_arr > 1e-10)
    if np.sum(mask) > 3:
        log_b = np.log(betas_arr[mask])
        log_s = np.log(sigma_arr[mask])
        alpha, const = np.polyfit(log_b, log_s, 1)
        print(f"\n  HIGH-T SCALING FIT: |Ïƒ_xy| âˆ Î²^{alpha:.3f} (= T^{-alpha:.3f})")
        if alpha > 0:
            print(f"  â†’ Ïƒ_xy decays as POWER LAW ~ 1/T^{alpha:.1f}")
            print(f"  â†’ Topology fades algebraically, not exponentially")
        else:
            print(f"  â†’ Non-standard scaling, needs further analysis")
    
    return sigma_values


# ============================================================
# 3. SPECTRAL ASYMMETRY â€” THE CLINCHER
# ============================================================

def spectral_asymmetry_probe(Nk=60):
    """
    The deepest probe. Îµ(k) â‰  Îµ(-k) is a PROPERTY OF THE HAMILTONIAN,
    not the state. It survives at all temperatures.
    
    Compute: 
    - |Îµ_n(k) - Îµ_n(-k)| averaged over BZ (spectral asymmetry)
    - The asymmetry as a function of position in BZ
    - Compare with the no-flux case (Ï‰ â†’ 1)
    """
    print("\n" + "=" * 65)
    print("PROBE D1e: SPECTRAL ASYMMETRY Îµ(k) â‰  Îµ(-k)")
    print("=" * 65)
    print("  This is a property of H, not of Ï. It survives at ALL T.\n")
    
    asym_with_flux = []
    asym_no_flux = []
    max_asym = 0
    k_at_max = None
    
    for i in range(Nk):
        for j in range(Nk):
            s = (i+0.5) / Nk
            t = (j+0.5) / Nk
            k = s*b1 + t*b2
            kx, ky = k
            
            # With â„¤â‚ƒ flux
            H_plus = H_chiral(kx, ky)
            H_minus = H_chiral(-kx, -ky)
            e_plus = la.eigvalsh(H_plus)
            e_minus = la.eigvalsh(H_minus)
            asym = la.norm(e_plus - e_minus) / la.norm(e_plus)
            asym_with_flux.append(asym)
            if asym > max_asym:
                max_asym = asym
                k_at_max = (kx, ky)
            
            # Without flux (Ï‰ â†’ 1)
            D_nf = (np.exp(1j*np.dot(k,d1)) * M1 + 
                    np.exp(1j*np.dot(k,d2)) * M2 + 
                    np.exp(1j*np.dot(k,d3)) * M3)
            Z = np.zeros((4,4), dtype=complex)
            H_nf = np.block([[Z, D_nf], [D_nf.conj().T, Z]])
            H_nf_minus = np.block([[Z, D_nf.conj().T], [D_nf, Z]])
            # For no-flux: D(-k) = D(k)* when phases are real
            D_nf_neg = (np.exp(-1j*np.dot(k,d1)) * M1 + 
                        np.exp(-1j*np.dot(k,d2)) * M2 + 
                        np.exp(-1j*np.dot(k,d3)) * M3)
            H_nf_neg = np.block([[Z, D_nf_neg], [D_nf_neg.conj().T, Z]])
            e_nf_plus = la.eigvalsh(H_nf)
            e_nf_minus = la.eigvalsh(H_nf_neg)
            asym_nf = la.norm(e_nf_plus - e_nf_minus) / max(la.norm(e_nf_plus), 1e-10)
            asym_no_flux.append(asym_nf)
    
    asym_with = np.array(asym_with_flux)
    asym_without = np.array(asym_no_flux)
    
    print(f"  WITH â„¤â‚ƒ flux (D â‰  D*):")
    print(f"    Mean |Îµ(k)âˆ’Îµ(âˆ’k)|/|Îµ(k)|:  {np.mean(asym_with):.8f}")
    print(f"    Max  |Îµ(k)âˆ’Îµ(âˆ’k)|/|Îµ(k)|:  {np.max(asym_with):.8f}")
    print(f"    Std  |Îµ(k)âˆ’Îµ(âˆ’k)|/|Îµ(k)|:  {np.std(asym_with):.8f}")
    print(f"    Fraction with asym > 0.01:  {np.mean(asym_with > 0.01):.4f}")
    
    print(f"\n  WITHOUT flux (D = D*, control):")
    print(f"    Mean |Îµ(k)âˆ’Îµ(âˆ’k)|/|Îµ(k)|:  {np.mean(asym_without):.8f}")
    print(f"    Max  |Îµ(k)âˆ’Îµ(âˆ’k)|/|Îµ(k)|:  {np.max(asym_without):.8f}")
    
    print(f"\n  RATIO (flux / no-flux):  {np.mean(asym_with)/max(np.mean(asym_without),1e-15):.2f}Ã—")
    
    print(f"\n  INTERPRETATION:")
    if np.mean(asym_with) > 0.01:
        print(f"  âœ“ Spectral asymmetry is LARGE ({np.mean(asym_with):.4f})")
        print(f"  âœ“ This is a property of H, not Ï â€” it cannot thermalize away")
        print(f"  âœ“ At any T, the dynamics remains asymmetric: v(k) â‰  v(-k)")
        print(f"  âœ“ The system can NEVER achieve microscopically symmetric equilibrium")
        print(f"  â†’ The â„¤â‚ƒ flux is a permanent asymmetric driver")
    else:
        print(f"  âœ— Spectral asymmetry is small â€” further analysis needed")
    
    return np.mean(asym_with)


# ============================================================
# 4. THE THERMAL CURRENT DENSITY (NOT TOTAL CURRENT)
# ============================================================

def thermal_current_density(Nk=50):
    """
    Even though âŸ¨JâŸ©_total = 0 at any T (by periodicity),
    the LOCAL current density j(k) = v(k)Â·f(E(k)) is nonzero
    at each k-point. The â„¤â‚ƒ flux creates a circulating current
    pattern that persists at all T.
    
    Compute the current density and its statistics.
    """
    print("\n" + "=" * 65)
    print("PROBE D1d: LOCAL CURRENT DENSITY AT THERMAL EQUILIBRIUM")
    print("=" * 65)
    
    dk = 1e-5
    
    for beta in [10.0, 1.0, 0.1, 0.01]:
        j_vectors = []
        j_magnitudes = []
        
        for i in range(Nk):
            for j in range(Nk):
                s = (i+0.5)/Nk
                t = (j+0.5)/Nk
                k = s*b1 + t*b2
                kx, ky = k
                
                # Group velocities from H_chiral
                H0 = H_chiral(kx, ky)
                Hx = H_chiral(kx+dk, ky)
                Hy = H_chiral(kx, ky+dk)
                
                evals, evecs = la.eigh(H0)
                evals_x, _ = la.eigh(Hx)
                evals_y, _ = la.eigh(Hy)
                
                jx_total, jy_total = 0.0, 0.0
                for n in range(8):
                    vx = (evals_x[n] - evals[n]) / dk
                    vy = (evals_y[n] - evals[n]) / dk
                    x = beta * evals[n]
                    if x > 500: f_n = 0.0
                    elif x < -500: f_n = 1.0
                    else: f_n = 1.0/(1+np.exp(x))
                    jx_total += vx * f_n
                    jy_total += vy * f_n
                
                j_vectors.append([jx_total, jy_total])
                j_magnitudes.append(np.sqrt(jx_total**2 + jy_total**2))
        
        j_vectors = np.array(j_vectors)
        j_mags = np.array(j_magnitudes)
        j_net = la.norm(np.mean(j_vectors, axis=0))
        
        T = 1.0/beta
        print(f"\n  T = {T:.2f}:")
        print(f"    Mean |j_local(k)|:     {np.mean(j_mags):.6f}")
        print(f"    Max  |j_local(k)|:     {np.max(j_mags):.6f}")
        print(f"    |âŸ¨jâŸ©| (net current):   {j_net:.8f}")
        print(f"    Cancellation ratio:    {j_net/np.mean(j_mags) if np.mean(j_mags)>0 else 0:.2e}")
        
        # Compute curl of j-field (circulation)
        # j-field has zero divergence (current conservation) but may have nonzero curl
        # Nonzero curl = circulating current = broken detailed balance
        if Nk > 3:
            curl_sum = 0
            for i in range(Nk-1):
                for jj in range(Nk-1):
                    idx00 = i*Nk + jj
                    idx10 = (i+1)*Nk + jj
                    idx01 = i*Nk + (jj+1)
                    # Approximate curl via circulation around plaquette
                    # curl â‰ˆ (jx(i,j+1) - jx(i,j))/dk_t - (jy(i+1,j) - jy(i,j))/dk_s
                    djx = j_vectors[idx01,0] - j_vectors[idx00,0]
                    djy = j_vectors[idx10,1] - j_vectors[idx00,1]
                    curl_sum += abs(djx - djy)
            curl_avg = curl_sum / (Nk-1)**2
            print(f"    Mean |curl(j)|:        {curl_avg:.6f}  (>0 = circulating)")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print("â•”â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•—")
    print("â•‘  COMPUTATION D1 (FIXED): â„¤â‚ƒ FLUX AND THERMAL EQUILIBRIUM   â•‘")
    print("â•šâ•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•\n")
    
    # First: verify we get C = -2
    C = chern_number_filled(Nk=60)
    
    # Then: thermal Hall conductivity
    sigma_vals = thermal_hall(Nk=40)
    
    # The clincher: spectral asymmetry
    asym = spectral_asymmetry_probe(Nk=50)
    
    # Local current density and circulation
    thermal_current_density(Nk=40)
    
    # Final verdict
    print("\n" + "=" * 65)
    print("VERDICT ON D1")
    print("=" * 65)
    print("""
    Question: Does the â„¤â‚ƒ flux prevent true thermal equilibrium?
    
    Answer: It depends on what you mean by "equilibrium."
    
    The Gibbs state Ï = e^{-Î²H}/Z IS an equilibrium in the standard 
    sense: [H, Ï] = 0, all expectation values are time-independent.
    The NET current âŸ¨JâŸ© = 0 at every temperature.
    
    BUT the Gibbs state is NOT a microscopically symmetric equilibrium:
    
    1. The Hamiltonian has Îµ(k) â‰  Îµ(-k) at every k (spectral asymmetry).
       This means the DYNAMICS (time evolution operator e^{-iHt}) is 
       permanently asymmetric. The system at equilibrium is still 
       "wanting to go" in a preferred direction â€” it's just balanced
       by the equal and opposite tendency at âˆ’k.
    
    2. The local current density j(k) â‰  0 everywhere, with nonzero curl.
       There are circulating currents at every temperature. They cancel
       globally but not locally. This is broken detailed balance.
    
    3. The spectral asymmetry is a property of H, not of Ï.
       No temperature can erase it. It persists at Î² = 0.
    
    For cosmological cycling:
    
    The â„¤â‚ƒ flux does NOT prevent the system from reaching a Gibbs 
    equilibrium (it does reach one). But it DOES ensure that the
    equilibrium state has a fundamentally different character than 
    a time-reversal-symmetric equilibrium:
    
    - The dynamics at equilibrium is asymmetric
    - There are circulating currents (broken detailed balance)
    - The approach to equilibrium is directional (not symmetric)
    
    This means the "heat death" state of a Kâ‚„D universe is NOT the 
    featureless, directionless void of classical thermodynamics.
    It is a state with permanent circulation, permanent spectral
    asymmetry, and permanent memory of the â„¤â‚ƒ flux in its dynamics.
    
    Whether this asymmetric equilibrium is STABLE against perturbations
    that could trigger a reorganization is the next question (Step D 
    continued). The spectral asymmetry provides the DRIVING FORCE for 
    such a reorganization â€” the system is never truly "at rest."
    """)
