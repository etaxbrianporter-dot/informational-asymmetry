"""
Kâ‚† Computational Tools (CORRECTED)
====================================
Verified against paper: aâ‚‚ = 3.3060, aâ‚„ = 4.0682, all 15 eigenvalues match.

Key corrections:
1. Gram matrix uses Re(Ï‰^dz), not complex Ï‰^dz
2. Quartic form uses Re(Ï‰^{dz_total}), not Re(Ï‰^{dz1}) Ã— Re(Ï‰^{dz2})
"""
import numpy as np
from scipy.linalg import eigh

omega = np.exp(2j * np.pi / 3)

# Sorted assignment from paper Appendix B
MATCHINGS_K6 = [
    {(0,1),(2,3),(4,5)}, {(0,1),(2,4),(3,5)}, {(0,1),(2,5),(3,4)},
    {(0,2),(1,3),(4,5)}, {(0,2),(1,4),(3,5)}, {(0,2),(1,5),(3,4)},
    {(0,3),(1,2),(4,5)}, {(0,3),(1,4),(2,5)}, {(0,3),(1,5),(2,4)},
    {(0,4),(1,2),(3,5)}, {(0,4),(1,3),(2,5)}, {(0,4),(1,5),(2,3)},
    {(0,5),(1,2),(3,4)}, {(0,5),(1,3),(2,4)}, {(0,5),(1,4),(2,3)},
]
DIR_K6 = [0,1,2, 0,1,2, 0,1,2, 0,1,2, 0,1,2]
Z3_K6  = [0,0,0, 1,1,1, 2,2,2, 1,1,1, 2,2,2]

def overlap_matrix(matchings):
    N = len(matchings)
    nv = max(max(max(e) for e in m) for m in matchings) + 1
    O = np.zeros((N,N))
    for i in range(N):
        for j in range(i,N):
            shared = len(matchings[i] & matchings[j])
            O[i,j] = O[j,i] = 2*shared if i!=j else nv
    return O

def gram_matrix(matchings, dirs, z3):
    """G_ij = Tr(M_i M_j) Â· Re(Ï‰^{z3_i - z3_j}) Â· Î´(dir_i, dir_j)"""
    O = overlap_matrix(matchings)
    N = len(matchings)
    G = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            if dirs[i] != dirs[j]: continue
            G[i,j] = O[i,j] * np.real(omega**((z3[i]-z3[j])%3))
    return G

def adjacency_matrices(matchings, nv):
    Ms = []
    for m in matchings:
        M = np.zeros((nv,nv))
        for (a,b) in m: M[a,b] = M[b,a] = 1
        Ms.append(M)
    return Ms

def quartic_form(v, Ms, dirs, z3):
    """aâ‚„ = H(v,v,v,v) with CORRECT phase: Re(Ï‰^{(z3_j-z3_i)+(z3_l-z3_k)})"""
    active = [i for i in range(len(v)) if abs(v[i]) > 1e-14]
    a4 = 0.0
    for i in active:
        for j in active:
            if dirs[j] != dirs[i]: continue
            MiMj = Ms[i] @ Ms[j]
            for k in active:
                if dirs[k] != dirs[i]: continue
                MiMjMk = MiMj @ Ms[k]
                for l in active:
                    if dirs[l] != dirs[i]: continue
                    tr = np.trace(MiMjMk @ Ms[l])
                    dz = (z3[j]-z3[i] + z3[l]-z3[k]) % 3
                    a4 += v[i]*v[j]*v[k]*v[l] * tr * np.real(omega**dz)
    return a4

def physical_projector(O, nv):
    """Projector onto V_{(2n-2,2)} from overlap matrix."""
    n = nv // 2
    def dfact(k): return k * dfact(k-2) if k > 0 else 1
    lam_max = 2*n*dfact(2*n-3)
    lam_mid = 4*(n-1)*dfact(2*n-5)
    N = O.shape[0]
    return O @ (O - lam_max*np.eye(N)) / (lam_mid*(lam_mid-lam_max)), lam_max, lam_mid

def physical_basis(Pi):
    evals, evecs = eigh(Pi)
    return evecs[:, np.abs(evals - 1.0) < 0.01]

def compute_k6():
    """Full Kâ‚† computation. Returns dict of all quantities."""
    G = gram_matrix(MATCHINGS_K6, DIR_K6, Z3_K6)
    O = overlap_matrix(MATCHINGS_K6)
    Ms = adjacency_matrices(MATCHINGS_K6, 6)
    
    evals, evecs = eigh(G)
    v0 = evecs[:,0]
    if v0[6] < 0: v0 = -v0
    
    a2 = evals[0]
    a4 = quartic_form(v0, Ms, DIR_K6, Z3_K6)
    
    Pi, lmax, lmid = physical_projector(O, 6)
    V_phys = physical_basis(Pi)
    G_phys = V_phys.T @ G @ V_phys
    evals_phys = np.sort(eigh(G_phys, eigvals_only=True))
    
    phys_frac = np.linalg.norm(Pi @ v0)**2
    
    return {
        'G': G, 'O': O, 'Ms': Ms,
        'evals': evals, 'evecs': evecs, 'v0': v0,
        'a2': a2, 'a4': a4, 'R': a4/a2**2,
        'Pi': Pi, 'V_phys': V_phys,
        'G_phys': G_phys, 'evals_phys': evals_phys,
        'phys_frac': phys_frac,
    }

if __name__ == "__main__":
    r = compute_k6()
    print(f"aâ‚‚ = {r['a2']:.10f}  âœ“" if abs(r['a2']-3.3060051829)<1e-6 else f"aâ‚‚ = {r['a2']:.10f}  âœ—")
    print(f"aâ‚„ = {r['a4']:.10f}  âœ“" if abs(r['a4']-4.0681621736)<1e-4 else f"aâ‚„ = {r['a4']:.10f}  âœ—")
    print(f"R  = {r['R']:.10f}  âœ“" if abs(r['R']-0.3722127085)<1e-6 else f"R  = {r['R']:.10f}  âœ—")
    print(f"\nPhysical sector: {r['V_phys'].shape[1]}D, vâ‚€ overlap = {r['phys_frac']:.4f}")
    mW = 80.379
    print(f"mH = {np.sqrt(2*r['a4']/r['a2'])*mW:.2f} GeV")
