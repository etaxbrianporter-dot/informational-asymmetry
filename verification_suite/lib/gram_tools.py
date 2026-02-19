"""
gram_tools.py — Gram matrix construction for K_{2n} on surfaces
================================================================
Handles direction assignments, Z₃ phases, and BZ averaging.
"""
import numpy as np
from scipy.linalg import eigh


def gram_matrix(matchings, dirs, z3, n_vertices):
    """
    Construct the Gram matrix G_{ij} = Tr(M_i M_j) · Re(ω^{z3_j - z3_i}) · δ(dir_i, dir_j)
    
    Parameters:
        matchings: list of sets of edges
        dirs: list of direction indices (one per matching)
        z3: list of Z₃ phase indices (one per matching)
        n_vertices: number of vertices
    
    Returns:
        G: N×N Gram matrix
    """
    omega = np.exp(2j * np.pi / 3)
    N = len(matchings)
    
    # Compute overlap matrix
    O = np.zeros((N, N))
    for i in range(N):
        for j in range(i, N):
            if i == j:
                O[i, j] = n_vertices
            else:
                shared = len(matchings[i] & matchings[j])
                O[i, j] = O[j, i] = 2 * shared
    
    # Apply direction filter and Z₃ phases
    G = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if dirs[i] != dirs[j]:
                continue
            dz = (z3[j] - z3[i]) % 3
            G[i, j] = O[i, j] * np.real(omega ** dz)
    
    return G


def quartic_form(v, matchings, dirs, z3, n_vertices):
    """
    Compute a₄ = H(v,v,v,v) the quartic form.
    Uses CORRECT phase: Re(ω^{(z3_j-z3_i) + (z3_l-z3_k)})
    """
    omega = np.exp(2j * np.pi / 3)
    N = len(v)
    
    # Build adjacency matrices
    Ms = []
    for m in matchings:
        M = np.zeros((n_vertices, n_vertices))
        for edge in m:
            a, b = tuple(sorted(edge))
            M[a, b] = M[b, a] = 1
        Ms.append(M)
    
    active = [i for i in range(N) if abs(v[i]) > 1e-14]
    a4 = 0.0
    for i in active:
        for j in active:
            if dirs[j] != dirs[i]:
                continue
            MiMj = Ms[i] @ Ms[j]
            for k in active:
                if dirs[k] != dirs[i]:
                    continue
                MiMjMk = MiMj @ Ms[k]
                for l in active:
                    if dirs[l] != dirs[i]:
                        continue
                    tr = np.trace(MiMjMk @ Ms[l])
                    dz = (z3[j] - z3[i] + z3[l] - z3[k]) % 3
                    a4 += v[i]*v[j]*v[k]*v[l] * tr * np.real(omega**dz)
    return a4


def diagonalize_gram(G):
    """
    Diagonalize the Gram matrix.
    Returns eigenvalues (sorted ascending) and eigenvectors.
    """
    evals, evecs = eigh(G)
    return evals, evecs


def vacuum_eigenvector(G):
    """
    Extract the vacuum (ground state) eigenvector.
    Returns (eigenvalue, eigenvector).
    """
    evals, evecs = eigh(G)
    v0 = evecs[:, 0]
    # Canonical sign convention
    max_idx = np.argmax(np.abs(v0))
    if v0[max_idx] < 0:
        v0 = -v0
    return evals[0], v0


# ========================================
# K₆ specific: sorted assignment from paper
# ========================================

K6_MATCHINGS = [
    frozenset({(0,1),(2,3),(4,5)}), frozenset({(0,1),(2,4),(3,5)}), frozenset({(0,1),(2,5),(3,4)}),
    frozenset({(0,2),(1,3),(4,5)}), frozenset({(0,2),(1,4),(3,5)}), frozenset({(0,2),(1,5),(3,4)}),
    frozenset({(0,3),(1,2),(4,5)}), frozenset({(0,3),(1,4),(2,5)}), frozenset({(0,3),(1,5),(2,4)}),
    frozenset({(0,4),(1,2),(3,5)}), frozenset({(0,4),(1,3),(2,5)}), frozenset({(0,4),(1,5),(2,3)}),
    frozenset({(0,5),(1,2),(3,4)}), frozenset({(0,5),(1,3),(2,4)}), frozenset({(0,5),(1,4),(2,3)}),
]

K6_DIRS = [0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2]
K6_Z3   = [0, 0, 0, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2]


def compute_k6_full():
    """
    Full K₆ computation returning all key quantities.
    """
    G = gram_matrix(K6_MATCHINGS, K6_DIRS, K6_Z3, 6)
    evals, evecs = eigh(G)
    v0 = evecs[:, 0]
    if v0[6] < 0:
        v0 = -v0
    
    a2 = evals[0]
    a4 = quartic_form(v0, K6_MATCHINGS, K6_DIRS, K6_Z3, 6)
    R = a4 / a2**2
    
    return {
        'G': G,
        'evals': evals,
        'v0': v0,
        'a2': a2,
        'a4': a4,
        'R': R,
    }
