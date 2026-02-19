"""
k4_tools.py — K₄ spectral action, Pfaffian mechanism, Hessian signature
========================================================================
Verifies Paper I claims: 13/15 Lorentzian, D² eigenvalues, s_crit.
"""
import numpy as np
from itertools import combinations


def enumerate_k4_subgraphs():
    """
    Enumerate all C(6,4) = 15 four-edge subgraphs of K₄.
    K₄ has 6 edges: (0,1),(0,2),(0,3),(1,2),(1,3),(2,3).
    Returns list of 4-edge subsets.
    """
    edges = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
    return list(combinations(edges, 4))


def build_D_matrix(subgraph, s=1.0, weights=None):
    """
    Build 4×4 skew-symmetric matrix D from a 4-edge subgraph.
    D_{ij} = s·w_e for edge e=(i,j) with i<j, D_{ji} = -D_{ij}.
    """
    if weights is None:
        weights = [1.0] * 4
    D = np.zeros((4, 4))
    for k, (i, j) in enumerate(subgraph):
        D[i, j] = s * weights[k]
        D[j, i] = -s * weights[k]
    return D


def pfaffian_4x4(D):
    """
    Pfaffian of a 4×4 skew-symmetric matrix.
    Pf(D) = D₁₂D₃₄ − D₁₃D₂₄ + D₁₄D₂₃
    """
    return D[0,1]*D[2,3] - D[0,2]*D[1,3] + D[0,3]*D[1,2]


def classify_subgraph(subgraph):
    """
    Classify a 4-edge subgraph of K₄.
    Returns: 'hub-spoke', 'sequential-HC', or 'scrambled-HC'
    """
    all_edges = {(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)}
    present = set(subgraph)
    absent = all_edges - present
    
    # Count vertex degrees
    degree = [0] * 4
    for (i, j) in present:
        degree[i] += 1
        degree[j] += 1
    
    # Hub-spoke: one vertex has degree 3
    if max(degree) == 3:
        return 'hub-spoke'
    
    # HC: all degrees 2 (forms a 4-cycle)
    # Sequential vs scrambled: check if matchings cancel
    D = build_D_matrix(subgraph, s=1.0)
    pf = pfaffian_4x4(D)
    if abs(pf) > 1e-10:
        return 'sequential-HC'
    else:
        return 'scrambled-HC'


def D_squared_eigenvalues(D):
    """Compute eigenvalues of D² (should all be ≤ 0 for skew-symmetric D)."""
    D2 = D @ D
    return np.sort(np.linalg.eigvalsh(D2))


def spectral_action(D):
    """I = Tr(exp(D²/2))"""
    from scipy.linalg import expm
    D2 = D @ D
    return np.trace(expm(D2 / 2))


def hessian_spectral_action(subgraph, s=1.0):
    """
    Compute the 4×4 Hessian ∂²I/∂w_a∂w_b of the spectral action
    at uniform weights w = (1,1,1,1).
    Uses finite differences (verified against analytic for hub-spokes).
    """
    h = 1e-5
    n_edges = 4
    H = np.zeros((n_edges, n_edges))
    
    w0 = np.ones(n_edges)
    I0 = spectral_action(build_D_matrix(subgraph, s, w0))
    
    for a in range(n_edges):
        for b in range(a, n_edges):
            wp = w0.copy(); wp[a] += h; wp[b] += h
            wm1 = w0.copy(); wm1[a] += h; wm1[b] -= h
            wm2 = w0.copy(); wm2[a] -= h; wm2[b] += h
            wmm = w0.copy(); wmm[a] -= h; wmm[b] -= h
            
            Ipp = spectral_action(build_D_matrix(subgraph, s, wp))
            Ipm = spectral_action(build_D_matrix(subgraph, s, wm1))
            Imp = spectral_action(build_D_matrix(subgraph, s, wm2))
            Imm = spectral_action(build_D_matrix(subgraph, s, wmm))
            
            H[a, b] = (Ipp - Ipm - Imp + Imm) / (4 * h**2)
            H[b, a] = H[a, b]
    
    return H


def hessian_signature(H, tol=1e-8):
    """Return (n_positive, n_negative) eigenvalue counts."""
    evals = np.linalg.eigvalsh(H)
    n_pos = np.sum(evals > tol)
    n_neg = np.sum(evals < -tol)
    return (n_pos, n_neg)


def k4_matchings():
    """The three perfect matchings of K₄."""
    return [
        frozenset({(0,1), (2,3)}),  # M₁
        frozenset({(0,2), (1,3)}),  # M₂
        frozenset({(0,3), (1,2)}),  # M₃
    ]


def matching_product(m1, m2, n_vertices=4):
    """
    Product of two matchings as permutations.
    Each matching is a fixed-point-free involution.
    Returns the composition as a permutation (list).
    """
    # Convert matchings to permutation form
    def to_perm(m):
        p = list(range(n_vertices))
        for edge in m:
            a, b = tuple(sorted(edge))
            p[a] = b
            p[b] = a
        return p
    
    p1 = to_perm(m1)
    p2 = to_perm(m2)
    # Composition: p1 ∘ p2
    return [p1[p2[i]] for i in range(n_vertices)]


def perm_to_matching(perm):
    """Convert an involution permutation back to a matching (frozenset of tuples)."""
    edges = set()
    for i in range(len(perm)):
        if perm[i] > i:
            edges.add((i, perm[i]))
    return frozenset(edges)


def heat_kernel_k4_internal():
    """
    K₄ internal spectral triple: eigenvalues λ_n = n, degeneracy D(n) = 2n.
    Chirality γ₅ assigns grading (-1)^{n+1}.
    
    Returns a₀ = Tr(1), b₀ = Tr(γ₅).
    """
    a0 = sum(2*n for n in range(1, 5))     # 2+4+6+8 = 20
    b0 = sum(2*n * (-1)**(n+1) for n in range(1, 5))  # 2-4+6-8 = -4
    return a0, b0
