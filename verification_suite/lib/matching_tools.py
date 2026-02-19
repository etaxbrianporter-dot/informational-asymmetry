"""
matching_tools.py — Perfect matching enumeration for K_{2n}
============================================================
No external dependencies beyond numpy.
All matchings are sets of frozensets of 2-element tuples.
"""
import numpy as np
from itertools import combinations
from math import factorial


def double_factorial(n):
    """Compute (2n-1)!! = 1·3·5·...·(2n-1)"""
    if n <= 0:
        return 1
    result = 1
    for k in range(1, 2*n, 2):
        result *= k
    return result


def enumerate_matchings(n_vertices):
    """
    Enumerate all perfect matchings of K_{n_vertices}.
    n_vertices must be even.
    Returns list of frozensets of 2-tuples.
    """
    assert n_vertices % 2 == 0, "Need even number of vertices"
    vertices = list(range(n_vertices))
    matchings = []
    _enumerate_recursive(vertices, frozenset(), matchings)
    return matchings


def _enumerate_recursive(remaining, current_matching, all_matchings):
    if len(remaining) == 0:
        all_matchings.append(current_matching)
        return
    first = remaining[0]
    rest = remaining[1:]
    for i, partner in enumerate(rest):
        edge = frozenset((first, partner))
        new_remaining = rest[:i] + rest[i+1:]
        _enumerate_recursive(new_remaining, current_matching | {edge}, all_matchings)


def matching_to_sorted_tuples(matching):
    """Convert frozenset matching to sorted list of sorted tuples."""
    return sorted(tuple(sorted(e)) for e in matching)


def overlap(m1, m2):
    """Number of shared edges between two matchings."""
    return len(m1 & m2)


def overlap_matrix(matchings):
    """
    Compute the overlap matrix O_{ij} = 2 * |M_i ∩ M_j|.
    Diagonal: O_{ii} = n_vertices (= 2n for K_{2n}).
    """
    N = len(matchings)
    if N == 0:
        return np.zeros((0, 0))
    # Determine n_vertices from first matching
    nv = 2 * len(list(matchings[0]))
    O = np.zeros((N, N))
    for i in range(N):
        for j in range(i, N):
            if i == j:
                O[i, j] = nv
            else:
                shared = overlap(matchings[i], matchings[j])
                O[i, j] = O[j, i] = 2 * shared
    return O


def overlap_distribution(matchings):
    """
    Count pairs with each overlap value.
    Returns dict: overlap_value -> count
    """
    N = len(matchings)
    dist = {}
    for i in range(N):
        for j in range(i+1, N):
            ov = 2 * overlap(matchings[i], matchings[j])
            dist[ov] = dist.get(ov, 0) + 1
    return dist


def incidence_matrix(matchings, n_vertices):
    """
    Construct the N × E matching-edge incidence matrix A.
    A_{i,e} = 1 if matching i contains edge e.
    Edges ordered lexicographically.
    """
    edges = [(a, b) for a in range(n_vertices) for b in range(a+1, n_vertices)]
    edge_idx = {e: i for i, e in enumerate(edges)}
    N = len(matchings)
    E = len(edges)
    A = np.zeros((N, E))
    for i, m in enumerate(matchings):
        for edge in m:
            e = tuple(sorted(edge))
            A[i, edge_idx[e]] = 1
    return A, edges


def adjacency_matrix_from_matching(matching, n_vertices):
    """Convert a matching to its n×n adjacency matrix."""
    M = np.zeros((n_vertices, n_vertices))
    for edge in matching:
        a, b = tuple(sorted(edge))
        M[a, b] = M[b, a] = 1
    return M
