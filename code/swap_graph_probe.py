# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Swap Graph Probe (Sense 5)
============================
Two matchings are neighbors if they differ by exactly one edge-swap:
remove two non-adjacent edges, reconnect the four freed vertices the
other way (there are exactly 2 reconnections; the one that isn't the
original gives the swap neighbor).

This probe builds the swap graph on matchings and asks:
  5a: Degree distribution — how many swap neighbors per matching?
  5b: Swap graph Laplacian spectrum — spectral gap, mixing time
  5c: Sign grading — does swapping ALWAYS flip the Pfaffian sign? (bipartite?)
  5d: Orbit connectivity — how many swaps cross orbit boundaries?
  5e: Swap × overlap correlation — are swap-neighbors high or low overlap?
  5f: Block structure — does the swap graph respect direction blocks?
  5g: Connected components — is the swap graph connected?

Usage:
  python swap_graph_probe.py [--levels 6 8 10 12]
"""

import numpy as np
from collections import defaultdict
import time
import argparse

# ======================================================================
# CORE FUNCTIONS (from coupling probe)
# ======================================================================

def gen_matchings(nv):
    verts = list(range(nv))
    matchings = []
    def _build(rem, cur):
        if not rem:
            matchings.append(frozenset(cur))
            return
        first = rem[0]
        for i in range(1, len(rem)):
            _build([v for j, v in enumerate(rem) if j != 0 and j != i],
                   cur + [(first, rem[i])])
    _build(verts, [])
    return matchings


def edge_dir(a, b, p):
    a, b = min(a, b), max(a, b)
    if a >= p or b >= p:
        return -1
    d = (b - a) % p
    return min(d, p - d) - 1


def torus_ms(matching, p):
    dirs = []
    for e in sorted(matching):
        d = edge_dir(e[0], e[1], p)
        if d >= 0:
            dirs.append(d)
    return tuple(sorted(dirs))


def compute_dir_perm(p, g=2):
    n_dirs = (p - 1) // 2
    def dc(d):
        d = d % p
        return min(d, p - d) - 1
    perm = {}
    for k in range(n_dirs):
        nd = dc((g * (k + 1)) % p)
        if nd < 0:
            return None
        perm[k] = nd
    return perm


def get_orbits(blocks, perm):
    visited = set()
    orbits = []
    for ms in sorted(blocks.keys()):
        if ms in visited:
            continue
        orbit_ms = []
        cur = ms
        while cur not in visited:
            visited.add(cur)
            orbit_ms.append(cur)
            cur = tuple(sorted(perm[d] for d in cur))
        orbits.append(orbit_ms)
    return orbits


def compress_spectrum(evals, tol=1e-6):
    groups = []
    i = 0
    while i < len(evals):
        v = evals[i]
        c = 1
        while i + c < len(evals) and abs(evals[i + c] - v) < tol:
            c += 1
        groups.append((v, c))
        i += c
    return groups


# ======================================================================
# PFAFFIAN SIGN (from pfaffian probe)
# ======================================================================

def perm_sign(perm):
    n = len(perm)
    visited = [False] * n
    n_cycles = 0
    for i in range(n):
        if visited[i]:
            continue
        n_cycles += 1
        j = i
        while not visited[j]:
            visited[j] = True
            j = perm[j]
    return (-1) ** (n - n_cycles)


def pfaffian_sign(matching, nv):
    edges = sorted(matching, key=lambda e: (min(e), max(e)))
    edges = [(min(a, b), max(a, b)) for a, b in edges]
    perm_list = []
    for a, b in edges:
        perm_list.append(a)
        perm_list.append(b)
    return perm_sign(perm_list)


# ======================================================================
# SWAP GRAPH CONSTRUCTION
# ======================================================================

def find_swap_neighbors(matching, nv, m_to_idx):
    """Find all swap neighbors of a matching.
    
    A swap takes two edges (a,b) and (c,d) with {a,b}∩{c,d}=∅
    and replaces them with either (a,c),(b,d) or (a,d),(b,c).
    Both alternatives are valid new matchings (if they differ from original).
    """
    edges = list(matching)
    neighbors = set()
    
    for i in range(len(edges)):
        for j in range(i + 1, len(edges)):
            e1 = edges[i]
            e2 = edges[j]
            a, b = min(e1), max(e1)
            c, d = min(e2), max(e2)
            
            # Vertices {a,b,c,d} must be disjoint
            if len({a, b, c, d}) != 4:
                continue
            
            # Two alternative pairings of {a,b,c,d}:
            # Original: (a,b), (c,d)
            # Alt 1: (a,c), (b,d)
            # Alt 2: (a,d), (b,c)
            remaining = matching - {e1, e2}
            
            for alt in [
                frozenset([(min(a,c), max(a,c)), (min(b,d), max(b,d))]),
                frozenset([(min(a,d), max(a,d)), (min(b,c), max(b,c))])
            ]:
                new_m = remaining | alt
                if new_m != matching and new_m in m_to_idx:
                    neighbors.add(m_to_idx[new_m])
    
    return neighbors


def build_swap_graph(raw_matchings, nv):
    """Build adjacency list for the swap graph."""
    m_to_idx = {m: i for i, m in enumerate(raw_matchings)}
    N = len(raw_matchings)
    
    adj = [set() for _ in range(N)]
    for i, m in enumerate(raw_matchings):
        nbrs = find_swap_neighbors(m, nv, m_to_idx)
        for j in nbrs:
            adj[i].add(j)
            adj[j].add(i)
    
    return adj


def swap_graph_laplacian(adj, N):
    """Build Laplacian matrix L = D - A."""
    L = np.zeros((N, N))
    for i in range(N):
        L[i, i] = len(adj[i])
        for j in adj[i]:
            L[i, j] = -1
    return L


# ======================================================================
# SENSE 5 PROBE
# ======================================================================

class SwapGraphProbe:
    def __init__(self, nv):
        self.nv = nv
        self.p = nv - 1
        self.n = nv // 2
        self.raw = None
        self.match_sets = None
        self.blocks = None
        self.orbits = None
        self.vac_orbit_idx = None
        self.fp_orbit_idx = None
        self.orbit_size = None
        self.signs = None
        self.adj = None
        self.results = {}
    
    def generate(self):
        """Generate matchings, orbits, signs, and swap graph."""
        t0 = time.time()
        self.raw = gen_matchings(self.nv)
        self.match_sets = [set(m) for m in self.raw]
        N = len(self.raw)
        
        # Signs
        self.signs = np.array([pfaffian_sign(m, self.nv) for m in self.raw], dtype=float)
        
        # Orbits
        self.blocks = defaultdict(list)
        for i, m in enumerate(self.raw):
            self.blocks[torus_ms(m, self.p)].append(i)
        
        perm = compute_dir_perm(self.p)
        raw_orbits = get_orbits(self.blocks, perm)
        
        self.orbits = []
        self.fp_orbit_idx = None
        for oi, oms in enumerate(raw_orbits):
            idx = []
            for ms in oms:
                idx.extend(self.blocks[ms])
            self.orbits.append((oms, idx))
            if len(oms) == 1:
                self.fp_orbit_idx = oi
        
        self.vac_orbit_idx = None
        for oi, (oms, idx) in enumerate(self.orbits):
            if len(self.blocks[oms[0]]) == 2 * self.p:
                self.vac_orbit_idx = oi
                break
        if self.vac_orbit_idx is None:
            for oi, (oms, idx) in enumerate(self.orbits):
                if len(oms) > 1:
                    self.vac_orbit_idx = oi
                    break
        
        if self.vac_orbit_idx is not None:
            self.orbit_size = len(self.orbits[self.vac_orbit_idx][0])
        
        # Build swap graph
        t1 = time.time()
        self.adj = build_swap_graph(self.raw, self.nv)
        t2 = time.time()
        
        total_edges = sum(len(a) for a in self.adj) // 2
        degrees = [len(a) for a in self.adj]
        
        print(f"\n  {N} matchings, {total_edges} swap edges ({t2-t1:.1f}s)")
        print(f"  Degree: min={min(degrees)}, max={max(degrees)}, "
              f"mean={np.mean(degrees):.1f}, median={np.median(degrees):.0f}")
        
        self.results['N'] = N
        self.results['n_edges'] = total_edges
        self.results['degrees'] = degrees
    
    def probe_5a_degree_distribution(self):
        """5a: Degree distribution by orbit."""
        print(f"\n  ── PROBE 5a: DEGREE DISTRIBUTION ──")
        
        degrees = self.results['degrees']
        
        # Degree histogram
        deg_counts = defaultdict(int)
        for d in degrees:
            deg_counts[d] += 1
        print(f"    Degree histogram:")
        for d in sorted(deg_counts.keys()):
            print(f"      deg={d}: {deg_counts[d]} matchings")
        
        # By orbit
        print(f"\n    Degree by orbit:")
        orbit_deg_data = []
        for oi, (oms, idx) in enumerate(self.orbits):
            degs = [degrees[i] for i in idx]
            tag = ""
            if oi == self.vac_orbit_idx:
                tag = " << VACUUM"
            elif oi == self.fp_orbit_idx:
                tag = " (FP)"
            
            all_same = len(set(degs)) == 1
            deg_str = f"{degs[0]}" if all_same else f"{min(degs)}-{max(degs)} (mean={np.mean(degs):.1f})"
            print(f"    Orbit {oi:>2} (ms={len(idx):>6}): deg = {deg_str}"
                  f"  {'[uniform]' if all_same else ''}{tag}")
            
            orbit_deg_data.append({
                'orbit': oi, 'uniform': all_same,
                'min': min(degs), 'max': max(degs), 'mean': np.mean(degs)
            })
        
        n_uniform = sum(1 for d in orbit_deg_data if d['uniform'])
        print(f"\n    Uniform degree orbits: {n_uniform}/{len(orbit_deg_data)}")
        
        # Is the whole graph regular?
        is_regular = len(deg_counts) == 1
        print(f"    Graph is {'REGULAR' if is_regular else 'NOT regular'}")
        
        self.results['orbit_degrees'] = orbit_deg_data
        self.results['is_regular'] = is_regular
    
    def probe_5b_laplacian_spectrum(self):
        """5b: Laplacian spectrum — spectral gap and mixing time."""
        print(f"\n  ── PROBE 5b: LAPLACIAN SPECTRUM ──")
        
        N = self.results['N']
        
        # For large N, only compute if feasible
        if N > 2000:
            print(f"    [SKIP] N={N} too large for full eigendecomposition")
            return
        
        t0 = time.time()
        L = swap_graph_laplacian(self.adj, N)
        evals = np.sort(np.linalg.eigvalsh(L))
        dt = time.time() - t0
        
        groups = compress_spectrum(evals)
        
        print(f"    Laplacian spectrum ({dt:.1f}s):")
        for v, c in groups:
            if c <= 3 or abs(v) < 1e-6 or abs(v - evals[-1]) < 1e-6:
                print(f"      {v:>10.4f} (x{c})")
            elif groups.index((v, c)) < 8:
                print(f"      {v:>10.4f} (x{c})")
        if len(groups) > 10:
            print(f"      ... ({len(groups)} distinct eigenvalues total)")
            print(f"      {evals[-1]:>10.4f} (x{groups[-1][1]})  [max]")
        
        # Key quantities
        lambda_1 = evals[1] if N > 1 else 0  # spectral gap (Fiedler value)
        lambda_max = evals[-1]
        
        # Connected components = multiplicity of eigenvalue 0
        n_zero = sum(1 for e in evals if abs(e) < 1e-8)
        connected = (n_zero == 1)
        
        print(f"\n    Connected components: {n_zero} "
              f"({'CONNECTED' if connected else 'DISCONNECTED'})")
        print(f"    Spectral gap (λ₁): {lambda_1:.6f}")
        print(f"    Max eigenvalue: {lambda_max:.4f}")
        
        if lambda_1 > 1e-8:
            # Mixing time ~ 1/λ₁ * ln(N)
            mixing_time = np.log(N) / lambda_1
            print(f"    Estimated mixing time: {mixing_time:.2f} steps")
        
        # Algebraic connectivity ratio
        if lambda_max > 0:
            ratio = lambda_1 / lambda_max
            print(f"    λ₁/λ_max = {ratio:.6f} "
                  f"({'expander-like' if ratio > 0.1 else 'not expander'})")
        
        # Also compute adjacency spectrum
        degrees = self.results['degrees']
        A = L * -1
        for i in range(N):
            A[i, i] = 0
        adj_evals = np.sort(np.linalg.eigvalsh(A))[::-1]
        adj_groups = compress_spectrum(adj_evals)
        
        print(f"\n    Adjacency spectrum (top/bottom):")
        for v, c in adj_groups[:5]:
            print(f"      {v:>10.4f} (x{c})")
        if len(adj_groups) > 10:
            print(f"      ...")
        for v, c in adj_groups[-3:]:
            print(f"      {v:>10.4f} (x{c})")
        
        self.results['laplacian_evals'] = evals
        self.results['adj_evals'] = adj_evals
        self.results['spectral_gap'] = lambda_1
        self.results['connected'] = connected
        self.results['n_components'] = n_zero
    
    def probe_5c_sign_grading(self):
        """5c: Does every swap flip the Pfaffian sign? (Bipartite test)"""
        print(f"\n  ── PROBE 5c: PFAFFIAN SIGN × SWAP GRAPH ──")
        
        N = self.results['N']
        n_same = 0
        n_flip = 0
        
        for i in range(N):
            for j in self.adj[i]:
                if j > i:
                    if self.signs[i] == self.signs[j]:
                        n_same += 1
                    else:
                        n_flip += 1
        
        total_edges = n_same + n_flip
        
        print(f"    Same-sign edges: {n_same}")
        print(f"    Sign-flip edges: {n_flip}")
        print(f"    Total edges: {total_edges}")
        
        if n_same == 0:
            print(f"    *** EVERY SWAP FLIPS THE PFAFFIAN SIGN ***")
            print(f"    *** The swap graph is BIPARTITE with Pfaffian sign as the grading ***")
            self.results['sign_bipartite'] = True
        elif n_flip == 0:
            print(f"    *** NO SWAP FLIPS THE SIGN — all edges same-sign ***")
            self.results['sign_bipartite'] = False
        else:
            frac_flip = n_flip / total_edges
            print(f"    Flip fraction: {frac_flip:.6f}")
            self.results['sign_bipartite'] = False
        
        # If bipartite, verify with adjacency spectrum
        if 'adj_evals' in self.results and n_same == 0:
            adj_evals = self.results['adj_evals']
            # Bipartite iff spectrum is symmetric about 0
            n_adj = len(adj_evals)
            symmetric = True
            for i in range(n_adj):
                if abs(adj_evals[i] + adj_evals[n_adj - 1 - i]) > 1e-6:
                    symmetric = False
                    break
            print(f"    Adjacency spectrum symmetric about 0: {symmetric}")
        
        self.results['n_same_sign_edges'] = n_same
        self.results['n_flip_sign_edges'] = n_flip
    
    def probe_5d_orbit_connectivity(self):
        """5d: How many swap edges cross orbit boundaries?"""
        print(f"\n  ── PROBE 5d: ORBIT CONNECTIVITY ──")
        
        N = self.results['N']
        
        # Build matching → orbit map
        m_to_orbit = {}
        m_to_block = {}
        for oi, (oms, idx) in enumerate(self.orbits):
            for ms_key in oms:
                for i in self.blocks[ms_key]:
                    m_to_orbit[i] = oi
                    m_to_block[i] = ms_key
        
        # Count edges: intra-orbit vs inter-orbit
        n_orbits = len(self.orbits)
        orbit_cross = np.zeros((n_orbits, n_orbits), dtype=int)
        
        for i in range(N):
            for j in self.adj[i]:
                if j > i:
                    oi = m_to_orbit[i]
                    oj = m_to_orbit[j]
                    orbit_cross[oi, oj] += 1
                    if oi != oj:
                        orbit_cross[oj, oi] += 1
        
        # Diagonal = intra-orbit, off-diagonal = inter-orbit
        intra = sum(orbit_cross[i, i] for i in range(n_orbits))
        inter = self.results['n_edges'] - intra
        
        print(f"    Intra-orbit edges: {intra}")
        print(f"    Inter-orbit edges: {inter}")
        print(f"    Inter fraction: {inter / self.results['n_edges']:.4f}")
        
        # Which orbit pairs communicate?
        print(f"\n    Orbit communication matrix (nonzero entries):")
        for i in range(n_orbits):
            for j in range(i, n_orbits):
                if orbit_cross[i, j] > 0:
                    tag_i = "V" if i == self.vac_orbit_idx else ("F" if i == self.fp_orbit_idx else "")
                    tag_j = "V" if j == self.vac_orbit_idx else ("F" if j == self.fp_orbit_idx else "")
                    diag = " [self]" if i == j else ""
                    print(f"      ({i}{tag_i},{j}{tag_j}): {orbit_cross[i,j]} edges{diag}")
        
        # Intra-block vs inter-block within orbits
        n_intra_block = 0
        n_inter_block = 0
        for i in range(N):
            for j in self.adj[i]:
                if j > i and m_to_orbit[i] == m_to_orbit[j]:
                    if m_to_block[i] == m_to_block[j]:
                        n_intra_block += 1
                    else:
                        n_inter_block += 1
        
        print(f"\n    Within orbits: {n_intra_block} intra-block, {n_inter_block} inter-block")
        if n_intra_block + n_inter_block > 0:
            print(f"    Intra-block fraction: {n_intra_block / (n_intra_block + n_inter_block):.4f}")
        
        self.results['orbit_cross'] = orbit_cross
        self.results['intra_orbit'] = intra
        self.results['inter_orbit'] = inter
    
    def probe_5e_overlap_correlation(self):
        """5e: Are swap neighbors high or low overlap?"""
        print(f"\n  ── PROBE 5e: SWAP × OVERLAP CORRELATION ──")
        
        N = self.results['N']
        
        # Overlap of swap neighbors
        swap_overlaps = []
        for i in range(N):
            for j in self.adj[i]:
                if j > i:
                    ov = 2 * len(self.match_sets[i] & self.match_sets[j])
                    swap_overlaps.append(ov)
        
        # Random pair overlap for comparison
        rng = np.random.RandomState(42)
        n_sample = min(len(swap_overlaps), 10000)
        random_overlaps = []
        for _ in range(n_sample):
            i, j = rng.randint(0, N, 2)
            while i == j:
                j = rng.randint(0, N)
            random_overlaps.append(2 * len(self.match_sets[i] & self.match_sets[j]))
        
        swap_arr = np.array(swap_overlaps)
        rand_arr = np.array(random_overlaps)
        
        print(f"    Swap neighbors overlap:  mean={np.mean(swap_arr):.4f}, "
              f"std={np.std(swap_arr):.4f}")
        print(f"    Random pairs overlap:    mean={np.mean(rand_arr):.4f}, "
              f"std={np.std(rand_arr):.4f}")
        
        # Overlap distribution for swap neighbors
        ov_counts = defaultdict(int)
        for ov in swap_overlaps:
            ov_counts[ov] += 1
        print(f"\n    Swap-neighbor overlap histogram:")
        for ov in sorted(ov_counts.keys()):
            print(f"      overlap={ov}: {ov_counts[ov]} edges")
        
        # Key structural question: does a swap change exactly 2 edges,
        # so overlap = nv - 4 always? (n edges minus the 2 that changed)
        expected_overlap = self.nv - 4  # n-2 shared edges × 2
        all_expected = all(ov == expected_overlap for ov in swap_overlaps)
        print(f"\n    Expected overlap per swap: {expected_overlap} "
              f"(nv-4 = {self.nv}-4)")
        print(f"    All swap neighbors have overlap {expected_overlap}: {all_expected}")
        
        self.results['swap_overlap_mean'] = np.mean(swap_arr)
        self.results['expected_overlap'] = expected_overlap
        self.results['all_expected_overlap'] = all_expected
    
    def probe_5f_restricted_spectra(self):
        """5f: Laplacian spectrum restricted to each orbit."""
        if self.results['N'] > 2000:
            return
        
        print(f"\n  ── PROBE 5f: ORBIT-RESTRICTED LAPLACIAN ──")
        
        for oi, (oms, idx) in enumerate(self.orbits):
            if len(idx) < 2:
                continue
            
            tag = ""
            if oi == self.vac_orbit_idx:
                tag = " << VACUUM"
            elif oi == self.fp_orbit_idx:
                tag = " (FP)"
            
            bs = len(idx)
            idx_set = set(idx)
            
            # Build induced subgraph
            local_map = {v: i for i, v in enumerate(idx)}
            L_orb = np.zeros((bs, bs))
            
            for i_loc, i_glob in enumerate(idx):
                for j_glob in self.adj[i_glob]:
                    if j_glob in idx_set:
                        j_loc = local_map[j_glob]
                        L_orb[i_loc, j_loc] = -1
                        L_orb[i_loc, i_loc] += 1
            
            evals_orb = np.sort(np.linalg.eigvalsh(L_orb))
            n_components_orb = sum(1 for e in evals_orb if abs(e) < 1e-8)
            gap_orb = evals_orb[n_components_orb] if n_components_orb < bs else 0
            
            # Internal degree
            int_degrees = [int(L_orb[i, i]) for i in range(bs)]
            
            print(f"    Orbit {oi:>2} (bs={bs:>5}){tag}: "
                  f"components={n_components_orb}, gap={gap_orb:.4f}, "
                  f"deg={min(int_degrees)}-{max(int_degrees)}")
    
    def probe_5g_sign_components(self):
        """5g: If swap graph is bipartite, analyze the two sign components."""
        if not self.results.get('sign_bipartite', False):
            return
        
        print(f"\n  ── PROBE 5g: BIPARTITE SIGN STRUCTURE ──")
        
        N = self.results['N']
        plus_idx = [i for i in range(N) if self.signs[i] > 0]
        minus_idx = [i for i in range(N) if self.signs[i] < 0]
        
        # Build the "double-step" graph: neighbors at distance 2 in swap graph
        # (same-sign matchings connected by two swaps)
        print(f"    (+) partition: {len(plus_idx)} matchings")
        print(f"    (-) partition: {len(minus_idx)} matchings")
        
        # For each partition, compute the distance-2 graph (bipartite square)
        for label, part in [("(+)", plus_idx), ("(-)", minus_idx)]:
            part_set = set(part)
            local_map = {v: i for i, v in enumerate(part)}
            bs = len(part)
            
            # Distance-2 neighbors: j reachable from i via i->k->j where k is opposite sign
            d2_adj = [set() for _ in range(bs)]
            for i_loc, i_glob in enumerate(part):
                for k in self.adj[i_glob]:  # k has opposite sign
                    for j_glob in self.adj[k]:
                        if j_glob in part_set and j_glob != i_glob:
                            d2_adj[i_loc].add(local_map[j_glob])
            
            d2_edges = sum(len(a) for a in d2_adj) // 2
            d2_degrees = [len(a) for a in d2_adj]
            
            print(f"\n    {label} double-step graph:")
            print(f"      Edges: {d2_edges}")
            print(f"      Degree: min={min(d2_degrees)}, max={max(d2_degrees)}, "
                  f"mean={np.mean(d2_degrees):.1f}")
            
            # Is it connected?
            visited = set()
            queue = [0]
            visited.add(0)
            while queue:
                cur = queue.pop(0)
                for nb in d2_adj[cur]:
                    if nb not in visited:
                        visited.add(nb)
                        queue.append(nb)
            connected = len(visited) == bs
            print(f"      Connected: {connected}")
    
    def probe_5h_direction_changes(self):
        """5h: How does a swap change the direction multiset?"""
        print(f"\n  ── PROBE 5h: DIRECTION CHANGE UNDER SWAP ──")
        
        N = self.results['N']
        
        # For each swap edge, compute direction change
        dir_changes = defaultdict(int)
        n_same_dir = 0
        n_diff_dir = 0
        
        for i in range(N):
            ms_i = torus_ms(self.raw[i], self.p)
            for j in self.adj[i]:
                if j > i:
                    ms_j = torus_ms(self.raw[j], self.p)
                    if ms_i == ms_j:
                        n_same_dir += 1
                    else:
                        n_diff_dir += 1
                        # What changed?
                        diff = tuple(sorted(set(ms_j) - set(ms_i)))
                        dir_changes[diff] += 1
        
        print(f"    Same direction multiset: {n_same_dir} edges")
        print(f"    Different direction multiset: {n_diff_dir} edges")
        print(f"    Direction-preserving fraction: "
              f"{n_same_dir / (n_same_dir + n_diff_dir):.4f}")
        
        # How many swaps stay within the same orbit?
        self.results['n_same_dir'] = n_same_dir
        self.results['n_diff_dir'] = n_diff_dir
    
    def run_all(self):
        print(f"\n{'=' * 72}")
        print(f"  SWAP GRAPH PROBE: K_{self.nv} (p={self.p}, n={self.n})")
        print(f"{'=' * 72}")
        
        self.generate()
        self.probe_5a_degree_distribution()
        self.probe_5b_laplacian_spectrum()
        self.probe_5c_sign_grading()
        self.probe_5d_orbit_connectivity()
        self.probe_5e_overlap_correlation()
        self.probe_5f_restricted_spectra()
        self.probe_5g_sign_components()
        self.probe_5h_direction_changes()


# ======================================================================
# CROSS-LEVEL SUMMARY
# ======================================================================

def cross_level_summary(probes):
    print(f"\n{'=' * 72}")
    print(f"  SWAP GRAPH PROBE: CROSS-LEVEL SUMMARY")
    print(f"{'=' * 72}")
    
    print(f"\n  GRAPH STRUCTURE:")
    print(f"  {'Level':<8} {'N':>7} {'Edges':>8} {'Regular?':>9} {'Degree':>8} "
          f"{'Connected':>10}")
    print(f"  {'-' * 55}")
    for pr in probes:
        r = pr.results
        deg_str = f"{r['degrees'][0]}" if r.get('is_regular') else \
                  f"{min(r['degrees'])}-{max(r['degrees'])}"
        conn = 'YES' if r.get('connected', '?') else f"NO ({r.get('n_components', '?')})"
        print(f"  K_{pr.nv:<5} {r['N']:>7} {r['n_edges']:>8} "
              f"{'YES' if r.get('is_regular') else 'NO':>9} {deg_str:>8} {conn:>10}")
    
    print(f"\n  PFAFFIAN BIPARTITENESS:")
    print(f"  {'Level':<8} {'Same-sign':>10} {'Flip':>10} {'Bipartite?':>11}")
    print(f"  {'-' * 45}")
    for pr in probes:
        r = pr.results
        print(f"  K_{pr.nv:<5} {r.get('n_same_sign_edges', '?'):>10} "
              f"{r.get('n_flip_sign_edges', '?'):>10} "
              f"{'YES' if r.get('sign_bipartite') else 'NO':>11}")
    
    print(f"\n  SPECTRAL GAP:")
    print(f"  {'Level':<8} {'λ₁':>10} {'λ_max':>10} {'λ₁/λ_max':>10}")
    print(f"  {'-' * 42}")
    for pr in probes:
        r = pr.results
        gap = r.get('spectral_gap')
        if gap is not None:
            lmax = r['laplacian_evals'][-1]
            print(f"  K_{pr.nv:<5} {gap:>10.4f} {lmax:>10.4f} {gap/lmax:>10.6f}")
    
    print(f"\n  ORBIT PERMEABILITY:")
    print(f"  {'Level':<8} {'Intra':>8} {'Inter':>8} {'Inter%':>8} "
          f"{'Dir-same':>9} {'Dir-diff':>9}")
    print(f"  {'-' * 55}")
    for pr in probes:
        r = pr.results
        total = r.get('intra_orbit', 0) + r.get('inter_orbit', 0)
        inter_pct = r.get('inter_orbit', 0) / total if total > 0 else 0
        print(f"  K_{pr.nv:<5} {r.get('intra_orbit', '?'):>8} "
              f"{r.get('inter_orbit', '?'):>8} {inter_pct:>7.2%} "
              f"{r.get('n_same_dir', '?'):>9} {r.get('n_diff_dir', '?'):>9}")
    
    print(f"\n  SWAP OVERLAP:")
    print(f"  {'Level':<8} {'Expected':>9} {'Actual':>9} {'Exact?':>7}")
    print(f"  {'-' * 38}")
    for pr in probes:
        r = pr.results
        print(f"  K_{pr.nv:<5} {r.get('expected_overlap', '?'):>9} "
              f"{r.get('swap_overlap_mean', 0):>9.4f} "
              f"{'YES' if r.get('all_expected_overlap') else 'NO':>7}")


# ======================================================================
# MAIN
# ======================================================================

def main():
    parser = argparse.ArgumentParser(description='Swap Graph Probe (Sense 5)')
    parser.add_argument('--levels', nargs='+', type=int, default=[6, 8, 10, 12],
                        help='K_{2n} levels to analyze (default: 6 8 10 12)')
    args = parser.parse_args()
    
    t_total = time.time()
    print("SWAP GRAPH PROBE (SENSE 5)")
    print(f"Levels: {['K_' + str(lv) for lv in args.levels]}")
    print("=" * 72)
    
    probes = []
    for nv in args.levels:
        pr = SwapGraphProbe(nv)
        pr.run_all()
        probes.append(pr)
    
    cross_level_summary(probes)
    
    print(f"\n{'=' * 72}")
    print(f"  TOTAL TIME: {time.time() - t_total:.1f}s")
    print(f"{'=' * 72}")


if __name__ == '__main__':
    main()
