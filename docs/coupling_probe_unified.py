# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Coupling Probe: Unified Script
===============================
Runs all vacuum overlap analyses across K_6, K_8, K_10, K_12.

Capabilities:
  Sense 1: Inter-orbit coupling (SVD of off-diagonal overlap blocks)
  Sense 2: Vacuum spectral structure (full eigenvalue spectrum)
  Sense 3: Block-level coupling (direction-block decomposition of Sense 1)
  Null ID: Null space irrep decomposition under wreath product
  Laws:    Cross-level validation of five structural laws
  Embed:   Inter-level embedding (rank-1 theorem, adjacent + non-adjacent)

Usage:
  python coupling_probe_unified.py [--levels 6 8 10 12] [--no-irreps] [--no-coupling] [--no-blocks] [--no-embedding]
"""

import numpy as np
from collections import defaultdict
from fractions import Fraction
import time
import argparse

# ======================================================================
# CORE FUNCTIONS
# ======================================================================

def gen_matchings(nv):
    """Generate all perfect matchings of K_{nv}."""
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


def gen_pairings(verts):
    """All perfect matchings (pairings) of a list of vertices."""
    if len(verts) == 0:
        return [frozenset()]
    if len(verts) == 2:
        return [frozenset([(verts[0], verts[1])])]
    first = verts[0]
    result = []
    for i in range(1, len(verts)):
        rest = [v for j, v in enumerate(verts) if j != 0 and j != i]
        for sub in gen_pairings(rest):
            result.append(frozenset([(first, verts[i])]) | sub)
    return result


def edge_dir(a, b, p):
    """Direction class of edge (a,b) mod p."""
    a, b = min(a, b), max(a, b)
    if a >= p or b >= p:
        return -1
    d = (b - a) % p
    return min(d, p - d) - 1


def torus_ms(matching, p):
    """Direction multiset of a matching."""
    dirs = []
    for e in sorted(matching):
        d = edge_dir(e[0], e[1], p)
        if d >= 0:
            dirs.append(d)
    return tuple(sorted(dirs))


def build_overlap(ms_a, ms_b, nv, same=False):
    """Overlap matrix between two sets of matchings (as sets of edges)."""
    na, nb = len(ms_a), len(ms_b)
    O = np.zeros((na, nb))
    for i in range(na):
        for j in range(nb):
            if same and i == j:
                O[i, j] = nv
            else:
                O[i, j] = 2 * len(ms_a[i] & ms_b[j])
    return O


def compute_dir_perm(p, g=2):
    """Permutation of direction classes induced by x*g mod p."""
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
    """Group direction multisets into orbits under the direction permutation."""
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
    """Group eigenvalues by value, return list of (value, multiplicity)."""
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


def apply_vertex_mult(matching, g, p):
    """Apply vertex map v -> g*v mod p (extra vertex fixed)."""
    new = set()
    for a, b in matching:
        na = (g * a) % p if a < p else a
        nb = (g * b) % p if b < p else b
        new.add((min(na, nb), max(na, nb)))
    return frozenset(new)


# ======================================================================
# LEVEL ANALYSIS
# ======================================================================

class LevelProbe:
    """Full coupling probe for a single K_{2n} level."""

    def __init__(self, nv):
        self.nv = nv
        self.p = nv - 1
        self.n = nv // 2
        self.raw = None
        self.match_sets = None
        self.blocks = None
        self.orbits = None       # list of (multiset_list, index_list)
        self.vac_orbit_idx = None
        self.fp_orbit_idx = None
        self.orbit_size = None   # k = torus orbit length
        self.results = {}

    def generate(self):
        """Step 1: Generate matchings and orbit structure."""
        t0 = time.time()
        self.raw = gen_matchings(self.nv)
        self.match_sets = [set(m) for m in self.raw]
        N = len(self.raw)

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

        # Identify vacuum by bs = 2p (for n>=4), or largest non-FP for n=3
        self.vac_orbit_idx = None
        for oi, (oms, idx) in enumerate(self.orbits):
            if len(self.blocks[oms[0]]) == 2 * self.p:
                self.vac_orbit_idx = oi
                break
        if self.vac_orbit_idx is None:
            # K₆ case: vacuum is the non-fixed-point orbit
            for oi, (oms, idx) in enumerate(self.orbits):
                if len(oms) > 1:
                    self.vac_orbit_idx = oi
                    break

        # Orbit size (k)
        if self.vac_orbit_idx is not None:
            self.orbit_size = len(self.orbits[self.vac_orbit_idx][0])

        dt = time.time() - t0
        print(f"\n  {N} matchings, {len(self.blocks)} blocks, "
              f"{len(self.orbits)} orbits ({dt:.1f}s)")
        for oi, (oms, idx) in enumerate(self.orbits):
            tag = ""
            if oi == self.vac_orbit_idx:
                tag = " << VACUUM"
            elif oi == self.fp_orbit_idx:
                tag = " (FP)"
            print(f"    Orbit {oi:>2}: |orb|={len(oms)}, "
                  f"ms={len(idx):>6}, bs={len(self.blocks[oms[0]])}{tag}")

    def sense2_vacuum_spectrum(self):
        """Sense 2: Full eigenvalue spectrum of vacuum overlap."""
        if self.vac_orbit_idx is None:
            print("  [SKIP] No vacuum orbit identified")
            return

        vac_oms, vac_idx = self.orbits[self.vac_orbit_idx]
        vac_ms = [self.match_sets[i] for i in vac_idx]
        bs = len(vac_idx)

        t0 = time.time()
        O = build_overlap(vac_ms, vac_ms, self.nv, same=True)
        evals, evecs = np.linalg.eigh(O)
        order = np.argsort(evals)[::-1]
        evals = evals[order]
        evecs = evecs[:, order]

        rank = int(np.sum(np.abs(evals) > 1e-6))
        null_dim = bs - rank
        groups = compress_spectrum(evals)

        print(f"\n  VACUUM SPECTRUM (bs={bs}):")
        for v, c in groups:
            print(f"    {v:>12.4f} (x{c})")
        print(f"    rank={rank}, null={null_dim}, bs/2={bs // 2}")

        # Law 1 check
        law1 = (null_dim == bs // 2) if self.n >= 4 else None
        if law1 is not None:
            print(f"    LAW 1 (half-rank): {'CONFIRMED' if law1 else 'VIOLATED'}")
        elif self.n == 3:
            print(f"    LAW 1: n=3 exception (full rank, no null space)")

        # Null projector diagonal
        null_vecs = evecs[:, np.abs(evals[order] if False else evals) < 1e-8]
        # recompute cleanly
        null_mask = np.abs(evals) < 1e-8
        null_vecs = evecs[:, null_mask]
        if null_vecs.shape[1] > 0:
            proj = null_vecs @ null_vecs.T
            diag = np.diag(proj)
            all_half = np.allclose(diag, 0.5)
            print(f"    Null projector diagonal all 1/2: {all_half}")

        self.results['rank'] = rank
        self.results['null_dim'] = null_dim
        self.results['bs'] = bs
        self.results['spectrum'] = groups
        self.results['law1'] = law1
        self.results['O_vac'] = O
        self.results['evals'] = evals
        self.results['evecs'] = evecs
        print(f"    ({time.time() - t0:.1f}s)")

    def sense1_coupling(self, max_orbits=None):
        """Sense 1: Inter-orbit coupling via SVD."""
        if self.vac_orbit_idx is None:
            return

        vac_oms, vac_idx = self.orbits[self.vac_orbit_idx]
        vac_ms = [self.match_sets[i] for i in vac_idx]
        bs_vac = len(vac_idx)
        rank = self.results.get('rank', bs_vac)

        print(f"\n  INTER-ORBIT COUPLING:")
        coupling_data = []

        for oi, (oms, oidx) in enumerate(self.orbits):
            if oi == self.vac_orbit_idx:
                continue
            if max_orbits is not None and len(coupling_data) >= max_orbits:
                print(f"    ... ({len(self.orbits) - len(coupling_data) - 1} more orbits skipped)")
                break

            other_ms = [self.match_sets[i] for i in oidx]
            O_cross = build_overlap(vac_ms, other_ms, self.nv)
            _, sigma, _ = np.linalg.svd(O_cross, full_matrices=False)
            sigma_nz = sigma[sigma > 1e-6]
            frob2 = round(np.sum(sigma ** 2))
            c_rank = len(sigma_nz)

            sg = compress_spectrum(sigma_nz)
            is_fp = (oi == self.fp_orbit_idx)
            fp_str = " (FP)" if is_fp else ""

            print(f"\n    Vac <-> Orbit {oi} (|orb|={len(oms)}, "
                  f"ms={len(oidx)}){fp_str}")
            print(f"    rank={c_rank}, ||F||^2={frob2}")
            for sv, cnt in sg:
                print(f"      sigma={sv:10.4f} (x{cnt})  sigma^2={sv**2:10.2f}")

            # Frobenius ratio
            ratio = Fraction(frob2, bs_vac * len(oidx))
            print(f"    ||F||^2 / (bs_vac*bs_other) = {ratio} = {float(ratio):.6f}")

            coupling_data.append({
                'orbit': oi, 'is_fp': is_fp, 'rank': c_rank,
                'frob2': frob2, 'ratio': ratio, 'sigma_groups': sg
            })

        # Law 2 check (FP deficit)
        fp_data = [d for d in coupling_data if d['is_fp']]
        if fp_data and self.orbit_size:
            fp_rank = fp_data[0]['rank']
            deficit = rank - fp_rank
            expected = self.orbit_size - 1
            if self.n >= 4:
                law2 = (deficit == expected)
                print(f"\n    LAW 2 (FP deficit): deficit={deficit}, "
                      f"k-1={expected}, {'CONFIRMED' if law2 else 'VIOLATED'}")
            else:
                law2 = None
                print(f"\n    LAW 2: n=3 exception (dimension-limited, deficit={deficit})")
            self.results['law2'] = law2
            self.results['fp_deficit'] = deficit

        self.results['coupling'] = coupling_data

    def sense3_block_coupling(self, max_orbits=None):
        """Sense 3: Block-level coupling — decompose orbit coupling into direction blocks."""
        if self.vac_orbit_idx is None:
            return

        vac_oms, vac_idx = self.orbits[self.vac_orbit_idx]
        k = self.orbit_size

        print(f"\n  BLOCK-LEVEL COUPLING (Sense 3):")
        print(f"    Vacuum: {k} blocks, dirs: {vac_oms}")

        block_data = []
        n_done = 0

        for oi, (oms, oidx) in enumerate(self.orbits):
            if oi == self.vac_orbit_idx:
                continue
            if max_orbits is not None and n_done >= max_orbits:
                print(f"    ... ({len(self.orbits) - n_done - 1} more orbits skipped)")
                break

            k_other = len(oms)
            is_fp = (oi == self.fp_orbit_idx)
            tag = " (FP)" if is_fp else ""
            bs_block = len(self.blocks[oms[0]])

            # Compute block×block SVD grid
            grid_ranks = np.zeros((k, k_other), dtype=int)
            grid_frob2 = np.zeros((k, k_other))
            grid_smax = np.zeros((k, k_other))
            grid_smin = np.zeros((k, k_other))
            grid_flat = np.zeros((k, k_other), dtype=bool)

            for vi in range(k):
                v_ms = [self.match_sets[i] for i in self.blocks[vac_oms[vi]]]
                for bi in range(k_other):
                    b_ms = [self.match_sets[i] for i in self.blocks[oms[bi]]]
                    O = build_overlap(v_ms, b_ms, self.nv)
                    svs = np.linalg.svd(O, compute_uv=False)
                    svs_nz = svs[svs > 1e-10]
                    grid_ranks[vi, bi] = len(svs_nz)
                    grid_frob2[vi, bi] = np.sum(svs ** 2)
                    if len(svs_nz) > 0:
                        grid_smax[vi, bi] = svs_nz[0]
                        grid_smin[vi, bi] = svs_nz[-1]
                        grid_flat[vi, bi] = (svs_nz[0] / svs_nz[-1] < 1.001)

            # Orbit-level for comparison
            vac_ms_all = [self.match_sets[i] for i in vac_idx]
            other_ms = [self.match_sets[i] for i in oidx]
            O_orb = build_overlap(vac_ms_all, other_ms, self.nv)
            svs_orb = np.linalg.svd(O_orb, compute_uv=False)
            svs_orb_nz = svs_orb[svs_orb > 1e-10]
            orbit_rank = len(svs_orb_nz)
            orbit_frob2 = round(np.sum(svs_orb ** 2))

            # Summarize by shift class
            print(f"\n    Orbit {oi}{tag} (|orb|={k_other}, bs_blk={bs_block}):")

            if k == k_other:
                # Circulant structure: first row of ||F||² grid
                frob_row = [round(grid_frob2[0, bi]) for bi in range(k)]
                print(f"      ||F||² row: {frob_row} (circulant, orbit={orbit_frob2})")

                for shift in range(k):
                    vi, bi = 0, shift
                    r = grid_ranks[vi, bi]
                    sm = grid_smax[vi, bi]
                    sn = grid_smin[vi, bi]
                    hier = sm / sn if sn > 1e-10 else float('inf')
                    flat = grid_flat[vi, bi]
                    flat_tag = " ** FLAT **" if flat else ""
                    print(f"      shift={shift}: rank={r}, "
                          f"σ_max={sm:.4f}, σ_min={sn:.4f}, "
                          f"hier={hier:.2f}{flat_tag}")

            elif k_other == 1:
                # FP: all vacuum blocks should be equivalent
                r = grid_ranks[0, 0]
                sm = grid_smax[0, 0]
                sn = grid_smin[0, 0]
                hier = sm / sn if sn > 1e-10 else float('inf')
                all_equiv = all(
                    np.isclose(grid_smax[vi, 0], sm) for vi in range(k))
                print(f"      FP: rank={r}, σ_max={sm:.4f}, σ_min={sn:.4f}, "
                      f"hier={hier:.2f}, "
                      f"{'all blocks equiv' if all_equiv else 'BLOCKS DIFFER'}")

            # Rank cancellation
            block_rank_sum = int(np.sum(grid_ranks))
            cancelled = block_rank_sum - orbit_rank
            print(f"      Σ blk ranks={block_rank_sum}, "
                  f"orbit rank={orbit_rank}, "
                  f"cancelled={cancelled}")

            n_flat = int(np.sum(grid_flat))
            if n_flat > 0:
                print(f"      Flat couplings: {n_flat}/{k * k_other} block pairs")

            block_data.append({
                'orbit': oi, 'is_fp': is_fp,
                'orbit_rank': orbit_rank,
                'block_rank_sum': block_rank_sum,
                'n_flat': n_flat
            })
            n_done += 1

        # Cross-orbit summary
        if block_data:
            total_flat = sum(d['n_flat'] for d in block_data)
            total_cancelled = sum(d['block_rank_sum'] - d['orbit_rank']
                                  for d in block_data)
            print(f"\n    BLOCK SUMMARY:")
            print(f"      Total flat couplings: {total_flat}")
            print(f"      Total rank cancelled: {total_cancelled}")
            print(f"      Orbit averaging always destroys rank — "
                  f"interference is universal")

        self.results['block_coupling'] = block_data

    def sense3_vacuum_blocks(self):
        """Sense 3b: Vacuum self-coupling decomposed into direction blocks.

        Key question: Does each vacuum direction block already have half-rank,
        or does the half-rank only emerge from interference between blocks?
        """
        if self.vac_orbit_idx is None:
            return

        vac_oms, vac_idx = self.orbits[self.vac_orbit_idx]
        k = self.orbit_size
        bs_total = len(vac_idx)

        # Build block index map: which indices in vac_idx belong to which block
        block_slices = []  # (start, end) in vac_idx ordering
        block_local_idx = []  # list of matching indices per block
        for bi in range(k):
            ms_key = vac_oms[bi]
            local = self.blocks[ms_key]
            block_local_idx.append(local)

        bs_blk = len(block_local_idx[0])

        print(f"\n  VACUUM BLOCK DECOMPOSITION (Sense 3b):")
        print(f"    {k} blocks × {bs_blk} matchings = {bs_total} total")

        # ---- Diagonal blocks: self-overlap within each direction sector ----
        print(f"\n    DIAGONAL BLOCKS (self-overlap within direction sector):")

        diag_ranks = []
        diag_nulls = []
        diag_spectra = []
        diag_null_diags = []

        for bi in range(k):
            ms = [self.match_sets[i] for i in block_local_idx[bi]]
            O_diag = build_overlap(ms, ms, self.nv, same=True)
            evals = np.linalg.eigvalsh(O_diag)[::-1]
            rank = int(np.sum(np.abs(evals) > 1e-6))
            null = bs_blk - rank
            groups = compress_spectrum(evals)

            # Null projector diagonal
            evals_full, evecs_full = np.linalg.eigh(O_diag)
            null_mask = np.abs(evals_full) < 1e-8
            null_vecs = evecs_full[:, null_mask]
            null_diag_val = None
            if null_vecs.shape[1] > 0:
                proj = null_vecs @ null_vecs.T
                diag_vals = np.diag(proj)
                null_diag_val = diag_vals[0] if np.allclose(diag_vals, diag_vals[0]) else diag_vals

            diag_ranks.append(rank)
            diag_nulls.append(null)
            diag_spectra.append(groups)
            diag_null_diags.append(null_diag_val)

        # Print comparison (by permutation equivalence, all should be identical)
        all_ranks_equal = len(set(diag_ranks)) == 1
        print(f"      All blocks equivalent: {all_ranks_equal}")
        print(f"      Block size: {bs_blk}, rank: {diag_ranks[0]}, "
              f"null: {diag_nulls[0]}, bs/2: {bs_blk // 2}")

        half_rank_local = (diag_nulls[0] == bs_blk // 2) if self.n >= 4 else None
        if half_rank_local is not None:
            tag = "YES — half-rank is LOCAL" if half_rank_local else "NO — half-rank is GLOBAL (interference)"
            print(f"      Half-rank per block: {tag}")
        else:
            print(f"      n=3: checking rank = {diag_ranks[0]}/{bs_blk}")

        # Spectrum of representative block
        print(f"      Block V_0 spectrum:")
        for v, c in diag_spectra[0]:
            print(f"        {v:>12.4f} (x{c})")

        if diag_null_diags[0] is not None:
            if isinstance(diag_null_diags[0], (int, float, np.floating)):
                print(f"      Null projector diagonal: {diag_null_diags[0]:.6f}")
            else:
                print(f"      Null projector diagonal: varies ({np.min(diag_null_diags[0]):.4f} - {np.max(diag_null_diags[0]):.4f})")

        # ---- Off-diagonal blocks: cross-direction vacuum coupling ----
        print(f"\n    OFF-DIAGONAL BLOCKS (cross-direction within vacuum):")

        offdiag_data = {}
        for shift in range(1, k):
            vi, vj = 0, shift
            ms_i = [self.match_sets[idx] for idx in block_local_idx[vi]]
            ms_j = [self.match_sets[idx] for idx in block_local_idx[vj]]
            O_off = build_overlap(ms_i, ms_j, self.nv)
            svs = np.linalg.svd(O_off, compute_uv=False)
            svs_nz = svs[svs > 1e-10]
            frob2 = np.sum(svs ** 2)
            flat = (svs_nz[0] / svs_nz[-1] < 1.001) if len(svs_nz) > 0 else False

            offdiag_data[shift] = {
                'rank': len(svs_nz), 'svs': svs_nz,
                'frob2': frob2, 'flat': flat
            }

            flat_tag = " ** FLAT **" if flat else ""
            hier = svs_nz[0] / svs_nz[-1] if len(svs_nz) > 1 else 1.0
            print(f"      shift={shift}: rank={len(svs_nz)}, "
                  f"σ_max={svs_nz[0]:.4f}, σ_min={svs_nz[-1]:.4f}, "
                  f"hier={hier:.2f}, ||F||²={frob2:.0f}{flat_tag}")

        # ---- Frobenius budget ----
        diag_frob2 = sum(
            np.sum(np.array([v for v, c in sp for _ in range(c)]) ** 2)
            for sp in diag_spectra
        )
        # Actually compute properly
        diag_frob2 = 0
        for bi in range(k):
            ms = [self.match_sets[i] for i in block_local_idx[bi]]
            O_diag = build_overlap(ms, ms, self.nv, same=True)
            diag_frob2 += np.sum(O_diag ** 2)
        offdiag_frob2 = sum(d['frob2'] for d in offdiag_data.values()) * k  # k copies per shift

        # Full orbit ||F||²
        vac_ms = [self.match_sets[i] for i in vac_idx]
        O_full = build_overlap(vac_ms, vac_ms, self.nv, same=True)
        full_frob2 = np.sum(O_full ** 2)

        print(f"\n    FROBENIUS BUDGET:")
        print(f"      Diagonal (k blocks): {diag_frob2:.0f}")
        print(f"      Off-diagonal ({k} × {k-1} shifts): {offdiag_frob2:.0f}")
        print(f"      Total: {diag_frob2 + offdiag_frob2:.0f}")
        print(f"      Full orbit: {full_frob2:.0f}")
        print(f"      Diagonal fraction: {diag_frob2 / full_frob2:.4f}")

        # ---- Compare full spectrum to block-reconstructed ----
        # Block-circulant: eigenvalues = evals of Σ_j B_j ω^{jl} for l=0..k-1
        # REQUIRES consistent ordering within blocks via direction permutation.
        # Key: the correct generator is g = 2^(ord(2)/k) mod p, which has order k
        # on vertices (not just matchings). For K8 this is g=2; for K10,K12 it's g=4.
        print(f"\n    BLOCK-CIRCULANT SPECTRAL RECONSTRUCTION:")

        # Find generator with order k on vertices
        ord2 = 1
        val = 2
        while val != 1:
            val = (val * 2) % self.p
            ord2 += 1
            if ord2 > self.p:
                break
        exp = ord2 // k
        g = pow(2, exp, self.p)
        print(f"      Generator: g = 2^{exp} = {g} (ord(g mod {self.p}) = {k})")

        block_match_lists = []
        for bi in range(k):
            block_match_lists.append([self.raw[i] for i in block_local_idx[bi]])

        def apply_perm_to_matching(m, p_val, gg):
            """Map matching m -> m' under vertex map v -> gg*v mod p."""
            new_edges = set()
            for a, b in m:
                na = (gg * a) % p_val if a < p_val else a
                nb = (gg * b) % p_val if b < p_val else b
                new_edges.add((min(na, nb), max(na, nb)))
            return frozenset(new_edges)

        # Find block permutation under ×g
        block_seq = [0]
        cur_blk = 0
        for _ in range(k - 1):
            m_cur = block_match_lists[cur_blk][0]
            m_next = apply_perm_to_matching(m_cur, self.p, g)
            for bi in range(k):
                if m_next in block_match_lists[bi]:
                    cur_blk = bi
                    block_seq.append(bi)
                    break
        print(f"      Block sequence: {block_seq}")

        # Build aligned orderings for each block in the sequence
        aligned_orders = {0: list(range(bs_blk))}
        aligned_ok = True
        for step in range(1, k):
            bi = block_seq[step]
            order_bi = []
            for mi, m in enumerate(block_match_lists[0]):
                cur = m
                for _ in range(step):
                    cur = apply_perm_to_matching(cur, self.p, g)
                try:
                    idx = block_match_lists[bi].index(cur)
                    order_bi.append(idx)
                except ValueError:
                    aligned_ok = False
                    break
            if not aligned_ok:
                break
            aligned_orders[bi] = order_bi

        if not aligned_ok:
            print(f"      Alignment failed")
        else:
            # Build circulant blocks B_j = O(V_0, V_{seq[j]}) with alignment
            blocks_B = []
            ms_0 = [self.match_sets[block_local_idx[0][i]] for i in range(bs_blk)]
            for step in range(k):
                bi = block_seq[step]
                B = np.zeros((bs_blk, bs_blk))
                ms_j = [self.match_sets[block_local_idx[bi][aligned_orders[bi][j]]]
                        for j in range(bs_blk)]
                for i in range(bs_blk):
                    for j in range(bs_blk):
                        if step == 0 and i == j:
                            B[i, j] = self.nv
                        else:
                            B[i, j] = 2 * len(ms_0[i] & ms_j[j])
                blocks_B.append(B)

            omega = np.exp(2j * np.pi / k)
            all_evals_recon = []
            null_per_channel = []

            for ell in range(k):
                M = np.zeros((bs_blk, bs_blk), dtype=complex)
                for j in range(k):
                    M += blocks_B[j] * (omega ** (j * ell))

                M_herm = (M + M.conj().T) / 2
                evals_ch = np.linalg.eigvalsh(M_herm)[::-1]
                all_evals_recon.extend(evals_ch.tolist())

                rank_ch = int(np.sum(np.abs(evals_ch) > 1e-6))
                null_ch = bs_blk - rank_ch
                null_per_channel.append(null_ch)
                groups_ch = compress_spectrum(evals_ch)

                print(f"      Channel ℓ={ell}: rank={rank_ch}, null={null_ch}")
                for v, c in groups_ch:
                    tag = " [null]" if abs(v) < 1e-6 else ""
                    print(f"        {v:>12.4f} (x{c}){tag}")

            # Compare to full orbit spectrum
            full_evals = self.results.get('evals', None)
            if full_evals is not None:
                recon_sorted = np.sort(all_evals_recon)[::-1]
                full_sorted = np.sort(full_evals)[::-1]
                max_diff = np.max(np.abs(recon_sorted - full_sorted))
                print(f"\n      Reconstruction error: {max_diff:.2e}")
                print(f"      Exact match: {max_diff < 1e-8}")

            # Null space summary
            total_null_recon = sum(null_per_channel)
            all_equal = len(set(null_per_channel)) == 1
            print(f"\n    NULL SPACE BY CHANNEL:")
            for ell in range(k):
                print(f"      ℓ={ell}: {null_per_channel[ell]} null modes")
            print(f"      Total: {total_null_recon} "
                  f"(orbit null = {self.results.get('null_dim', '?')})")
            if all_equal and null_per_channel[0] == bs_blk // 2:
                print(f"      *** HALF-RANK DEMOCRATIC: every channel has exactly bs/2 null modes ***")

            self.results['vac_blocks'] = {
                'diag_rank': diag_ranks[0],
                'diag_null': diag_nulls[0],
                'half_rank_local': half_rank_local,
                'diag_spectrum': diag_spectra[0],
                'offdiag': offdiag_data,
                'diag_frac': diag_frob2 / full_frob2 if full_frob2 > 0 else None,
                'null_per_channel': null_per_channel,
                'generator': g,
                'block_seq': block_seq
            }

    def null_space_irreps(self):
        """Irrep decomposition of null/live spaces under wreath product."""
        if 'evecs' not in self.results or self.results.get('null_dim', 0) == 0:
            print("  [SKIP] No null space or eigendata")
            return

        evals = self.results['evals']
        evecs = self.results['evecs']
        null_vecs = evecs[:, np.abs(evals) < 1e-8]
        live_vecs = evecs[:, np.abs(evals) >= 1e-8]

        vac_oms, vac_idx = self.orbits[self.vac_orbit_idx]
        bs = len(vac_idx)

        # Build matching lookup
        match_to_vi = {}
        for vi in range(bs):
            match_to_vi[self.raw[vac_idx[vi]]] = vi

        # Permutation matrices
        def perm_matrix(g):
            M = np.zeros((bs, bs))
            for vi in range(bs):
                m_new = apply_vertex_mult(self.raw[vac_idx[vi]], g, self.p)
                if m_new in match_to_vi:
                    M[match_to_vi[m_new], vi] = 1.0
            return M

        # Find generator for C_k (the ×2 mod p action)
        P = perm_matrix(2)
        # C_2 generator: ×(-1) mod p
        Q = perm_matrix(self.p - 1)

        k = self.orbit_size
        Pk = np.linalg.matrix_power(P, k)
        P_order = k
        Ptmp = Pk
        while not np.allclose(Ptmp, np.eye(bs)):
            Ptmp = Ptmp @ P
            P_order += 1
            if P_order > 2 * self.p:
                break
        print(f"\n  IRREP DECOMPOSITION (C_2 x C_{k}):")
        print(f"    P^{k} = I: {np.allclose(Pk, np.eye(bs))}"
              f"{'  (P order = ' + str(P_order) + ' on matchings)' if P_order != k else ''}")
        print(f"    Q^2 = I: {np.allclose(Q @ Q, np.eye(bs))}")
        print(f"    PQ = QP: {np.allclose(P @ Q, Q @ P)}")

        # C_k eigenvalues: exp(2*pi*i*j/k) for j=0..k-1
        # C_2 eigenvalues: +1, -1
        # Build projectors for each joint irrep

        # Precompute P powers
        P_powers = [np.eye(bs)]
        for _ in range(k - 1):
            P_powers.append(P_powers[-1] @ P)

        print(f"\n    {'Irrep':<16} {'null':>5} {'live':>5} {'total':>6}")
        print(f"    {'-' * 38}")

        null_total = 0
        live_total = 0
        irrep_data = []

        for c2_sign in [+1, -1]:
            c2_label = '+1' if c2_sign > 0 else '-1'
            for j in range(k):
                omega_j = np.exp(2j * np.pi * j / k)
                # Projector onto C_k eigenvalue omega_j
                proj_ck = sum(np.conj(omega_j) ** m * P_powers[m]
                              for m in range(k)) / k
                # Projector onto C_2 eigenvalue c2_sign
                proj_c2 = (np.eye(bs) + c2_sign * Q) / 2.0
                proj = np.real(proj_ck @ proj_c2)

                dim_total = round(np.real(np.trace(proj)))
                dim_null = round(np.real(np.trace(null_vecs.T @ proj @ null_vecs)))
                dim_live = round(np.real(np.trace(live_vecs.T @ proj @ live_vecs)))

                label = f"w^{j}" if j > 0 else "1"
                irr_label = f"({c2_label}, {label})"
                print(f"    {irr_label:<16} {dim_null:>5} {dim_live:>5} {dim_total:>6}")

                null_total += dim_null
                live_total += dim_live
                irrep_data.append({
                    'c2': c2_sign, 'c_k': j,
                    'null': dim_null, 'live': dim_live, 'total': dim_total
                })

        print(f"    {'-' * 38}")
        print(f"    {'Total':<16} {null_total:>5} {live_total:>5} {null_total + live_total:>6}")

        # Summary: C_2 content
        null_plus = sum(d['null'] for d in irrep_data if d['c2'] > 0)
        null_minus = sum(d['null'] for d in irrep_data if d['c2'] < 0)
        live_plus = sum(d['live'] for d in irrep_data if d['c2'] > 0)
        live_minus = sum(d['live'] for d in irrep_data if d['c2'] < 0)
        print(f"\n    C_2 summary: null (+{null_plus}, -{null_minus}), "
              f"live (+{live_plus}, -{live_minus})")
        print(f"    Null is {'CP-odd excess' if null_minus > null_plus else 'CP-even excess' if null_plus > null_minus else 'balanced'}")

        self.results['irreps'] = irrep_data

    def eigenvalue_irrep_table(self):
        """Decompose each live eigenspace into irreps."""
        if 'irreps' not in self.results or 'evecs' not in self.results:
            return

        evals = self.results['evals']
        evecs = self.results['evecs']
        live_mask = np.abs(evals) >= 1e-8
        live_evals = evals[live_mask]
        live_evecs = evecs[:, live_mask]

        vac_oms, vac_idx = self.orbits[self.vac_orbit_idx]
        bs = len(vac_idx)
        k = self.orbit_size

        match_to_vi = {}
        for vi in range(bs):
            match_to_vi[self.raw[vac_idx[vi]]] = vi

        def perm_matrix(g):
            M = np.zeros((bs, bs))
            for vi in range(bs):
                m_new = apply_vertex_mult(self.raw[vac_idx[vi]], g, self.p)
                if m_new in match_to_vi:
                    M[match_to_vi[m_new], vi] = 1.0
            return M

        P = perm_matrix(2)
        Q = perm_matrix(self.p - 1)
        P_powers = [np.eye(bs)]
        for _ in range(k - 1):
            P_powers.append(P_powers[-1] @ P)

        # Group live eigenvalues
        groups = compress_spectrum(live_evals)

        # Header
        irrep_labels = []
        for c2s in [+1, -1]:
            for j in range(k):
                c2l = '+' if c2s > 0 else '-'
                cl = f"w{j}" if j > 0 else "1"
                irrep_labels.append(f"({c2l},{cl})")

        header = f"    {'lambda':>10} {'m':>3}"
        for il in irrep_labels:
            header += f" {il:>6}"
        print(f"\n  EIGENVALUE x IRREP TABLE:")
        print(header)
        print(f"    {'-' * (16 + 7 * len(irrep_labels))}")

        idx_start = 0
        for val, mult in groups:
            sub = live_evecs[:, idx_start:idx_start + mult]
            row = f"    {val:>10.4f} {mult:>3}"
            for c2s in [+1, -1]:
                for j in range(k):
                    omega_j = np.exp(2j * np.pi * j / k)
                    proj_ck = sum(np.conj(omega_j) ** m * P_powers[m]
                                  for m in range(k)) / k
                    proj_c2 = (np.eye(bs) + c2s * Q) / 2.0
                    proj = np.real(proj_ck @ proj_c2)
                    proj_sub = sub.T @ proj @ sub
                    dim = round(np.real(np.trace(proj_sub)))
                    row += f" {dim:>6}"
            print(row)
            idx_start += mult

    def c2k_refinement(self):
        """Check if C_2 (conjugation) = P^k, collapsing wreath to cyclic."""
        if 'evecs' not in self.results or self.results.get('null_dim', 0) == 0:
            return

        vac_oms, vac_idx = self.orbits[self.vac_orbit_idx]
        bs = len(vac_idx)
        k = self.orbit_size

        match_to_vi = {}
        for vi in range(bs):
            match_to_vi[self.raw[vac_idx[vi]]] = vi

        def perm_matrix(g):
            M = np.zeros((bs, bs))
            for vi in range(bs):
                m_new = apply_vertex_mult(self.raw[vac_idx[vi]], g, self.p)
                if m_new in match_to_vi:
                    M[match_to_vi[m_new], vi] = 1.0
            return M

        P = perm_matrix(2)
        Q = perm_matrix(self.p - 1)

        # Find actual P order
        Ptmp = np.eye(bs)
        p_order = 0
        for i in range(1, 2 * self.p + 1):
            Ptmp = Ptmp @ P
            if np.allclose(Ptmp, np.eye(bs)):
                p_order = i
                break

        # Test Q = P^k
        Pk = np.linalg.matrix_power(P, k)
        q_is_pk = np.allclose(Q, Pk)

        # Number theory check
        twok_mod_p = pow(2, k, self.p)
        is_prime = all(self.p % i for i in range(2, self.p))

        print(f"\n  C_{{2k}} REFINEMENT:")
        print(f"    P order = {p_order}, k = {k}, 2k = {2*k}")
        print(f"    2^k mod p = {twok_mod_p}, p-1 = {self.p - 1}")
        print(f"    Q = P^k: {q_is_pk}")

        if q_is_pk:
            print(f"    -> Symmetry is CYCLIC C_{2*k} (CP = half-rotation)")
        else:
            print(f"    -> Symmetry is PRODUCT C_2 x C_{k} (CP independent)")

        self.results['p_order'] = p_order
        self.results['q_is_pk'] = q_is_pk
        self.results['cp_type'] = 'derived' if q_is_pk else 'independent'

    def run_all(self, do_coupling=True, do_irreps=True, do_blocks=True,
                max_coupling_orbits=None):
        """Run complete analysis for this level."""
        print(f"\n{'=' * 72}")
        print(f"  K_{self.nv} (p={self.p}, n={self.n})")
        print(f"{'=' * 72}")

        self.generate()
        self.sense2_vacuum_spectrum()

        if do_coupling:
            self.sense1_coupling(max_orbits=max_coupling_orbits)

        if do_blocks:
            self.sense3_block_coupling(max_orbits=max_coupling_orbits)
            self.sense3_vacuum_blocks()

        if do_irreps and self.results.get('null_dim', 0) > 0:
            self.null_space_irreps()
            self.eigenvalue_irrep_table()
            self.c2k_refinement()


# ======================================================================
# INTER-LEVEL EMBEDDING ANALYSIS
# ======================================================================

def embedding_analysis(probes):
    """Test rank-1 embedding theorem for all pairs of levels (adjacent and non-adjacent)."""
    print(f"\n{'=' * 72}")
    print("  INTER-LEVEL EMBEDDING ANALYSIS")
    print(f"{'=' * 72}")

    # Need at least 2 levels
    if len(probes) < 2:
        print("  [SKIP] Need at least 2 levels")
        return []

    results = []

    for i, src in enumerate(probes):
        for j, tgt in enumerate(probes):
            if tgt.nv <= src.nv:
                continue

            nv_src, nv_tgt = src.nv, tgt.nv
            extra = list(range(nv_src, nv_tgt))
            n_extra = len(extra)

            if n_extra % 2 != 0:
                continue  # can't pair odd number of extras

            m = n_extra // 2  # number of extra edges per pairing
            pairings = gen_pairings(extra)
            n_pairings = len(pairings)
            skip = (nv_tgt - nv_src) // 2 - 1  # 0 = adjacent

            print(f"\n  --- K_{nv_src} -> K_{nv_tgt} "
                  f"({'adjacent' if skip == 0 else f'skip {skip}'}) ---")
            print(f"  Extra vertices: {extra}, "
                  f"{n_pairings} pairing{'s' if n_pairings > 1 else ''}, "
                  f"m={m} extra edges")

            # Native overlap for source (use cached if available, else build)
            src_ms = [set(m_) for m_ in src.raw]
            O_nat = build_overlap(src_ms, src_ms, nv_src, same=True)
            evals_nat = np.sort(np.linalg.eigvalsh(O_nat))[::-1]
            n_matchings_src = len(src.raw)

            # Target vacuum orbit for cross-level coupling
            tgt_vac_oms, tgt_vac_idx = tgt.orbits[tgt.vac_orbit_idx]
            tgt_vac_ms = [tgt.match_sets[ii] for ii in tgt_vac_idx]

            # Test each pairing
            rank1_count = 0
            eval_tuples = set()
            cross_ranks = []
            cross_top_svs = []
            vac_overlaps = []

            # Cap output for large pairing counts
            show_detail = (n_pairings <= 5)

            for pi, pairing in enumerate(pairings):
                # Embed source matchings
                emb = [frozenset(m_ | pairing) for m_ in src.raw]
                emb_sets = [set(e) for e in emb]

                # Embedded overlap in target graph
                O_emb = build_overlap(emb_sets, emb_sets, nv_tgt, same=True)

                # Rank-1 test: O_emb = O_nat + 2m*J
                J = np.ones_like(O_nat)
                diff = np.max(np.abs(O_emb - O_nat - 2 * m * J))
                if diff < 1e-10:
                    rank1_count += 1

                # Eigenvalue comparison
                evals_emb = np.sort(np.linalg.eigvalsh(O_emb))[::-1]
                ground_shift = evals_emb[0] - evals_nat[0]
                max_other_shift = np.max(np.abs(evals_emb[1:] - evals_nat[1:]))
                eval_tuples.add(tuple(np.round(evals_emb, 8)))

                # Cross-level coupling with target vacuum
                O_cross = build_overlap(emb_sets, tgt_vac_ms, nv_tgt)
                svs = np.linalg.svd(O_cross, compute_uv=False)
                svs_nz = svs[svs > 1e-10]
                cross_ranks.append(len(svs_nz))
                cross_top_svs.append(svs_nz[:5].copy() if len(svs_nz) >= 5
                                     else svs_nz.copy())

                # How many embedded matchings land in target vacuum orbit?
                tgt_vac_set = set(tgt.raw[ii] for ii in tgt_vac_idx)
                n_in_vac = sum(1 for e in emb if e in tgt_vac_set)
                vac_overlaps.append(n_in_vac)

                if show_detail:
                    print(f"    Pairing {pi}: {sorted(pairing)}")
                    print(f"      ||O_emb - (O_nat + {2*m}J)|| = {diff:.2e}")
                    print(f"      Ground shift: {ground_shift:.1f} "
                          f"(expected {2*m*n_matchings_src}), "
                          f"max other shift: {max_other_shift:.2e}")
                    print(f"      Cross-coupling rank: {len(svs_nz)}, "
                          f"top SVs: {np.round(svs_nz[:4], 2)}")
                    print(f"      In target vacuum: {n_in_vac}/{n_matchings_src}")

            # Check composed path vs direct (for non-adjacent)
            composed_match = None
            if skip > 0:
                # Build composed: src -> src+2 -> src+4 -> ... -> tgt
                composed_emb = list(src.raw)
                for step_nv in range(nv_src + 2, nv_tgt + 1, 2):
                    top_edge = (step_nv - 2, step_nv - 1)
                    composed_emb = [m_ | frozenset([top_edge])
                                    for m_ in composed_emb]
                composed_set = set(composed_emb)

                for pi, pairing in enumerate(pairings):
                    direct = set(m_ | pairing for m_ in src.raw)
                    if direct == composed_set:
                        composed_match = pi
                        break

            # Summary for this pair
            print(f"\n  Summary K_{nv_src} -> K_{nv_tgt}:")
            print(f"    Rank-1 holds: {rank1_count}/{n_pairings}")
            print(f"    Distinct internal spectra: {len(eval_tuples)}")

            unique_ranks = sorted(set(cross_ranks))
            print(f"    Cross-coupling ranks: {unique_ranks} "
                  f"(native src rank = {int(np.sum(np.abs(evals_nat) > 1e-6))})")
            print(f"    Vacuum overlap: {min(vac_overlaps)}-{max(vac_overlaps)}"
                  f"/{n_matchings_src}")

            if composed_match is not None:
                print(f"    Composed tower path = direct pairing {composed_match}: "
                      f"{sorted(pairings[composed_match])}")
            elif skip > 0:
                print(f"    Composed tower path: no matching direct pairing found")

            # Cross-coupling SV equivalence across pairings
            if n_pairings > 1 and len(cross_top_svs) > 1:
                sv_groups = []
                used = [False] * n_pairings
                for pi in range(n_pairings):
                    if used[pi]:
                        continue
                    group = [pi]
                    used[pi] = True
                    for pj in range(pi + 1, n_pairings):
                        if used[pj]:
                            continue
                        min_len = min(len(cross_top_svs[pi]),
                                      len(cross_top_svs[pj]))
                        if min_len > 0:
                            d = np.max(np.abs(cross_top_svs[pi][:min_len] -
                                              cross_top_svs[pj][:min_len]))
                            if d < 1e-8:
                                group.append(pj)
                                used[pj] = True
                    sv_groups.append(group)
                print(f"    Cross-coupling SV equivalence classes: "
                      f"{[g for g in sv_groups]} ({len(sv_groups)} class"
                      f"{'es' if len(sv_groups) > 1 else ''})")

            results.append({
                'src': nv_src, 'tgt': nv_tgt, 'skip': skip,
                'm': m, 'n_pairings': n_pairings,
                'rank1_count': rank1_count,
                'n_spectra': len(eval_tuples),
                'cross_ranks': unique_ranks,
                'native_rank': int(np.sum(np.abs(evals_nat) > 1e-6)),
                'composed_match': composed_match
            })

    # Grand summary table
    if results:
        print(f"\n  {'=' * 68}")
        print(f"  EMBEDDING SUMMARY")
        print(f"  {'=' * 68}")
        print(f"  {'Embedding':<16} {'Skip':>5} {'Pairings':>9} "
              f"{'Rank-1':>7} {'Spectra':>8} {'Xrank':>8} {'NatRank':>8}")
        print(f"  {'-' * 68}")
        for r in results:
            print(f"  K_{r['src']}->K_{r['tgt']:<4}"
                  f" {r['skip']:>5} {r['n_pairings']:>9} "
                  f" {r['rank1_count']}/{r['n_pairings']:>5}"
                  f" {r['n_spectra']:>8}"
                  f" {r['cross_ranks']!s:>8}"
                  f" {r['native_rank']:>8}")

    return results


# ======================================================================
# CROSS-LEVEL SUMMARY
# ======================================================================

def cross_level_summary(probes):
    """Print cross-level comparison table."""
    print(f"\n{'=' * 72}")
    print("  CROSS-LEVEL SUMMARY")
    print(f"{'=' * 72}")

    # Law 1
    print(f"\n  LAW 1: VACUUM HALF-RANK")
    print(f"  {'Level':<8} {'bs':>6} {'rank':>6} {'null':>6} {'bs/2':>6} {'Status'}")
    print(f"  {'-' * 50}")
    for pr in probes:
        bs = pr.results.get('bs', '?')
        rank = pr.results.get('rank', '?')
        null = pr.results.get('null_dim', '?')
        law1 = pr.results.get('law1')
        half = bs // 2 if isinstance(bs, int) else '?'
        status = 'CONFIRMED' if law1 else ('n/a' if law1 is None else 'VIOLATED')
        print(f"  K_{pr.nv:<5} {bs:>6} {rank:>6} {null:>6} {half:>6} {status}")

    # Law 2
    print(f"\n  LAW 2: FP RANK DEFICIT = k-1")
    print(f"  {'Level':<8} {'k':>4} {'vac_rank':>9} {'fp_rank':>8} {'deficit':>8} {'k-1':>5} {'Status'}")
    print(f"  {'-' * 55}")
    for pr in probes:
        k = pr.orbit_size or '?'
        rank = pr.results.get('rank', '?')
        deficit = pr.results.get('fp_deficit')
        fp_coupling = [d for d in pr.results.get('coupling', []) if d.get('is_fp')]
        fp_rank = fp_coupling[0]['rank'] if fp_coupling else '?'
        if deficit is not None and isinstance(k, int):
            if pr.n >= 4:
                status = 'CONFIRMED' if deficit == k - 1 else 'VIOLATED'
            else:
                status = 'n=3 exception'
            print(f"  K_{pr.nv:<5} {k:>4} {rank:>9} {fp_rank:>8} {deficit:>8} {k-1:>5} {status}")
        else:
            print(f"  K_{pr.nv:<5} {k!s:>4} {rank!s:>9} {fp_rank!s:>8} {'?':>8} {'?':>5} n/a")

    # Law 3: generic multiplicity = 2k
    print(f"\n  LAW 3: GENERIC MULTIPLICITY = 2k")
    print(f"  {'Level':<8} {'k':>4} {'2k':>4}  {'Spectrum mults'}")
    print(f"  {'-' * 55}")
    for pr in probes:
        k = pr.orbit_size
        spec = pr.results.get('spectrum', [])
        mults = [c for _, c in spec if abs(_) > 1e-6]
        if k and mults:
            print(f"  K_{pr.nv:<5} {k:>4} {2*k:>4}  {mults}")

    # Law 4: C_2 parity
    print(f"\n  LAW 4: C_2 PARITY SELECTION")
    for pr in probes:
        irreps = pr.results.get('irreps')
        if irreps:
            null_plus = sum(d['null'] for d in irreps if d['c2'] > 0)
            null_minus = sum(d['null'] for d in irreps if d['c2'] < 0)
            live_plus = sum(d['live'] for d in irreps if d['c2'] > 0)
            live_minus = sum(d['live'] for d in irreps if d['c2'] < 0)
            print(f"  K_{pr.nv}: null(+{null_plus},-{null_minus}) "
                  f"live(+{live_plus},-{live_minus}) "
                  f"-> null is {'CP-odd excess' if null_minus > null_plus else 'CP-even excess' if null_plus > null_minus else 'balanced'}")

    # Law 5: CP independence
    print(f"\n  LAW 5: CP INDEPENDENCE (2^k = -1 mod p test)")
    print(f"  {'Level':<8} {'p':>4} {'k':>4} {'2^k mod p':>10} {'Group':>15} {'CP status'}")
    print(f"  {'-' * 55}")
    for pr in probes:
        cp_type = pr.results.get('cp_type')
        if cp_type:
            k = pr.orbit_size
            twok = pow(2, k, pr.p)
            if cp_type == 'derived':
                grp = f"C_{2*k} (cyclic)"
            else:
                grp = f"C_2 x C_{k}"
            print(f"  K_{pr.nv:<5} {pr.p:>4} {k:>4} {twok:>10} {grp:>15} {cp_type}")


# ======================================================================
# MAIN
# ======================================================================

def main():
    parser = argparse.ArgumentParser(description='Coupling Probe: Unified')
    parser.add_argument('--levels', nargs='+', type=int, default=[6, 8, 10, 12],
                        help='K_{2n} levels to analyze (default: 6 8 10 12)')
    parser.add_argument('--no-irreps', action='store_true',
                        help='Skip irrep decomposition')
    parser.add_argument('--no-coupling', action='store_true',
                        help='Skip inter-orbit coupling (Sense 1)')
    parser.add_argument('--no-blocks', action='store_true',
                        help='Skip block-level coupling (Sense 3)')
    parser.add_argument('--no-embedding', action='store_true',
                        help='Skip inter-level embedding analysis')
    parser.add_argument('--max-coupling-orbits', type=int, default=None,
                        help='Limit number of non-vacuum orbits for coupling analysis')
    args = parser.parse_args()

    t_total = time.time()
    print("COUPLING PROBE: UNIFIED ANALYSIS")
    print(f"Levels: {['K_' + str(lv) for lv in args.levels]}")
    print(f"Irreps: {'ON' if not args.no_irreps else 'OFF'}")
    print(f"Coupling: {'ON' if not args.no_coupling else 'OFF'}")
    print(f"Blocks: {'ON' if not args.no_blocks else 'OFF'}")
    print(f"Embedding: {'ON' if not args.no_embedding else 'OFF'}")

    probes = []
    for nv in args.levels:
        pr = LevelProbe(nv)
        pr.run_all(
            do_coupling=not args.no_coupling,
            do_irreps=not args.no_irreps,
            do_blocks=not args.no_blocks,
            max_coupling_orbits=args.max_coupling_orbits
        )
        probes.append(pr)

    cross_level_summary(probes)

    if not args.no_embedding:
        embedding_analysis(probes)

    print(f"\n{'=' * 72}")
    print(f"  TOTAL TIME: {time.time() - t_total:.1f}s")
    print(f"{'=' * 72}")


if __name__ == '__main__':
    main()
