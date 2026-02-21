# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Pfaffian Orientation Probe (Sense 8)
=====================================
Every perfect matching contributes to the Pfaffian with sign ±1.
This probe asks: how does the fermion sign correlate with everything
we already know about the matching landscape?

Probes:
  8a: Sign distribution per orbit (balanced? imbalanced?)
  8b: Sign vector projected onto null/live eigenspaces
  8c: Sign correlation with C_2 (CP) parity
  8d: Sign-stratified overlap: separate ± matchings, compare sub-block spectra
  8e: Signed overlap matrix (element-wise sign product) vs unsigned
  8f: Sign flow under torus action (does ×2 mod p flip signs?)

Usage:
  python pfaffian_probe.py [--levels 6 8 10 12]
"""

import numpy as np
from collections import defaultdict
import time
import argparse

# ======================================================================
# IMPORT CORE FROM COUPLING PROBE
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


def build_overlap(ms_a, ms_b, nv, same=False):
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


def apply_vertex_mult(matching, g, p):
    new = set()
    for a, b in matching:
        na = (g * a) % p if a < p else a
        nb = (g * b) % p if b < p else b
        new.add((min(na, nb), max(na, nb)))
    return frozenset(new)


# ======================================================================
# PFAFFIAN SIGN COMPUTATION
# ======================================================================

def perm_sign(perm):
    """Sign of a permutation given as a list.
    
    Uses cycle decomposition: sign = (-1)^(n - number_of_cycles).
    """
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
    """Compute the Pfaffian sign of a perfect matching.
    
    Given matching M = {(a₁,b₁), ..., (aₙ,bₙ)} with aₖ < bₖ,
    sorted so a₁ < a₂ < ... < aₙ, the sign is the sign of the
    permutation π = (a₁, b₁, a₂, b₂, ..., aₙ, bₙ).
    """
    # Sort edges by first vertex
    edges = sorted(matching, key=lambda e: (min(e), max(e)))
    edges = [(min(a, b), max(a, b)) for a, b in edges]
    
    # Build the permutation as a list
    perm_list = []
    for a, b in edges:
        perm_list.append(a)
        perm_list.append(b)
    
    # perm_list is now a permutation of [0, 1, ..., nv-1]
    # We need the sign of this permutation
    return perm_sign(perm_list)


def compute_all_signs(raw_matchings, nv):
    """Compute Pfaffian signs for all matchings."""
    signs = []
    for m in raw_matchings:
        signs.append(pfaffian_sign(m, nv))
    return np.array(signs, dtype=float)


# ======================================================================
# SENSE 8 PROBE
# ======================================================================

class PfaffianProbe:
    """Full Pfaffian orientation analysis for a single K_{2n} level."""
    
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
        self.results = {}
    
    def generate(self):
        """Step 1: Generate matchings, orbits, and signs."""
        t0 = time.time()
        self.raw = gen_matchings(self.nv)
        self.match_sets = [set(m) for m in self.raw]
        N = len(self.raw)
        
        # Compute Pfaffian signs
        self.signs = compute_all_signs(self.raw, self.nv)
        n_plus = int(np.sum(self.signs > 0))
        n_minus = int(np.sum(self.signs < 0))
        
        # Build orbit structure
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
        
        # Identify vacuum
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
        
        dt = time.time() - t0
        print(f"\n  {N} matchings: {n_plus} (+), {n_minus} (-)")
        print(f"  Global sign balance: {'EXACT' if n_plus == n_minus else f'{n_plus}:{n_minus}'}")
        print(f"  {len(self.blocks)} blocks, {len(self.orbits)} orbits ({dt:.1f}s)")
        
        self.results['N'] = N
        self.results['n_plus'] = n_plus
        self.results['n_minus'] = n_minus
    
    def probe_8a_orbit_sign_distribution(self):
        """8a: Sign distribution within each orbit."""
        print(f"\n  ── PROBE 8a: SIGN DISTRIBUTION PER ORBIT ──")
        
        orbit_sign_data = []
        
        for oi, (oms, idx) in enumerate(self.orbits):
            s = self.signs[idx]
            n_p = int(np.sum(s > 0))
            n_m = int(np.sum(s < 0))
            total = len(idx)
            net = n_p - n_m
            balanced = (n_p == n_m)
            
            tag = ""
            if oi == self.vac_orbit_idx:
                tag = " << VACUUM"
            elif oi == self.fp_orbit_idx:
                tag = " (FP)"
            
            imbalance = abs(net) / total if total > 0 else 0
            print(f"    Orbit {oi:>2} (ms={total:>6}): "
                  f"(+){n_p:>5} (-){n_m:>5}  net={net:>+5}  "
                  f"{'BALANCED' if balanced else f'imbalance={imbalance:.4f}'}{tag}")
            
            # Block-level detail for vacuum
            if oi == self.vac_orbit_idx:
                print(f"      Block-level detail:")
                for bi, ms_key in enumerate(oms):
                    blk_idx = self.blocks[ms_key]
                    bs = self.signs[blk_idx]
                    bp = int(np.sum(bs > 0))
                    bm = int(np.sum(bs < 0))
                    print(f"        Block {bi} (dir={ms_key}): (+){bp} (-){bm} net={bp-bm:+d}")
            
            orbit_sign_data.append({
                'orbit': oi, 'total': total,
                'n_plus': n_p, 'n_minus': n_m,
                'net': net, 'balanced': balanced
            })
        
        n_balanced = sum(1 for d in orbit_sign_data if d['balanced'])
        print(f"\n    Balanced orbits: {n_balanced}/{len(orbit_sign_data)}")
        
        # The Pfaffian is the sum of signed matchings (for the adjacency matrix)
        # Net sign per orbit tells us which orbits contribute to Pf(A)
        net_contributions = [(d['orbit'], d['net']) for d in orbit_sign_data if d['net'] != 0]
        if net_contributions:
            print(f"    Orbits with net Pfaffian contribution: {net_contributions}")
        else:
            print(f"    ALL orbits sign-balanced → each orbit contributes zero to Pfaffian")
        
        self.results['orbit_signs'] = orbit_sign_data
    
    def probe_8b_sign_projection(self):
        """8b: Project sign vector onto eigenspaces of vacuum overlap."""
        if self.vac_orbit_idx is None:
            return
        
        print(f"\n  ── PROBE 8b: SIGN VECTOR EIGENSPACE PROJECTION ──")
        
        vac_oms, vac_idx = self.orbits[self.vac_orbit_idx]
        vac_ms = [self.match_sets[i] for i in vac_idx]
        bs = len(vac_idx)
        
        # Build vacuum overlap and diagonalize
        O = build_overlap(vac_ms, vac_ms, self.nv, same=True)
        evals, evecs = np.linalg.eigh(O)
        order = np.argsort(evals)[::-1]
        evals = evals[order]
        evecs = evecs[:, order]
        
        # Sign vector restricted to vacuum orbit
        sign_vec = self.signs[vac_idx]
        sign_norm = np.linalg.norm(sign_vec)
        
        print(f"    Vacuum sign vector: {int(np.sum(sign_vec > 0))} (+), "
              f"{int(np.sum(sign_vec < 0))} (-), ||s|| = {sign_norm:.4f}")
        
        # Null/live split
        null_mask = np.abs(evals) < 1e-8
        live_mask = ~null_mask
        
        null_vecs = evecs[:, null_mask]
        live_vecs = evecs[:, live_mask]
        
        # Project sign vector
        if null_vecs.shape[1] > 0:
            sign_null = null_vecs.T @ sign_vec
            sign_null_norm = np.linalg.norm(sign_null)
        else:
            sign_null_norm = 0.0
        
        sign_live = live_vecs.T @ sign_vec
        sign_live_norm = np.linalg.norm(sign_live)
        
        total_norm_sq = sign_null_norm**2 + sign_live_norm**2
        null_frac = sign_null_norm**2 / total_norm_sq if total_norm_sq > 0 else 0
        live_frac = sign_live_norm**2 / total_norm_sq if total_norm_sq > 0 else 0
        
        print(f"    ||s_null||² = {sign_null_norm**2:.6f} ({null_frac:.4%})")
        print(f"    ||s_live||² = {sign_live_norm**2:.6f} ({live_frac:.4%})")
        print(f"    ||s||² = {total_norm_sq:.6f} (check: {sign_norm**2:.6f})")
        
        if null_frac > 0.99:
            print(f"    *** SIGN IS PURELY NULL-SPACE: fermion phase is a silent mode ***")
        elif null_frac < 0.01:
            print(f"    *** SIGN IS PURELY LIVE-SPACE: fermion phase is fully dynamical ***")
        elif abs(null_frac - 0.5) < 0.01:
            print(f"    *** SIGN IS HALF-NULL, HALF-LIVE: equidistributed across sectors ***")
        
        # Project onto individual eigenspaces
        print(f"\n    Eigenspace-resolved sign content:")
        groups = compress_spectrum(evals)
        idx_start = 0
        for val, mult in groups:
            sub_vecs = evecs[:, idx_start:idx_start + mult]
            proj_coeff = sub_vecs.T @ sign_vec
            proj_norm_sq = np.sum(proj_coeff**2)
            frac = proj_norm_sq / total_norm_sq if total_norm_sq > 0 else 0
            tag = " [null]" if abs(val) < 1e-6 else ""
            if frac > 0.001:
                print(f"      λ = {val:>10.4f} (x{mult}): "
                      f"||proj||² = {proj_norm_sq:>8.4f} ({frac:>6.2%}){tag}")
            idx_start += mult
        
        # Store for later probes
        self.results['O_vac'] = O
        self.results['evals'] = evals
        self.results['evecs'] = evecs
        self.results['sign_vec'] = sign_vec
        self.results['null_frac'] = null_frac
        self.results['live_frac'] = live_frac
        self.results['bs'] = bs
    
    def probe_8c_sign_cp_correlation(self):
        """8c: Does the Pfaffian sign correlate with C_2 (CP) parity?"""
        if 'evecs' not in self.results:
            return
        if self.results.get('bs', 0) < 4:
            return
        
        print(f"\n  ── PROBE 8c: SIGN × CP CORRELATION ──")
        
        vac_oms, vac_idx = self.orbits[self.vac_orbit_idx]
        bs = self.results['bs']
        k = self.orbit_size
        
        # Build permutation matrices
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
        
        P = perm_matrix(2)  # Torus generator ×2
        Q = perm_matrix(self.p - 1)  # CP: ×(-1) mod p
        
        sign_vec = self.results['sign_vec']
        
        # Does CP flip the sign? Q·s = ±s?
        Qs = Q @ sign_vec
        
        # Check if Qs = +s or Qs = -s
        plus_match = np.allclose(Qs, sign_vec)
        minus_match = np.allclose(Qs, -sign_vec)
        
        if plus_match:
            print(f"    Q·s = +s: sign vector is CP-EVEN")
            self.results['sign_cp'] = 'even'
        elif minus_match:
            print(f"    Q·s = -s: sign vector is CP-ODD")
            self.results['sign_cp'] = 'odd'
        else:
            # Decompose into CP-even and CP-odd parts
            s_even = (sign_vec + Qs) / 2
            s_odd = (sign_vec - Qs) / 2
            norm_even = np.linalg.norm(s_even)**2
            norm_odd = np.linalg.norm(s_odd)**2
            total = norm_even + norm_odd
            print(f"    Sign vector is MIXED CP:")
            print(f"      CP-even component: ||s_+||² = {norm_even:.4f} ({norm_even/total:.2%})")
            print(f"      CP-odd component:  ||s_-||² = {norm_odd:.4f} ({norm_odd/total:.2%})")
            self.results['sign_cp'] = 'mixed'
            self.results['sign_cp_even_frac'] = norm_even / total
        
        # Does the torus action preserve the sign?
        Ps = P @ sign_vec
        torus_preserves = np.allclose(Ps, sign_vec)
        torus_flips = np.allclose(Ps, -sign_vec)
        
        if torus_preserves:
            print(f"    P·s = +s: sign is TORUS-INVARIANT")
        elif torus_flips:
            print(f"    P·s = -s: torus FLIPS sign (×2 reverses orientation)")
        else:
            # How much does torus change the sign?
            dot = sign_vec @ Ps
            cos_angle = dot / (np.linalg.norm(sign_vec) * np.linalg.norm(Ps))
            print(f"    P·s ≠ ±s: sign partially rotated by torus")
            print(f"      cos(angle(s, P·s)) = {cos_angle:.6f}")
            
            # Track sign through full torus orbit
            cur = sign_vec.copy()
            torus_signs = [1.0]  # dot product with original
            for step in range(1, k):
                cur = P @ cur
                torus_signs.append(sign_vec @ cur / (np.linalg.norm(sign_vec)**2))
            print(f"      Torus orbit of sign: {[f'{x:.4f}' for x in torus_signs]}")
        
        self.results['torus_preserves_sign'] = torus_preserves
        self.results['torus_flips_sign'] = torus_flips
        
        # Joint sign × CP × null analysis
        if 'evecs' in self.results:
            evals = self.results['evals']
            evecs = self.results['evecs']
            
            null_mask = np.abs(evals) < 1e-8
            null_vecs = evecs[:, null_mask]
            live_vecs = evecs[:, ~null_mask]
            
            if null_vecs.shape[1] > 0:
                # Project CP-even and CP-odd parts of sign into null/live
                s_even = (sign_vec + Qs) / 2
                s_odd = (sign_vec - Qs) / 2
                
                for label, svec in [("CP-even", s_even), ("CP-odd", s_odd)]:
                    if np.linalg.norm(svec) < 1e-10:
                        continue
                    null_proj = np.linalg.norm(null_vecs.T @ svec)**2
                    live_proj = np.linalg.norm(live_vecs.T @ svec)**2
                    total = null_proj + live_proj
                    if total > 1e-10:
                        print(f"    {label} sign: null={null_proj/total:.2%}, live={live_proj/total:.2%}")
    
    def probe_8d_sign_stratified_overlap(self):
        """8d: Separate ± matchings within vacuum, compare sub-block structure."""
        if self.vac_orbit_idx is None:
            return
        
        print(f"\n  ── PROBE 8d: SIGN-STRATIFIED OVERLAP ──")
        
        vac_oms, vac_idx = self.orbits[self.vac_orbit_idx]
        sign_vec = self.results.get('sign_vec')
        if sign_vec is None:
            return
        
        # Split vacuum matchings by sign
        plus_mask = sign_vec > 0
        minus_mask = sign_vec < 0
        n_plus = int(np.sum(plus_mask))
        n_minus = int(np.sum(minus_mask))
        
        plus_local = np.where(plus_mask)[0]
        minus_local = np.where(minus_mask)[0]
        
        plus_ms = [self.match_sets[vac_idx[i]] for i in plus_local]
        minus_ms = [self.match_sets[vac_idx[i]] for i in minus_local]
        
        print(f"    Vacuum: {n_plus} (+) and {n_minus} (-) matchings")
        
        # Build 4 sub-blocks
        if n_plus > 0 and n_minus > 0:
            O_pp = build_overlap(plus_ms, plus_ms, self.nv, same=True)
            O_mm = build_overlap(minus_ms, minus_ms, self.nv, same=True)
            O_pm = build_overlap(plus_ms, minus_ms, self.nv)
            O_mp = build_overlap(minus_ms, plus_ms, self.nv)
            
            # Spectra of diagonal blocks
            evals_pp = np.sort(np.linalg.eigvalsh(O_pp))[::-1]
            evals_mm = np.sort(np.linalg.eigvalsh(O_mm))[::-1]
            
            rank_pp = int(np.sum(np.abs(evals_pp) > 1e-6))
            rank_mm = int(np.sum(np.abs(evals_mm) > 1e-6))
            null_pp = n_plus - rank_pp
            null_mm = n_minus - rank_mm
            
            print(f"\n    O(+,+): {n_plus}×{n_plus}, rank={rank_pp}, null={null_pp}")
            groups_pp = compress_spectrum(evals_pp)
            for v, c in groups_pp:
                if abs(v) > 1e-6 or c > 1:
                    print(f"      {v:>12.4f} (x{c})")
            
            print(f"\n    O(-,-): {n_minus}×{n_minus}, rank={rank_mm}, null={null_mm}")
            groups_mm = compress_spectrum(evals_mm)
            for v, c in groups_mm:
                if abs(v) > 1e-6 or c > 1:
                    print(f"      {v:>12.4f} (x{c})")
            
            # Are O(+,+) and O(-,-) isospectral?
            if n_plus == n_minus:
                max_diff = np.max(np.abs(evals_pp - evals_mm))
                isospectral = max_diff < 1e-6
                print(f"\n    O(+,+) ≅ O(-,-) isospectral: {isospectral} (max diff = {max_diff:.2e})")
                
                if isospectral:
                    print(f"    *** Sign sectors have IDENTICAL internal geometry ***")
                
                # Is O(+,+) = O(-,-) as matrices (after reordering)?
                # Check Frobenius norms
                frob_pp = np.sum(O_pp**2)
                frob_mm = np.sum(O_mm**2)
                frob_pm = np.sum(O_pm**2)
                frob_mp = np.sum(O_mp**2)
                print(f"    ||O(+,+)||² = {frob_pp:.0f}")
                print(f"    ||O(-,-)||² = {frob_mm:.0f}")
                print(f"    ||O(+,-)||² = {frob_pm:.0f}")
                print(f"    ||O(-,+)||² = {frob_mp:.0f}")
            
            # Cross-sign coupling via SVD
            svs = np.linalg.svd(O_pm, compute_uv=False)
            svs_nz = svs[svs > 1e-6]
            print(f"\n    Cross-sign coupling O(+,-): rank={len(svs_nz)}")
            sg = compress_spectrum(svs_nz)
            for sv, cnt in sg:
                print(f"      σ = {sv:10.4f} (x{cnt})")
            
            self.results['sign_strat'] = {
                'n_plus': n_plus, 'n_minus': n_minus,
                'rank_pp': rank_pp, 'rank_mm': rank_mm,
                'null_pp': null_pp, 'null_mm': null_mm,
                'isospectral': n_plus == n_minus and np.max(np.abs(evals_pp - evals_mm)) < 1e-6,
                'cross_rank': len(svs_nz)
            }
        else:
            print(f"    Only one sign class in vacuum — no stratification possible")
    
    def probe_8e_signed_overlap_matrix(self):
        """8e: The signed overlap matrix S_{ij} = sign(i)*sign(j)*O_{ij}.
        
        Key insight: this IS a similarity transformation (D*O*D with D²=I),
        so eigenvalues are identical. But the EIGENVECTORS rotate — and that
        rotation relative to the orbit structure is the information.
        
        More interesting: the "sign-weighted mean field" ⟨O⟩_sign = Σ_j s_j O_{ij}
        which is O·s — the overlap matrix acting on the sign vector.
        """
        if 'O_vac' not in self.results:
            return
        
        print(f"\n  ── PROBE 8e: SIGNED OVERLAP ANALYSIS ──")
        
        O = self.results['O_vac']
        sign_vec = self.results['sign_vec']
        bs = self.results['bs']
        
        # O·s: how the overlap matrix acts on the sign vector
        Os = O @ sign_vec
        
        # This tells us: for each matching i, what is the "signed mean field"
        # it experiences from all other matchings?
        # If O·s ∝ s, then the sign vector is an eigenvector of O.
        
        # Check if sign is an eigenvector
        if np.linalg.norm(sign_vec) > 0:
            # Rayleigh quotient
            ray = (sign_vec @ Os) / (sign_vec @ sign_vec)
            residual = Os - ray * sign_vec
            residual_norm = np.linalg.norm(residual) / np.linalg.norm(Os)
            
            print(f"    O·s Rayleigh quotient: {ray:.6f}")
            print(f"    ||O·s - λs|| / ||O·s|| = {residual_norm:.6f}")
            
            if residual_norm < 0.01:
                print(f"    *** SIGN IS APPROXIMATE EIGENVECTOR of O with λ ≈ {ray:.4f} ***")
            elif residual_norm < 0.05:
                print(f"    Sign is near-eigenvector (5% residual)")
            else:
                print(f"    Sign is NOT an eigenvector of O")
            
            # Decompose O·s into eigenspace components
            evals = self.results['evals']
            evecs = self.results['evecs']
            
            Os_components = evecs.T @ Os
            Os_norm_sq = np.sum(Os_components**2)
            
            groups = compress_spectrum(evals)
            idx_start = 0
            print(f"\n    O·s eigenspace decomposition:")
            for val, mult in groups:
                comp = Os_components[idx_start:idx_start + mult]
                comp_norm_sq = np.sum(comp**2)
                frac = comp_norm_sq / Os_norm_sq if Os_norm_sq > 0 else 0
                if frac > 0.001:
                    tag = " [null]" if abs(val) < 1e-6 else ""
                    print(f"      λ = {val:>10.4f} (x{mult}): {frac:>6.2%}{tag}")
                idx_start += mult
        
        # The "Pfaffian inner product": s^T O s
        # This is the total signed overlap — the Pfaffian squared (for adjacency)
        pf_inner = sign_vec @ Os
        print(f"\n    s^T · O · s = {pf_inner:.4f}")
        print(f"    s^T · O · s / ||s||² = {pf_inner / (sign_vec @ sign_vec):.6f}")
        print(f"    (This is the sign-weighted total overlap ≈ Pfaffian structure constant)")
    
    def probe_8f_sign_flow(self):
        """8f: How does the sign transform under vertex permutations?
        
        Key question: does ×g mod p preserve or flip Pfaffian signs?
        This connects the combinatorial sign to the number theory of the surface.
        """
        if self.vac_orbit_idx is None:
            return
        
        print(f"\n  ── PROBE 8f: SIGN FLOW UNDER VERTEX MAPS ──")
        
        # Test various vertex maps: ×2, ×3, ..., ×(p-1) mod p
        # For each, check if it preserves/flips/mixes signs
        
        for g in range(2, self.p):
            # Apply ×g to all matchings and track sign changes
            sign_changes = []
            valid = True
            
            for i, m in enumerate(self.raw):
                m_new = apply_vertex_mult(m, g, self.p)
                # Find m_new in raw
                try:
                    j = self.raw.index(m_new)
                except ValueError:
                    valid = False
                    break
                sign_changes.append(self.signs[j] / self.signs[i])
            
            if not valid:
                continue
            
            sign_changes = np.array(sign_changes)
            all_preserve = np.all(sign_changes > 0)
            all_flip = np.all(sign_changes < 0)
            
            if all_preserve:
                det = "+1"
            elif all_flip:
                det = "-1"
            else:
                n_preserve = int(np.sum(sign_changes > 0))
                n_flip = int(np.sum(sign_changes < 0))
                det = f"mixed ({n_preserve}+, {n_flip}-)"
            
            # Is g a quadratic residue mod p?
            qr = any(pow(x, 2, self.p) == g % self.p for x in range(1, self.p))
            qr_tag = "QR" if qr else "NR"
            
            # Only print the interesting ones
            if g <= 5 or g == self.p - 1 or not (all_preserve or all_flip):
                print(f"    ×{g} mod {self.p} ({qr_tag}): sign → {det}")
        
        # Special focus: does the sign flow factorize as (−1)^{something combinatorial}?
        # Test: is sign(×g · m) = legendre(g,p)^n · sign(m)?
        print(f"\n    Testing Legendre symbol factorization:")
        
        def legendre(a, p):
            """Legendre symbol (a/p)."""
            ls = pow(a, (p - 1) // 2, p)
            return -1 if ls == p - 1 else ls
        
        for g in range(2, min(self.p, 8)):
            leg = legendre(g, self.p)
            predicted_det = leg ** self.n  # n = nv/2
            
            # Check actual
            actual_det_list = []
            for i, m in enumerate(self.raw):
                m_new = apply_vertex_mult(m, g, self.p)
                try:
                    j = self.raw.index(m_new)
                    actual_det_list.append(self.signs[j] / self.signs[i])
                except ValueError:
                    break
            
            if actual_det_list:
                all_same = all(np.isclose(x, actual_det_list[0]) for x in actual_det_list)
                actual = actual_det_list[0] if all_same else "mixed"
                matches = all_same and np.isclose(actual_det_list[0], predicted_det)
                print(f"    ×{g}: Legendre({g},{self.p})^{self.n} = ({leg})^{self.n} = {predicted_det}, "
                      f"actual = {actual}, {'MATCHES' if matches else 'DIFFERS'}")
    
    def probe_8g_direction_sign_decomposition(self):
        """8g: Sign decomposition by direction class.
        
        Within each direction block of the vacuum, is the sign distribution
        uniform or does it depend on which directions a matching uses?
        """
        if self.vac_orbit_idx is None:
            return
        
        print(f"\n  ── PROBE 8g: DIRECTION × SIGN STRUCTURE ──")
        
        vac_oms, vac_idx = self.orbits[self.vac_orbit_idx]
        k = self.orbit_size
        
        for bi, ms_key in enumerate(vac_oms):
            blk_idx = self.blocks[ms_key]
            blk_signs = self.signs[blk_idx]
            n_p = int(np.sum(blk_signs > 0))
            n_m = int(np.sum(blk_signs < 0))
            
            # Within this block, do + and - matchings have different
            # overlap statistics?
            plus_idx = [blk_idx[j] for j in range(len(blk_idx)) if blk_signs[j] > 0]
            minus_idx = [blk_idx[j] for j in range(len(blk_idx)) if blk_signs[j] < 0]
            
            if plus_idx and minus_idx:
                # Average intra-sign overlap vs inter-sign overlap
                plus_ms = [self.match_sets[i] for i in plus_idx]
                minus_ms = [self.match_sets[i] for i in minus_idx]
                
                O_pp = build_overlap(plus_ms, plus_ms, self.nv, same=True)
                O_mm = build_overlap(minus_ms, minus_ms, self.nv, same=True)
                O_pm = build_overlap(plus_ms, minus_ms, self.nv)
                
                # Average off-diagonal overlaps
                np_count = len(plus_idx)
                nm_count = len(minus_idx)
                
                avg_pp = (np.sum(O_pp) - np.trace(O_pp)) / max(np_count * (np_count - 1), 1)
                avg_mm = (np.sum(O_mm) - np.trace(O_mm)) / max(nm_count * (nm_count - 1), 1)
                avg_pm = np.mean(O_pm) if O_pm.size > 0 else 0
                
                assort = "ASSORTATIVE" if avg_pp > avg_pm and avg_mm > avg_pm else \
                         "DISASSORTATIVE" if avg_pp < avg_pm and avg_mm < avg_pm else \
                         "NEUTRAL"
                
                print(f"    Block {bi} (dir={ms_key}): (+){n_p} (-){n_m}  "
                      f"⟨O⟩(+,+)={avg_pp:.2f} ⟨O⟩(-,-)={avg_mm:.2f} "
                      f"⟨O⟩(+,-)={avg_pm:.2f}  {assort}")
            else:
                print(f"    Block {bi} (dir={ms_key}): (+){n_p} (-){n_m}  "
                      f"(single sign class)")
    
    def run_all(self):
        """Run complete Pfaffian probe."""
        print(f"\n{'=' * 72}")
        print(f"  PFAFFIAN PROBE: K_{self.nv} (p={self.p}, n={self.n})")
        print(f"{'=' * 72}")
        
        self.generate()
        self.probe_8a_orbit_sign_distribution()
        self.probe_8b_sign_projection()
        self.probe_8c_sign_cp_correlation()
        self.probe_8d_sign_stratified_overlap()
        self.probe_8e_signed_overlap_matrix()
        self.probe_8f_sign_flow()
        self.probe_8g_direction_sign_decomposition()


# ======================================================================
# CROSS-LEVEL SUMMARY
# ======================================================================

def cross_level_summary(probes):
    """Print cross-level comparison of Pfaffian structure."""
    print(f"\n{'=' * 72}")
    print(f"  PFAFFIAN PROBE: CROSS-LEVEL SUMMARY")
    print(f"{'=' * 72}")
    
    # Global sign balance
    print(f"\n  GLOBAL SIGN BALANCE:")
    print(f"  {'Level':<8} {'N':>7} {'(+)':>7} {'(-)':>7} {'Balance'}")
    print(f"  {'-' * 40}")
    for pr in probes:
        r = pr.results
        bal = 'EXACT' if r['n_plus'] == r['n_minus'] else f"{r['n_plus']}:{r['n_minus']}"
        print(f"  K_{pr.nv:<5} {r['N']:>7} {r['n_plus']:>7} {r['n_minus']:>7} {bal}")
    
    # Vacuum sign ↔ null/live
    print(f"\n  VACUUM SIGN → NULL/LIVE PROJECTION:")
    print(f"  {'Level':<8} {'null%':>8} {'live%':>8} {'Interpretation'}")
    print(f"  {'-' * 50}")
    for pr in probes:
        nf = pr.results.get('null_frac')
        lf = pr.results.get('live_frac')
        if nf is not None:
            if nf > 0.99:
                interp = "PURELY NULL"
            elif nf < 0.01:
                interp = "PURELY LIVE"
            elif abs(nf - 0.5) < 0.05:
                interp = "EQUIDISTRIBUTED"
            else:
                interp = f"mixed"
            print(f"  K_{pr.nv:<5} {nf:>7.2%} {lf:>7.2%}  {interp}")
    
    # Sign × CP
    print(f"\n  SIGN × CP PARITY:")
    print(f"  {'Level':<8} {'CP type':>12} {'Sign under CP':>15}")
    print(f"  {'-' * 40}")
    for pr in probes:
        cp = pr.results.get('sign_cp', 'n/a')
        print(f"  K_{pr.nv:<5} {'':>12} {cp:>15}")
    
    # Sign stratification
    print(f"\n  SIGN-STRATIFIED OVERLAP:")
    print(f"  {'Level':<8} {'n(+)':>6} {'n(-)':>6} {'iso?':>5} {'xrank':>6}")
    print(f"  {'-' * 40}")
    for pr in probes:
        ss = pr.results.get('sign_strat')
        if ss:
            iso = 'YES' if ss.get('isospectral') else 'NO'
            print(f"  K_{pr.nv:<5} {ss['n_plus']:>6} {ss['n_minus']:>6} "
                  f"{iso:>5} {ss['cross_rank']:>6}")


# ======================================================================
# MAIN
# ======================================================================

def main():
    parser = argparse.ArgumentParser(description='Pfaffian Orientation Probe (Sense 8)')
    parser.add_argument('--levels', nargs='+', type=int, default=[6, 8, 10, 12],
                        help='K_{2n} levels to analyze (default: 6 8 10 12)')
    args = parser.parse_args()
    
    t_total = time.time()
    print("PFAFFIAN ORIENTATION PROBE (SENSE 8)")
    print(f"Levels: {['K_' + str(lv) for lv in args.levels]}")
    print("=" * 72)
    
    probes = []
    for nv in args.levels:
        pr = PfaffianProbe(nv)
        pr.run_all()
        probes.append(pr)
    
    cross_level_summary(probes)
    
    print(f"\n{'=' * 72}")
    print(f"  TOTAL TIME: {time.time() - t_total:.1f}s")
    print(f"{'=' * 72}")


if __name__ == '__main__':
    main()
