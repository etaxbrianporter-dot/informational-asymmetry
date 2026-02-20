#!/usr/bin/env python3
"""
MATCHING EXTRACTOR — Capability 1 (v2): Correct Direction Permutation
=======================================================================

CRITICAL FIX: The direction permutation is NOT d → d+1 mod (n-1).
It's the action of ×2 mod p on direction classes, where p = 2n-1.

Why ×2? The Heawood embedding has cyclic Z_p symmetry (vertex rotation).
The MAP i → 2i mod p is an additional automorphism of the Heawood graph
that PERMUTES direction classes. Specifically:
  edge (a,b) with diff d → edge (2a, 2b) with diff 2d (mod p)

Since direction class of diff d is the same as diff (p-d), we work 
with unsigned residues: dir(d) = min(d, p-d) - 1.

The permutation ×2 mod p acting on direction classes:
  K₆ (p=5):  0→1→0  (2-cycle, happens to equal d+1 mod 2)
  K₈ (p=7):  0→1→2→0  (3-cycle, happens to equal d+1 mod 3)
  K₁₂ (p=11): 0→1→3→2→4→0  (5-cycle, NOT d+1 mod 5)
  K₁₄ (p=13): depends on ×2 mod 13

For composite p (K₁₀, p=9): ×2 mod 9 doesn't generate a full cycle
on direction classes because some residues share factors with 9.
"""

import numpy as np
from collections import defaultdict, Counter
from itertools import combinations_with_replacement
import time


# ========================================================================
# MATCHING GENERATION (same as v1)
# ========================================================================

def gen_matchings(n_vertices):
    vertices = list(range(n_vertices))
    def _gen(verts):
        if len(verts) == 2:
            return [[(verts[0], verts[1])]]
        results = []
        first = verts[0]
        for i, partner in enumerate(verts[1:]):
            remaining = [v for v in verts[1:] if v != partner]
            for sub in _gen(remaining):
                results.append([(first, partner)] + sub)
        return results
    return [frozenset((min(a,b), max(a,b)) for a,b in m) for m in _gen(vertices)]


def double_factorial(n):
    result = 1
    for k in range(1, n+1):
        result *= (2*k - 1)
    return result


def is_prime(n):
    if n < 2: return False
    for d in range(2, int(n**0.5) + 1):
        if n % d == 0: return False
    return True


# ========================================================================
# MULTIPLICATIVE DIRECTION PERMUTATION
# ========================================================================

def compute_dir_perm(p, g=2):
    """Compute the permutation of direction classes induced by ×g mod p.
    
    Direction classes are 0-indexed: class k = edges with diff (k+1).
    Diff d and diff (p-d) are the same direction class.
    
    ×g sends diff d → diff (g*d mod p), identifying d with p-d.
    
    Returns: 
      perm: dict mapping dir class → dir class
      orbits: list of direction class orbits
      cycle_structure: list of orbit lengths
    """
    n_dirs = (p - 1) // 2
    
    def unsigned_dir(d):
        """Map a nonzero residue mod p to its direction class (0-indexed)."""
        d_mod = d % p
        if d_mod == 0:
            return None
        d_unsigned = min(d_mod, p - d_mod)
        return d_unsigned - 1  # 0-indexed
    
    # Compute permutation
    perm = {}
    for k in range(n_dirs):
        diff = k + 1  # direction class k has diff (k+1)
        new_diff = (g * diff) % p
        new_class = unsigned_dir(new_diff)
        if new_class is None:
            return None, [], []  # ×g sends a direction to 0 mod p (composite p)
        perm[k] = new_class
    
    # Check it's actually a permutation (injective on n_dirs classes)
    if len(set(perm.values())) != n_dirs:
        return None, [], []
    
    # Compute orbits
    visited = set()
    orbits = []
    for k in range(n_dirs):
        if k in visited:
            continue
        orbit = []
        current = k
        while current not in visited:
            visited.add(current)
            orbit.append(current)
            current = perm[current]
        orbits.append(orbit)
    
    return perm, orbits, [len(o) for o in orbits]


def find_best_generator(p):
    """Find generator g giving longest orbits on direction classes.
    For prime p, ×g is always a valid permutation (no zero images).
    For composite p, some g values are invalid."""
    n_dirs = (p - 1) // 2
    best_g = 2
    best_max_orbit = 0
    
    for g in range(2, p):
        perm, orbits, lengths = compute_dir_perm(p, g)
        if perm is None:
            continue  # invalid for this g
        max_len = max(lengths) if lengths else 0
        if max_len > best_max_orbit:
            best_max_orbit = max_len
            best_g = g
            if max_len == n_dirs:
                break  # found a full cycle
    
    return best_g, best_max_orbit


# ========================================================================
# DIRECTION STRUCTURE (corrected)
# ========================================================================

class DirectionStructure:
    def __init__(self, n):
        self.n = n
        self.n_vertices = 2 * n
        self.handle_vertex = 2*n - 1
        self.p = 2*n - 1
        self.n_torus_dirs = n - 1
        self.p_is_prime = is_prime(self.p)
        
        # Find the correct multiplicative generator
        self.generator, self.max_orbit_len = find_best_generator(self.p)
        self.dir_perm, self.dir_orbits, self.dir_orbit_lengths = \
            compute_dir_perm(self.p, self.generator)
        
        # Fallback: if no valid perm found, use identity
        if self.dir_perm is None:
            self.dir_perm = {k: k for k in range(self.n_torus_dirs)}
            self.dir_orbits = [[k] for k in range(self.n_torus_dirs)]
            self.dir_orbit_lengths = [1] * self.n_torus_dirs
            self.max_orbit_len = 1
        
        # Full cycle means all direction classes are in one orbit
        self.has_full_cycle = (self.max_orbit_len == self.n_torus_dirs)
        
    def edge_dir(self, i, j):
        a, b = min(i,j), max(i,j)
        if a >= self.p or b >= self.p:
            return -1
        d = (b - a) % self.p
        d = min(d, self.p - d)
        return d - 1
    
    def torus_multiset(self, matching):
        dirs = []
        for a, b in matching:
            d = self.edge_dir(a, b)
            if d >= 0:
                dirs.append(d)
        return tuple(sorted(dirs))
    
    def z3_phase(self, multiset):
        return sum(multiset) % 3
    
    def rotate_multiset(self, ms):
        """Apply the CORRECT direction permutation (×g mod p)."""
        return tuple(sorted(self.dir_perm[d] for d in ms))


# ========================================================================
# ORBIT + EXTRACTOR (same logic as v1 but with corrected permutation)
# ========================================================================

class Orbit:
    def __init__(self, multisets, phase, block_size, orbit_id):
        self.id = orbit_id
        self.multisets = multisets
        self.size = len(multisets)
        self.phase = phase
        self.block_size = block_size
        self.total_matchings = self.size * block_size
        self.full_orbit = multisets
        self.full_orbit_size = len(multisets)
        self.uniform = True
        self.sizes = [block_size] * len(multisets)
        self.empty_count = 0


class K2nExtractor:
    def __init__(self, n, compute_matchings=True):
        self.n = n
        self.n_vertices = 2 * n
        self.n_matchings_expected = double_factorial(n)
        self.dirs = DirectionStructure(n)
        
        self.matchings = None
        self.blocks = None
        self.block_id = None
        self.matching_block = None
        self.orbits = None
        
        if compute_matchings:
            self._compute()
    
    def _compute(self):
        t0 = time.time()
        print(f"K_{self.n_vertices}: generating matchings...", end=" ", flush=True)
        self.matchings = gen_matchings(self.n_vertices)
        N = len(self.matchings)
        print(f"{N} matchings", end="")
        
        self.multisets = [self.dirs.torus_multiset(m) for m in self.matchings]
        self.blocks = defaultdict(list)
        for i, ms in enumerate(self.multisets):
            self.blocks[ms].append(i)
        
        all_ms = sorted(self.blocks.keys())
        self.block_id = {ms: i for i, ms in enumerate(all_ms)}
        self.matching_block = [self.block_id[ms] for ms in self.multisets]
        
        self._compute_orbits()
        t1 = time.time()
        print(f" ({t1-t0:.2f}s)")
    
    def _compute_orbits(self):
        visited = set()
        self.orbits = []
        orbit_id = 0
        
        all_possible = list(combinations_with_replacement(
            range(self.dirs.n_torus_dirs), self.n - 1))
        
        for ms in sorted(all_possible):
            if ms in visited:
                continue
            full_orbit = []
            current = ms
            while current not in visited:
                visited.add(current)
                full_orbit.append(current)
                current = self.dirs.rotate_multiset(current)
            
            nonempty = [m for m in full_orbit if m in self.blocks and len(self.blocks[m]) > 0]
            if not nonempty:
                continue
            
            phase = self.dirs.z3_phase(full_orbit[0])
            sizes = [len(self.blocks[m]) for m in nonempty]
            block_size = sizes[0]
            
            uniform = len(set(sizes)) == 1
            orb = Orbit(nonempty, phase, block_size, orbit_id)
            orb.full_orbit = full_orbit
            orb.full_orbit_size = len(full_orbit)
            orb.uniform = uniform
            orb.sizes = sizes
            orb.empty_count = len(full_orbit) - len(nonempty)
            if not uniform:
                orb.block_size = max(sizes)
            
            self.orbits.append(orb)
            orbit_id += 1
    
    def check_zp_preserves_blocks(self):
        p = self.dirs.p
        N = len(self.matchings)
        m_to_idx = {m: i for i, m in enumerate(self.matchings)}
        
        def apply_zp(m):
            return frozenset(
                (min((a+1)%p if a < p else p, (b+1)%p if b < p else p),
                 max((a+1)%p if a < p else p, (b+1)%p if b < p else p))
                for a, b in m)
        
        for i in range(N):
            j = m_to_idx[apply_zp(self.matchings[i])]
            if self.matching_block[i] != self.matching_block[j]:
                return False
        return True
    
    def check_mult_preserves_blocks(self):
        """Check if ×g mod p (as vertex map) preserves direction blocks."""
        p = self.dirs.p
        g = self.dirs.generator
        N = len(self.matchings)
        m_to_idx = {m: i for i, m in enumerate(self.matchings)}
        
        def apply_mult(m):
            """Apply i → g*i mod p on torus vertices, fix handle vertex."""
            new_edges = set()
            for a, b in m:
                na = (g * a) % p if a < p else p
                nb = (g * b) % p if b < p else p
                new_edges.add((min(na, nb), max(na, nb)))
            return frozenset(new_edges)
        
        # Check that every image matching exists
        for i in range(N):
            img = apply_mult(self.matchings[i])
            if img not in m_to_idx:
                return False, "image matching not in set"
        
        # Check block preservation in the orbit sense:
        # ×g should map multiset ms → rotate_multiset(ms)
        mismatches = 0
        for i in range(N):
            img = apply_mult(self.matchings[i])
            j = m_to_idx[img]
            ms_i = self.multisets[i]
            ms_j = self.multisets[j]
            expected = self.dirs.rotate_multiset(ms_i)
            if ms_j != expected:
                mismatches += 1
        
        return mismatches == 0, f"{mismatches} mismatches out of {N}"
    
    def report(self):
        p = self.dirs.p
        nd = self.dirs.n_torus_dirs
        N = len(self.matchings)
        g = self.dirs.generator
        
        print(f"\n{'='*72}")
        print(f"K_{self.n_vertices} ORBIT DECOMPOSITION (×{g} mod {p})")
        print(f"{'='*72}")
        
        print(f"\n  Vertices: {self.n_vertices}, p = {p} ({'PRIME' if self.dirs.p_is_prime else 'COMPOSITE'})")
        print(f"  Matchings: {N}")
        print(f"  Torus directions: {nd}")
        print(f"  Generator: ×{g} mod {p}")
        
        # Direction permutation
        perm_str = ", ".join(f"{k}→{v}" for k,v in sorted(self.dirs.dir_perm.items()))
        print(f"  Dir perm: {perm_str}")
        print(f"  Dir orbits: {self.dirs.dir_orbits}")
        print(f"  Full cycle: {self.dirs.has_full_cycle}")
        
        # Verify multiplicative map preserves blocks
        mult_ok, mult_msg = self.check_mult_preserves_blocks()
        print(f"  ×{g} preserves block structure: {mult_ok} ({mult_msg})")
        
        # Blocks
        n_blocks = len(self.blocks)
        print(f"\n  Direction blocks: {n_blocks}")
        
        # Orbits summary
        n_uniform = sum(1 for o in self.orbits if o.uniform)
        n_broken = sum(1 for o in self.orbits if not o.uniform)
        print(f"\n  Orbits: {len(self.orbits)} total ({n_uniform} uniform, {n_broken} broken)")
        
        for orb in self.orbits:
            fixed = " (fixed pt)" if orb.full_orbit_size == 1 else ""
            broken = " ⚠ BROKEN" if not orb.uniform else ""
            ms_display = [str(m) for m in orb.multisets[:3]]
            if len(orb.multisets) > 3:
                ms_display.append("...")
            ms_str = ", ".join(ms_display)
            
            if orb.uniform:
                print(f"    Orbit {orb.id:>2}: [{ms_str}]{fixed}")
                print(f"             |orbit|={orb.size}, block_size={orb.block_size}, "
                      f"total={orb.total_matchings}")
                if orb.size > 1:
                    print(f"             → Gal(deg-{orb.block_size}) ≀ C_{orb.size}")
            else:
                print(f"    Orbit {orb.id:>2}: [{ms_str}]{broken}")
                print(f"             sizes={orb.sizes}, {orb.empty_count} empty blocks")
        
        # Summary
        print(f"\n  ORBIT SIZE DISTRIBUTION:")
        orbit_sizes = Counter(orb.size for orb in self.orbits if orb.uniform)
        broken_sizes = Counter(orb.size for orb in self.orbits if not orb.uniform)
        for s, c in sorted(orbit_sizes.items()):
            print(f"    {c} uniform orbit(s) of size {s}")
        for s, c in sorted(broken_sizes.items()):
            print(f"    {c} BROKEN orbit(s) of size {s}")
        
        return self


# ========================================================================
# RUN
# ========================================================================

if __name__ == "__main__":
    
    # First: show the direction permutations
    print("=" * 72)
    print("DIRECTION PERMUTATIONS (×2 mod p)")
    print("=" * 72)
    for n in [3, 4, 5, 6, 7]:
        p = 2*n - 1
        g, max_orb = find_best_generator(p)
        perm, orbits, lengths = compute_dir_perm(p, g)
        full = "FULL CYCLE" if max_orb == n-1 else f"max orbit {max_orb}/{n-1}"
        perm_str = ", ".join(f"{k}→{v}" for k,v in sorted(perm.items()))
        print(f"  K_{2*n} (p={p:>2}): ×{g}, [{perm_str}], orbits {orbits} — {full}")
    
    print()
    
    # K₆ verify
    e6 = K2nExtractor(3)
    e6.report()
    
    # K₈ verify
    e8 = K2nExtractor(4)
    e8.report()
    
    # K₁₀ (composite p=9)
    e10 = K2nExtractor(5)
    e10.report()
    
    # K₁₂ (prime p=11)
    e12 = K2nExtractor(6)
    e12.report()
    
    # K₁₄ (prime p=13)
    if double_factorial(7) < 150000:
        print(f"\n\nK_14 has {double_factorial(7)} matchings")
        e14 = K2nExtractor(7)
        e14.report()
        all_ext = [e6, e8, e10, e12, e14]
    else:
        print(f"\n\nK_14 has {double_factorial(7)} matchings — SKIPPING (too large)")
        e14 = None
        all_ext = [e6, e8, e10, e12]
    
    # Cross-level
    print(f"\n{'='*72}")
    print("CROSS-LEVEL SUMMARY")
    print(f"{'='*72}")
    print(f"\n  {'Level':<6} {'p':<4} {'prime?':<8} {'gen':<4} {'full?':<7} "
          f"{'blocks':<8} {'orbits':<8} {'uniform':<8} {'broken':<8}")
    print(f"  {'-'*65}")
    for ext in all_ext:
        n_u = sum(1 for o in ext.orbits if o.uniform)
        n_b = sum(1 for o in ext.orbits if not o.uniform)
        fc = "YES" if ext.dirs.has_full_cycle else "NO"
        pr = "YES" if ext.dirs.p_is_prime else "NO"
        print(f"  K_{ext.n_vertices:<4} {ext.dirs.p:<4} {pr:<8} ×{ext.dirs.generator:<3} {fc:<7} "
              f"{len(ext.blocks):<8} {len(ext.orbits):<8} {n_u:<8} {n_b:<8}")