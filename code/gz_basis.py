"""
Gelfand-Tsetlin Basis and Bordered Eigenproblem
================================================
Constructs the GZ basis for V_{(2n-2,2)} via SYT enumeration,
implements Young's seminormal form, and solves bordered eigenproblems.
"""
import numpy as np
from scipy.linalg import eigh
from itertools import combinations
import time

def standard_young_tableaux(m, row2_len=2):
    """All SYT of shape (m-row2_len, row2_len)."""
    entries = list(range(1, m+1))
    tableaux = []
    for combo in combinations(entries, row2_len):
        row2 = list(combo)
        row1 = [x for x in entries if x not in combo]
        valid = all(row2[j] > row1[j] for j in range(row2_len))
        if valid:
            tableaux.append((tuple(row1), tuple(row2)))
    return tableaux

def content(tableau, val):
    """Content of entry val in tableau: col - row (0-indexed)."""
    for c, v in enumerate(tableau[0]):
        if v == val: return c - 0  # row 0
    for c, v in enumerate(tableau[1]):
        if v == val: return c - 1  # row 1
    raise ValueError(f"{val} not in tableau")

def axial_distance(tableau, i, j):
    return content(tableau, j) - content(tableau, i)

def swap_entries(tableau, a, b):
    """Swap a,b in tableau; return new SYT or None if invalid."""
    rows = [list(tableau[0]), list(tableau[1])]
    pa = pb = None
    for r in range(2):
        for c in range(len(rows[r])):
            if rows[r][c] == a: pa = (r,c)
            if rows[r][c] == b: pb = (r,c)
    if not pa or not pb: return None
    rows[pa[0]][pa[1]], rows[pb[0]][pb[1]] = b, a
    # Check standard
    for r in range(2):
        for c in range(len(rows[r])-1):
            if rows[r][c] >= rows[r][c+1]: return None
    for c in range(min(len(rows[0]), len(rows[1]))):
        if rows[0][c] >= rows[1][c]: return None
    return (tuple(rows[0]), tuple(rows[1]))

def rep_matrix_transposition(tableaux, k):
    """Matrix of adjacent transposition (k,k+1) in seminormal form."""
    d = len(tableaux)
    idx = {t:i for i,t in enumerate(tableaux)}
    M = np.zeros((d,d))
    for i, T in enumerate(tableaux):
        a = axial_distance(T, k, k+1)
        M[i,i] = 1.0/a
        Ts = swap_entries(T, k, k+1)
        if Ts and Ts in idx:
            M[i, idx[Ts]] = np.sqrt(1.0 - 1.0/a**2)
    return M

def verify_seminormal(m):
    """Verify seminormal form satisfies S_m relations."""
    tabs = standard_young_tableaux(m)
    d = len(tabs)
    mats = {k: rep_matrix_transposition(tabs, k) for k in range(1, m)}
    
    ok = True
    # T_kÂ² = I
    for k in range(1, m):
        if np.max(np.abs(mats[k] @ mats[k] - np.eye(d))) > 1e-10:
            ok = False; print(f"  T_{k}Â² â‰  I!")
    # Braid: T_k T_{k+1} T_k = T_{k+1} T_k T_{k+1}
    for k in range(1, m-1):
        lhs = mats[k] @ mats[k+1] @ mats[k]
        rhs = mats[k+1] @ mats[k] @ mats[k+1]
        if np.max(np.abs(lhs - rhs)) > 1e-10:
            ok = False; print(f"  Braid {k},{k+1} fails!")
    # Far commutativity
    for j in range(1, m):
        for k in range(j+2, m):
            if np.max(np.abs(mats[j]@mats[k] - mats[k]@mats[j])) > 1e-10:
                ok = False; print(f"  [{j},{k}] â‰  0!")
    return ok

def branching_structure(n):
    """
    Identify the 4-piece branching of V_{(2n-2,2)} under S_{2n} â†’ S_{2n-2}.
    Returns indices of SYT belonging to each piece:
    - inherited: SYT of shape (2n-4,2) (entries 1..2n-2 in same positions)
    - standard_1, standard_2: two copies of (2n-3,1)
    - trivial: the single (2n-2) direction
    """
    m = 2*n
    tabs = standard_young_tableaux(m)
    
    # A SYT of shape (2n-2,2) restricted to entries {1,...,2n-2} gives
    # a skew tableau. The branching corresponds to where entries 2n-1 and 2n sit.
    
    inherited = []  # entries 2n-1, 2n both in row 1
    std_from_row1 = []  # 2n-1 in row 1, 2n in row 2 (or path through (2n-3,2)â†’(2n-4,2))
    std_from_row2 = []  # 2n-1 in row 2, 2n in row 1 (path through (2n-2,1)â†’(2n-3,1))
    trivial = []  # This is harder to identify from positions alone
    
    for i, T in enumerate(tabs):
        row1, row2 = T
        pos_2nm1 = 'r1' if (m-1) in row1 else 'r2'
        pos_2n = 'r1' if m in row1 else 'r2'
        
        if pos_2nm1 == 'r1' and pos_2n == 'r1':
            inherited.append(i)
        elif pos_2nm1 == 'r1' and pos_2n == 'r2':
            std_from_row1.append(i)
        elif pos_2nm1 == 'r2' and pos_2n == 'r1':
            std_from_row2.append(i)
        else:  # both in row 2
            trivial.append(i)
    
    return {
        'inherited': inherited,
        'std_1': std_from_row1,
        'std_2': std_from_row2,
        'trivial': trivial,
        'all_tabs': tabs,
    }

if __name__ == "__main__":
    print("GZ BASIS VERIFICATION")
    print("=" * 60)
    
    for n in [3, 4, 5]:
        m = 2*n
        tabs = standard_young_tableaux(m)
        d = len(tabs)
        d_expected = n*(2*n-3)
        print(f"\nS_{m}, shape ({m-2},2): {d} SYT (expected {d_expected}) {'âœ“' if d==d_expected else 'âœ—'}")
        
        ok = verify_seminormal(m)
        print(f"  Seminormal relations: {'âœ“' if ok else 'âœ—'}")
        
        br = branching_structure(n)
        ni = len(br['inherited'])
        ns1 = len(br['std_1'])
        ns2 = len(br['std_2'])
        nt = len(br['trivial'])
        d_prev = (n-1)*(2*n-5) if n >= 3 else 0
        print(f"  Branching: {ni} inherited + {ns1}+{ns2} std + {nt} triv = {ni+ns1+ns2+nt}")
        print(f"    Expected: {d_prev} + 2Ã—{2*n-3} + 1 = {d_prev + 2*(2*n-3) + 1}")
