---
title: Why the Vacuum Selects the Doublet
section: 
status: active
---

# Why the Vacuum Selects the Doublet

## The Phase Scan Result

For the sorted vacuum's direction group (matchings 0,1,2,4,9), all 81 Z₃ phase assignments give:

| Phase assignment | a₂ | a₄/a₂ | mH (GeV) |
|---|---|---|---|
| [0,0,0,0,0] (min a₂) | 2.877 | 1.498 | 139.1 |
| [0,1,2,1,2] (sorted) | 3.306 | 1.231 | **126.1** |
| [0,1,0,0,0] (max ratio) | 4.000 | 2.000 | 160.8 |
| Range across all 81 | [2.88, 4.00] | [1.17, 2.00] | [122.8, 160.8] |

Compare with the non-doublet champion (matchings 0,1,4,13,14):

| | a₂ | a₄/a₂ | mH (GeV) |
|---|---|---|---|
| Best phases | 1.101 | 0.232 | **54.7** |

## The Mechanism

**The doublet group maintains a₄/a₂ > 1.17 for ALL phase choices.** The non-doublet champion crushes this ratio to 0.23. The ratio, not the individual eigenvalues, determines physics.

The sorted vacuum's phases [0,1,2,1,2] don't minimize a₂ (that's 2.877). They don't maximize a₄/a₂ (that's 2.0). They produce the **specific ratio** a₄/a₂ = 1.2305 that gives mH = 126.1 GeV.

## What This Means for the "4"

The formula mH² = 2(a₄/a₂)mW² works because:

1. **K₆ combinatorics** forces a₄/a₂ into a bounded range for doublet configurations
2. **The doublet structure** (hub edge shared by 3 matchings) is required to maintain quartic stiffness relative to quadratic — non-doublet configurations are too "floppy"
3. **Z₃ phases** tune within the doublet range, and the sorted assignment lands at 1.2305

The factor of 4 in λ/g² = a₄/(4a₂) is not a separate structural constant to be derived. It's absorbed into the ratio a₄/a₂ = 1.2305, which IS the physical content. The "4" appears when you FACTOR this ratio through the SM relations mH² = 2λv², v = 2mW/g.

**The 4 is the SM's 4, not the graph's 4.** It comes from:
- mH² = 2λv² (Higgs potential shape, factor 2)  
- v = 2mW/g (SU(2) doublet, factor 2)
- Combined: mH² = 8(λ/g²)mW²

Writing λ/g² = a₄/(4a₂) just means 8/4 = 2, i.e., mH² = 2(a₄/a₂)mW².

## The Real Achievement

K₆ doesn't explain the 4. K₆ explains the **1.2305**. This pure combinatorial number, times 2 (from SM), times mW², gives the Higgs mass to 0.8%.

The doublet selection is forced by needing a₄/a₂ in the right range. The specific value 1.2305 comes from the interplay of:
- 120:90:15 overlap distribution (K₆ graph invariant)
- 5-matching direction group with doublet hub structure
- Z₃ phase pattern [0,1,2,1,2]
- Ground eigenvector of the resulting 5×5 Gram matrix
