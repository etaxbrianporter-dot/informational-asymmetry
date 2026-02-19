---
title: The Spectral Budget Identity
section: Mathematical Structure
status: active
---

# The Spectral Budget Identity

## The Discovery

The hub memory fraction in V_physical equals the trace spectral weight in V_physical:

$$\frac{\|\Pi_{\text{phys}}(h)\|^2}{\|h\|^2} = \frac{\lambda_{\text{mid}} \cdot d_{\text{phys}}}{\text{Tr}(O)} = \frac{2n-2}{2n-1}$$

This is not a coincidence. It is forced by a single algebraic structure.

## The Incidence Matrix

The overlap matrix factors as **O = 2AA^T**, where A is the N Ã— E matching-edge incidence matrix:

- A_{i,e} = 1 if matching i contains edge e, else 0
- N = (2n-1)!! matchings, E = n(2n-1) edges

**A is doubly regular:**
- Row sums: r = n (each matching has n edges)
- Column sums: c = (2n-3)!! (each edge appears in (2n-3)!! matchings = N_hub)

Both regularities follow from S_{2n} acting transitively on matchings and edges respectively.

## Why Hub = Trace

From double regularity alone:

**Trace side:** Î»_max = 2rc (row sum of O), Tr(O) = 2Nr (Frobenius norm of A), so Î»_max/Tr(O) = c/N.

**Hub side:** h = Ae_{eâ‚€} (column of A), so ||h||Â² = c and hÂ·ðŸ = c. Projection: ||Î _triv(h)||Â²/||h||Â² = cÂ²/(Nc) = c/N.

**Kernel side:** ker(O) = ker(A^T), and h is a column of A, so h âŠ¥ ker(O) (Part A).

**Therefore:** Physical fraction = 1 - c/N for both. QED.

The two fractions are the same because both reduce to the **aspect ratio c/N** of the doubly-regular incidence matrix. This ratio is the single number 1/(2n-1) that controls everything.

## The Dual: Edge Co-occurrence and the Kneser Graph

The dual matrix A^T A is the E Ã— E edge co-occurrence matrix:

$$(A^T A)_{e,e'} = |\{\text{matchings containing both } e \text{ and } e'\}|$$

This has a clean structure:
- Diagonal: (2n-3)!! (= N_hub, matchings per edge)
- Adjacent edges (share vertex): 0 (a matching has disjoint edges only)
- Disjoint edges: (2n-5)!! (matchings of K_{2n-4} on remaining vertices)

So: **A^T A = (2n-3)!! Â· I + (2n-5)!! Â· Adj(KG(2n,2))**

where KG(2n,2) is the Kneser graph on edges of K_{2n} (adjacent iff disjoint).

## The Johnson Scheme Forces Three Eigenvalues

The Kneser graph KG(2n,2) lives in the Johnson association scheme J(2n,2), which has exactly **3 classes**:
- Class 0: identity (same edge)
- Class 1: edges sharing one vertex (Johnson graph adjacency)
- Class 2: disjoint edges (Kneser adjacency)

Three classes â†’ three simultaneous eigenspaces â†’ A^T A has three eigenvalues â†’ O = 2AA^T has three eigenvalues.

**This is why the overlap matrix has exactly three eigenvalues.** It's not a numerical accident â€” it's forced by the Johnson scheme.

## The Universal Kneser Eigenvalue

The Kneser graph eigenvalues on each eigenspace of J(2n,2):

| Eigenspace | Dimension | Kneser eigenvalue | Identity |
|:---|:---|:---|:---|
| Vâ‚€ (trivial) | 1 | C(2n-2, 2) | (2n-2)(2n-3)/2 |
| Vâ‚ (standard) | 2n-1 | -(2n-3) | kills A^T A eigenvalue â†’ kernel |
| Vâ‚‚ (symmetric traceless) | n(2n-3) | **1** | universal, independent of n |

**The Kneser eigenvalue on Vâ‚‚ is always 1.**

Proof: The Johnson adjacency J and Kneser adjacency KG satisfy J + KG = J_{\text{allones}} - I. On Vâ‚‚: J has eigenvalue -2, J_{\text{allones}} has eigenvalue 0, I has eigenvalue 1. So KG eigenvalue = 0 - 1 - (-2) = **1**. â–¡

## Closed Forms for Everything

From A^T A = (2n-3)!! I + (2n-5)!! Â· Adj(KG):

**Eigenvalues of A^T A:**
- Vâ‚€: (2n-3)!! + (2n-5)!! Â· C(2n-2,2) â†’ Î»_max(O)/2
- Vâ‚: (2n-3)!! + (2n-5)!! Â· (-(2n-3)) = (2n-3)!! - (2n-3)!! = **0** â†’ kernel
- Vâ‚‚: (2n-3)!! + (2n-5)!! Â· 1 = (2n-3)!! + (2n-5)!! â†’ Î»_mid(O)/2

**Therefore:**
- Î»_max = 2n(2n-3)!!
- **Î»_mid = 2[(2n-3)!! + (2n-5)!!] = 4(n-1)(2n-5)!!**
- d_phys = n(2n-3) (dimension of Vâ‚‚ in Johnson scheme = hook length for S_{2n} irrep (2n-2,2))

**The identity Î»_mid Ã— d_phys = 2(n-1) Ã— Î»_max** is the budget constraint: physical spectral weight = (2n-2) trivial weights.

**The ratio Î»_mid/Î»_max = 2(n-1)/[n(2n-3)]** â†’ 0 as n â†’ âˆž, but the total physical weight (2n-2)/(2n-1) â†’ 1.

## The Physical Eigenvalue Decomposition

Î»_mid/2 = (2n-3)!! + (2n-5)!! Ã— **1**

This decomposes as:
- **(2n-3)!!** = matchings per edge (self-correlation contribution)
- **(2n-5)!!** = co-occurrences of disjoint edge pairs Ã— Kneser eigenvalue **1**

The Kneser eigenvalue 1 means: in the physical eigenspace, each edge correlates with each disjoint edge by exactly the minimum unit. The Johnson scheme forces this universally.

## Self-Focusing of the Chain

The physical eigenspace is a **vanishing** fraction of total matching space:

| K_{2n} | N = (2n-1)!! | d_phys = n(2n-3) | Fraction | Hub in phys |
|:---|:---|:---|:---|:---|
| Kâ‚† | 15 | 9 | 60% | 4/5 = 80% |
| Kâ‚ˆ | 105 | 20 | 19% | 6/7 = 86% |
| Kâ‚â‚€ | 945 | 35 | 3.7% | 8/9 = 89% |
| Kâ‚â‚‚ | 10395 | 54 | 0.5% | 10/11 = 91% |
| Kâ‚â‚† | 2027025 | 104 | 0.005% | 14/15 = 93% |

The physical sector shrinks super-exponentially as a fraction of the full space, yet captures a growing fraction of the hub memory. The chain **self-focuses**: each step concentrates the previous level's information more precisely into the degrees of freedom that matter.

The budget identity Î»_max = 2n Ã— N_hub is the **lens equation** governing this focusing. It ensures that despite the kernel dominating the dimension count, the hub aims almost entirely into V_physical.

## Cross-Block Structure (Part E)

The hub/non-hub partition O = [[H, B], [B^T, C]] yields:

**BB^T has two nonzero eigenvalues:**
- ÏƒÂ²_max = bÂ² N_hub/N_non (multiplicity 1, trivial coupling)
- ÏƒÂ²_mid = Î»_mid^{prev} Ã— Î”Î»_mid (multiplicity d_phys^{prev})

where Î”Î»_mid = Î»_mid(K_{2n}) - Î»_mid(K_{2n-2}) = 4(2n-5)!!(2nÂ²-4n+1).

**The cross-block bandwidth = d_phys of the previous level.** The coupling between hub and non-hub sectors has exactly as many independent channels as there are physical degrees of freedom at the previous level.

## The Complete Deductive Chain

```
S_{2n} acts on K_{2n}
    â”‚
    â”œâ”€â†’ transitive on edges â”€â”€â†’ A doubly regular (c = (2n-3)!!)
    â”‚                              â”‚
    â”œâ”€â†’ transitive on matchings â”€â†’ A doubly regular (r = n)
    â”‚                              â”‚
    â”‚                              â”œâ”€â†’ O = 2AA^T
    â”‚                              â”‚     â”œâ”€â†’ Î»_max = 2rc = 2n(2n-3)!!
    â”‚                              â”‚     â”œâ”€â†’ Tr(O) = 2Nr = 2n(2n-1)!!
    â”‚                              â”‚     â””â”€â†’ ker(O) = ker(A^T)
    â”‚                              â”‚
    â”‚                              â””â”€â†’ h âŠ¥ ker(O)           [Part A]
    â”‚                                   h fractions = c/N    [Part B]
    â”‚
    â”œâ”€â†’ transitive on edges â”€â”€â†’ P_e all conjugate
    â”‚                              â””â”€â†’ Tr(P_e|V_Î»)/dim = c/N [Part C]
    â”‚
    â””â”€â†’ Johnson scheme J(2n,2) has 3 classes
           â”‚
           â”œâ”€â†’ Kneser eigenvalues: {C(2n-2,2), -(2n-3), 1}
           â”‚     â”œâ”€â†’ THREE eigenvalues for O
           â”‚     â”œâ”€â†’ Î»_mid = 4(n-1)(2n-5)!!  [from universal 1]
           â”‚     â””â”€â†’ d_phys = n(2n-3)
           â”‚
           â””â”€â†’ Hub embedding: O_hub = O_prev + 2J [Part D]
                 â””â”€â†’ Cross-block mult = d_phys^prev [Part E]
                      â””â”€â†’ ÏƒÂ²_mid = Î»_mid^prev Ã— Î”Î»_mid
```

**Everything** â€” the three eigenvalues, their values, the hub fractions, the isotropy, the eigenvalue inheritance, the cross-block multiplicity, and the budget identity â€” follows from:

1. O = 2AA^T (overlap factors through edge incidence)
2. A is doubly regular (S_{2n} transitive on edges and matchings)
3. J(2n,2) has 3 classes (Johnson association scheme)
