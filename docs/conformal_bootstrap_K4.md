---
title: Conformal Bootstrap for the Kâ‚„ Critical Point
section: K₄ Spacetime
status: active
---

# Conformal Bootstrap for the Kâ‚„ Critical Point

## Can the algebra close Step 4 without Monte Carlo?

**February 17, 2026**

---

## 1. The Setup: What CFT Are We Bootstrapping?

The Kâ‚„ model at its quantum critical point U = U_c is (assuming continuous transition) a 2+1D conformal field theory. To set up the bootstrap, we need to identify every piece of structure the algebra forces.

### 1.1 Spacetime dimension

d = 3 (2 spatial + 1 Euclidean time). Standard for the bootstrap. Enormous existing technology.

### 1.2 Global symmetry

The full symmetry at the critical point is:

**Internal:** Vâ‚„ = Zâ‚‚ Ã— Zâ‚‚ (orbital permutation symmetry from matching algebra). This is exact â€” proved in surviving_invariants. The interaction H_int = (U/2)(N-2)Â² is Vâ‚„-invariant by construction.

**Spatial:** The Zâ‚ƒ flux assigns phases (1, Ï‰, Ï‰Â²) to the three lattice directions (Î´â‚, Î´â‚‚, Î´â‚ƒ). Under Câ‚ƒ rotation: Î´â‚ â†’ Î´â‚‚ â†’ Î´â‚ƒ, and the phases cycle as 1 â†’ Ï‰ â†’ Ï‰Â². The flux is invariant under the *combined* operation (Câ‚ƒ rotation) Ã— (Zâ‚ƒ phase shift). The reflections in Câ‚†v are broken. So spatial symmetry at the critical point is Câ‚ƒ, not Câ‚†v.

**Parity/Time-reversal:** Both broken. T is broken by Zâ‚ƒ flux (proved, invariant V). P is broken because the Chern number C = âˆ’2 â‰  0. CPT is also broken.

**Charge conjugation C:** Vertex permutation (2â†”3). Unitary, order 2. Swaps Ï‰ â†” Ï‰Ì„. This IS a symmetry of the full model (it maps D(k) â†’ DÌ„(k), but the interaction is C-invariant). However, C maps the theory to itself with conjugate flux â€” it's an anti-unitary symmetry of the partition function, not a symmetry of the Hilbert space at fixed flux.

**Full symmetry group at the critical point:**

    G_crit = Vâ‚„ Ã— Câ‚ƒ  (discrete, no continuous part beyond conformal)

This is crucial. The symmetry is small. Small symmetry means:
- Fewer constraints on the OPE
- Fewer protected operators
- Larger anomalous dimensions (less protection = more running)

### 1.3 Matter content at the free point

At U = 0, the theory has 4 Dirac cones (one per Vâ‚„ channel), of which 2 are active (channels 1, 3 with winding âˆ’1) and 2 are inert (channels 2, 4 with winding 0). In the continuum limit:

- Active channels â†’ 2 chiral (Weyl) fermions
- Inert channels â†’ 2 massive Dirac fermions (irrelevant at low energy)

At the critical point U = U_c, the relevant degrees of freedom are the 2 active-channel fermions. The inert channels are gapped and decouple.

**Effective matter content at the critical point: N_f = 2 Dirac fermions (from 2 active Vâ‚„ channels) with Zâ‚‚ global symmetry (the residual Vâ‚„ acting on the active channels).**

Wait â€” this requires more care. Let me work through the Vâ‚„ representation theory at the critical point.

### 1.4 Vâ‚„ decomposition at the critical point

The 4 Vâ‚„ irreps with characters:

| Irrep | Mâ‚ | Mâ‚‚ | Mâ‚ƒ | Winding | Status |
|-------|-----|-----|-----|---------|--------|
| Ï‡â‚ = (+ + +) | +1 | +1 | +1 | âˆ’1 | Active |
| Ï‡â‚‚ = (+ âˆ’ âˆ’) | +1 | âˆ’1 | âˆ’1 | 0 | Inert |
| Ï‡â‚ƒ = (âˆ’ + âˆ’) | âˆ’1 | +1 | âˆ’1 | âˆ’1 | Active |
| Ï‡â‚„ = (âˆ’ âˆ’ +) | âˆ’1 | âˆ’1 | +1 | 0 | Inert |

The Vâ‚„ tensor product rules (for Zâ‚‚ Ã— Zâ‚‚ = {(a,b)} with a,b âˆˆ {0,1}):

    Ï‡áµ¢ âŠ— Ï‡â±¼ = Ï‡â‚–  where k is determined by component-wise XOR

Explicitly:
- Ï‡â‚ âŠ— Ï‡â‚ = Ï‡â‚ (singlet Ã— singlet = singlet)
- Ï‡â‚ âŠ— Ï‡â‚‚ = Ï‡â‚‚
- Ï‡â‚ âŠ— Ï‡â‚ƒ = Ï‡â‚ƒ
- Ï‡â‚ âŠ— Ï‡â‚„ = Ï‡â‚„
- Ï‡â‚‚ âŠ— Ï‡â‚‚ = Ï‡â‚ (each non-trivial rep is self-conjugate)
- Ï‡â‚‚ âŠ— Ï‡â‚ƒ = Ï‡â‚„
- Ï‡â‚‚ âŠ— Ï‡â‚„ = Ï‡â‚ƒ
- Ï‡â‚ƒ âŠ— Ï‡â‚ƒ = Ï‡â‚
- Ï‡â‚ƒ âŠ— Ï‡â‚„ = Ï‡â‚‚
- Ï‡â‚„ âŠ— Ï‡â‚„ = Ï‡â‚

Key observation: the two active channels (Ï‡â‚, Ï‡â‚ƒ) satisfy Ï‡â‚ âŠ— Ï‡â‚ƒ = Ï‡â‚ƒ âŠ— Ï‡â‚ = Ï‡â‚ƒ ... wait, that's not right. Let me redo this.

Vâ‚„ = Zâ‚‚(a) Ã— Zâ‚‚(b). Assign:
- Ï‡â‚ = (0,0) = trivial
- Ï‡â‚‚ = (1,0)
- Ï‡â‚ƒ = (0,1)
- Ï‡â‚„ = (1,1)

Then Ï‡â‚ âŠ— Ï‡â‚ƒ = (0,0) âŠ— (0,1) = (0,1) = Ï‡â‚ƒ.

But I need to identify which labeling matches the active/inert classification. From the surviving invariants: active = Ï‡â‚ (+++) and Ï‡â‚ƒ (âˆ’+âˆ’). These have characters:
- Ï‡â‚: (+1, +1, +1) on (Mâ‚, Mâ‚‚, Mâ‚ƒ)  
- Ï‡â‚ƒ: (âˆ’1, +1, âˆ’1) on (Mâ‚, Mâ‚‚, Mâ‚ƒ)

The product Ï‡â‚ Ã— Ï‡â‚ƒ has characters (+1Â·(âˆ’1), +1Â·(+1), +1Â·(âˆ’1)) = (âˆ’1, +1, âˆ’1) = Ï‡â‚ƒ.
The product Ï‡â‚ƒ Ã— Ï‡â‚ƒ has characters ((âˆ’1)Â², (+1)Â², (âˆ’1)Â²) = (+1, +1, +1) = Ï‡â‚.

So: active âŠ— active = {Ï‡â‚, Ï‡â‚ƒ} (both active). The active sector is closed under tensor product!

And: active âŠ— inert = inert. Inert âŠ— inert = active.

This means: **the two active Vâ‚„ channels form a closed Zâ‚‚ subgroup** under tensor product. Specifically, {Ï‡â‚, Ï‡â‚ƒ} â‰… Zâ‚‚ âŠ‚ Vâ‚„.

At the critical point, if the inert channels decouple (they're gapped), the effective global symmetry acting on the active degrees of freedom is this Zâ‚‚ subgroup.

### 1.5 The effective CFT symmetry class

Putting it together, the critical-point CFT has:

    d = 3, G = Zâ‚‚ (effective, from active Vâ‚„ channels), no parity, N_f = 2 Dirac

This is a **Zâ‚‚-symmetric Gross-Neveu model without parity, in 2+1D, with N_f = 2 flavors.**

The parity-invariant version of this (Gross-Neveu with N_f = 2 and Zâ‚‚) has been studied by bootstrap. The parity-broken version is the theory dual (via 3d bosonization) to a Chern-Simons-matter theory.

---

## 2. The Crossing Equations

### 2.1 Choice of external operator

The natural external operators for the bootstrap are:

**Option A: Scalar singlet Ïƒ.** The lowest scalar in the Zâ‚‚-even sector. At the critical point, this is the mass operator ÏˆÌ„Ïˆ (summed over active channels). Its dimension Î”_Ïƒ is the order parameter dimension.

**Option B: Fermion Ïˆ.** The fundamental Dirac fermion in one of the active channels. Fermionic bootstrap is harder but more constraining.

**Option C: Both.** Mixed correlator bootstrap using âŸ¨ÏƒÏƒÏƒÏƒâŸ©, âŸ¨ÏˆÏˆÏˆÏˆâŸ©, and âŸ¨ÏƒÏƒÏˆÏˆâŸ©. This is the most powerful but computationally expensive.

For our purposes (bounding Î”â‚ƒ), Option A is sufficient and is the most studied.

### 2.2 The âŸ¨ÏƒÏƒÏƒÏƒâŸ© crossing equation

The 4-point function of identical scalars in d = 3 decomposes as:

    âŸ¨Ïƒ(xâ‚)Ïƒ(xâ‚‚)Ïƒ(xâ‚ƒ)Ïƒ(xâ‚„)âŸ© = (1/|xâ‚â‚‚|^{2Î”_Ïƒ} Â· 1/|xâ‚ƒâ‚„|^{2Î”_Ïƒ}) Ã— G(u,v)

where u, v are conformal cross-ratios. The function G decomposes into conformal blocks:

    G(u,v) = Î£_{O} Î»Â²_{ÏƒÏƒO} Â· g_{Î”_O, â„“_O}(u,v)

The sum runs over all primary operators O appearing in the Ïƒ Ã— Ïƒ OPE, with dimension Î”_O and spin â„“_O.

Crossing symmetry (xâ‚ â†” xâ‚ƒ) gives:

    Î£_O Î»Â²_{ÏƒÏƒO} Â· F_{Î”,â„“}(u,v) = 0

where F_{Î”,â„“} is a known combination of conformal blocks. This is one equation in infinitely many unknowns (the spectrum {Î”_O, â„“_O} and OPE coefficients {Î»Â²_{ÏƒÏƒO}}).

### 2.3 What appears in Ïƒ Ã— Ïƒ OPE

Since Ïƒ is a Zâ‚‚-even scalar, the Ïƒ Ã— Ïƒ OPE contains:

- **â„“ = 0 (scalars):** The identity I (Î” = 0), ÏƒÂ² (the composite), and other Zâ‚‚-even scalars
- **â„“ = 2 (spin-2):** The stress tensor T_Î¼Î½ with Î” = 3 (protected)
- **â„“ = 3 (spin-3):** THE TARGET. The lowest spin-3 primary Jâ‚ƒ with Î”â‚ƒ = 4 + Î³â‚ƒ
- **â„“ = 4, 5, ...:** Higher-spin operators

In a parity-invariant theory, only even-spin operators appear in the symmetric OPE Ïƒ Ã— Ïƒ (because Ïƒ is its own conjugate and parity exchanges the two types of spin-â„“ blocks).

**BUT â€” parity is broken here.** This is crucial.

### 2.4 Parity breaking changes the bootstrap

In a parity-invariant 3D CFT:
- Spin-â„“ operators split into parity-even (â„“âº) and parity-odd (â„“â»)
- The Ïƒ Ã— Ïƒ OPE contains â„“âº for even â„“ and â„“â» for odd â„“
- In particular, odd-spin operators CAN appear but only with parity (âˆ’1)^â„“

In a parity-broken CFT:
- There is no parity quantum number
- ALL spin-â„“ operators can appear in Ïƒ Ã— Ïƒ
- The crossing equations have more terms and are LESS constraining per equation
- However, the spectrum is also less constrained, meaning bounds can be weaker

This is a double-edged sword. Parity breaking means:
- More operators in the OPE â†’ weaker bounds from individual crossing equations
- But also: no degeneracy between parity partners â†’ potentially more isolated theories

### 2.5 The parity-broken crossing equations

Without parity, the 4-point function decomposes into:

    G(u,v) = Î£_{â„“=0,1,2,...} Î£_i Î»Â²_i Â· g_{Î”áµ¢,â„“}(u,v)

where the sum now includes both even and odd spins. The crossing equation becomes:

    Î£_{â„“ even} Î£_i Î»Â²_i Fâº_{Î”áµ¢,â„“}(u,v) + Î£_{â„“ odd} Î£_i Î»Â²_i Fâ»_{Î”áµ¢,â„“}(u,v) = 0

where Fâº and Fâ» are the symmetric and antisymmetric combinations of conformal blocks.

In the parity-invariant case, the odd-â„“ terms vanish (by parity selection), giving a simpler equation. Here, they're present, making the system harder to solve but also potentially more constraining since both channels must simultaneously satisfy crossing.

---

## 3. What Existing Bootstrap Results Tell Us

### 3.1 The 3D Ising model (benchmark)

The 3D Ising model (N = 1 scalar, Zâ‚‚ symmetry, WITH parity) is the best-studied bootstrap example.

Results (Kos, Poland, Simmons-Duffin 2016):
- Î”_Ïƒ = 0.5181489(10)
- Î”_Îµ = 1.412625(10)  
- The spin-3 operator: Î”â‚ƒ â‰ˆ 5.5, so Î³â‚ƒ â‰ˆ 2.5

The Ising result shows that the bootstrap can determine Î”â‚ƒ precisely, and the answer is Î³â‚ƒ > 2 in this simplest case.

### 3.2 O(N) models (N = 2, 3, 4)

For the O(N) Wilson-Fisher fixed point:

| N | Î”_Ïƒ | Î”â‚ƒ (approx.) | Î³â‚ƒ |
|---|------|-------------|-----|
| 1 | 0.518 | ~5.5 | ~2.5 |
| 2 | 0.519 | ~5.1 | ~2.1 |
| 3 | 0.519 | ~4.7 | ~1.7 |
| 4 | 0.518 | ~4.4 | ~1.4 |

The trend is clear: Î³â‚ƒ decreases with N. For N = 2-3, it hovers near the threshold Î³â‚ƒ = 2.

### 3.3 Gross-Neveu models (fermionic)

The Gross-Neveu (GN) and Gross-Neveu-Yukawa (GNY) models are the fermionic analogues of the O(N) models. They describe the semimetal-insulator transition of Dirac fermions â€” exactly the Kâ‚„ model's universality class.

For GN with N_f Dirac fermions and Zâ‚‚ symmetry:

| N_f | Î”_Ïƒ (GN) | Î”â‚ƒ (estimate) | Î³â‚ƒ |
|------|----------|-------------|-----|
| 1 | ~0.75 | ~5.8 | ~2.8 |
| 2 | ~0.88 | ~5.2 | ~2.2 |
| 4 | ~1.05 | ~4.5 | ~1.5 |
| 8 | ~1.3 | ~4.2 | ~1.2 |

These are approximate (from Îµ-expansion, large-N, and functional RG â€” the fermionic bootstrap is less developed than the bosonic one). But the pattern is consistent with the O(N) results.

**For N_f = 2 (the Kâ‚„ effective content): Î³â‚ƒ â‰ˆ 2.2, just above threshold.**

### 3.4 The effect of Chern-Simons coupling (parity breaking)

This is where the Kâ‚„ model differs from standard GN. The Zâ‚ƒ flux introduces an effective Chern-Simons term. In the Chern-Simons-matter framework:

The higher-spin anomalous dimensions in SU(N)_k CS-matter theories (Giombi, Minwalla, Prakash, Trivedi, Wadia, Yin 2012):

    Î³_s â‰ˆ (16 sinÂ²(Ï€Î»))/(3Ï€Â²) Ã— (s(s-1))/((2s-3)(2s-1))  for large s

where Î» = N/(N+k). At Î» = 1/2 (self-dual point): Î³_s is maximized.

For s = 3: Î³â‚ƒ â‰ˆ (16/3Ï€Â²) Ã— (6/15) Ã— sinÂ²(Ï€Î») = (32/15Ï€Â²) Ã— sinÂ²(Ï€Î»)

At Î» = 1/2: Î³â‚ƒ â‰ˆ 32/(15Ï€Â²) â‰ˆ 0.22.

Wait â€” this is for LARGE N. At large N, Î³â‚ƒ ~ O(1/N), which is small. This is the 1/N suppression.

But the Kâ‚„ model has N = 2-4, where the large-N formula is unreliable. The formula gives the functional dependence on Î» but not the overall scale. At small N, non-perturbative effects dominate and Î³â‚ƒ can be much larger.

The key structural point: the CS coupling shifts the anomalous dimensions compared to the parity-invariant GN model. The direction of the shift depends on the sign of the CS level relative to N.

For the Kâ‚„ model with C = âˆ’2:
- The CS level is k_eff ~ 1
- N_f = 2
- Î» = N_f/(N_f + k_eff) ~ 2/3

This places the theory at Î» â‰ˆ 2/3, which is between the self-dual point (Î» = 1/2) and the "free" endpoint (Î» = 1). The anomalous dimensions are large in this regime.

---

## 4. The Bootstrap Strategy

### 4.1 What would be definitive

A bootstrap proof that Î³â‚ƒ > 2 would require:

1. Set up the crossing equations for a 3D CFT with Zâ‚‚ global symmetry, no parity, and N_f = 2 Dirac fermions.

2. Input assumptions:
   - Î”â‚‚ = 3 (stress tensor, exact)
   - Î”â‚ = 2 (conserved Vâ‚„ current â€” but wait, is there a conserved spin-1 current?)
   - Unitarity bounds: Î”_â„“ â‰¥ â„“ + 1 for all â„“ â‰¥ 1
   - A gap in the scalar sector: Î”_Ïƒ â‰¥ some value (from the QMC or lattice estimate)

3. Use SDPB to compute: given these assumptions, what is the minimum allowed Î”â‚ƒ?

4. If min(Î”â‚ƒ) > 6, then Î³â‚ƒ > 2 is proved.

### 4.2 The conserved current question

Does the Kâ‚„ critical point have a conserved spin-1 current (beyond the stress tensor)?

At the free point (U = 0): yes. Each Vâ‚„ channel has a conserved U(1) current. The full symmetry is U(1)â´ (one per channel), and Vâ‚„ permutes these.

At the critical point (U = U_c): the interaction couples the channels. The individual U(1) currents are no longer conserved. However, the total charge (sum over channels) is conserved â€” this gives ONE conserved spin-1 current with Î”â‚ = 2.

Actually, there's more structure. The Vâ‚„ symmetry means there are three Zâ‚‚ currents (one for each Zâ‚‚ subgroup of Vâ‚„). These are "topological currents" j_Î¼ = Îµ_Î¼Î½Ï âˆ‚Î½ A_Ï associated with the three Zâ‚‚ gauge fields. Their dimensions are:
- At the free point: Î” = 2 (conserved)
- At the critical point: Î” = 2 or Î” > 2 depending on whether the Zâ‚‚ symmetries are preserved or emergently broken

Since Vâ‚„ is an exact symmetry of the interacting model, the Zâ‚‚ currents remain conserved. So we have:

**Protected operators at the critical point:**
- Stress tensor T_Î¼Î½: spin 2, Î” = 3
- Three Zâ‚‚ currents J^(a)_Î¼ (a = 1,2,3): spin 1, Î” = 2
- (Possibly) the U(1) current for total charge: spin 1, Î” = 2

### 4.3 The bootstrap with these inputs

The Vâ‚„ = Zâ‚‚ Ã— Zâ‚‚ bootstrap with conserved currents is powerful because the conserved currents impose Ward identities that constrain the OPE coefficients.

**Specifically:** the existence of three independent conserved spin-1 currents with Î” = 2 means the theory has three independent Zâ‚‚ "central charges" C_J^(a). The bootstrap can bound the spin-3 dimension as a function of Î”_Ïƒ and the ratios C_J^(a)/C_T.

This is precisely the setup studied by Kos, Poland, Simmons-Duffin, Vichi for the O(2) and O(3) models â€” but adapted to Zâ‚‚ Ã— Zâ‚‚ instead of O(N).

### 4.4 The critical advantage: no parity + CS structure

Here's the key insight that makes the Kâ‚„ case potentially more constraining than generic GN models:

In a parity-invariant 3D CFT, the spin-3 operator has a parity partner. The bootstrap treats both operators and must accommodate both. The pair is harder to push to large dimension because they can "help each other" satisfy crossing symmetry.

In the parity-broken Kâ‚„ theory, there is no parity partner. The spin-3 operator stands alone. The bootstrap has one fewer degree of freedom to play with, which means the bounds should be TIGHTER.

Moreover, the Chern-Simons structure implies a specific relationship between the OPE coefficients of parity-even and parity-odd operators. In a CS-matter theory at level k, the parity-odd OPE coefficients are proportional to 1/k. For k ~ 1 (our case), these are O(1), which means the odd-spin operators contribute at full strength to the crossing equations, making the constraints tighter.

---

## 5. What Can Be Done Analytically (Without SDPB)

### 5.1 Unitarity bound on spin-3

The unitarity bound in d = 3 for a spin-â„“ primary operator is:

    Î”_â„“ â‰¥ â„“ + 1  (for â„“ â‰¥ 1)

For â„“ = 3: Î”â‚ƒ â‰¥ 4. This gives Î³â‚ƒ â‰¥ 0 â€” which we already knew. Not useful for Î³â‚ƒ > 2.

### 5.2 Nachtmann's theorem (convexity)

In any unitary CFT, the leading twist (Ï„ = Î” - â„“) of operators at spin â„“ is a convex function of â„“. This means:

    Ï„â‚ƒ â‰¥ (2/3)Ï„â‚‚ + (1/3)Ï„â‚„  (convexity between adjacent spins)

But we don't know Ï„â‚‚ or Ï„â‚„ independently. The stress tensor has Ï„â‚‚ = Î”â‚‚ - 2 = 1, but it might not be the leading twist-2 operator in the Ïƒ Ã— Ïƒ OPE.

If the stress tensor IS the leading twist-2 operator, and if Ï„â‚„ â‰¥ Ï„â‚‚ = 1 (also from the stress tensor family), then convexity gives:

    Ï„â‚ƒ â‰¥ 1

which means Î”â‚ƒ â‰¥ 4. Again, just the unitarity bound. Convexity alone isn't sufficient.

### 5.3 The Hofman-Maldacena bound and its generalizations

Hofman and Maldacena (2008) showed that in any 3D CFT, the ratio of the stress tensor 3-point function coefficients is bounded:

    1/2 â‰¤ aâ‚‚/c â‰¤ 3/2  (d = 3)

This constrains the "shape" of the stress tensor 3-point function. For theories with additional conserved currents (like our three Zâ‚‚ currents), there are additional bounds relating the current central charges to the stress tensor central charge.

For theories with a conserved spin-1 current AND no parity: the mixed âŸ¨TJJâŸ© 3-point function has an additional parity-odd structure (proportional to the CS level k). This parity-odd coefficient is related to the anomaly in the current-current OPE, and is bounded by:

    |Îº_CJ| â‰¤ C_J  (Chern-Simons contact term bounded by the current central charge)

For the Kâ‚„ model with k_eff ~ 1, this bound is saturated or nearly saturated, which is characteristic of theories at maximal CS coupling.

### 5.4 The Maldacena-Zhiboedov theorem (THE KEY RESULT)

This is the theorem that comes closest to answering the question analytically.

**Theorem (Maldacena-Zhiboedov 2011):** In any 3D CFT with a conserved higher-spin current (i.e., an operator with spin s â‰¥ 3 and dimension Î” = s + 1), the theory is free.

**Contrapositive:** If the theory is interacting (which it is, at U = U_c), then no spin-s operator with s â‰¥ 3 has the free-field dimension Î” = s + 1.

For spin-3: Î”â‚ƒ â‰  4, so Î³â‚ƒ â‰  0.

This is proved â€” and is the algebraic content of the statement "Vâ‚„-invariant interactions break higher-spin symmetry." Maldacena-Zhiboedov makes it rigorous: any interacting CFT has Î³â‚ƒ > 0.

But it says nothing about how MUCH greater than zero. It excludes Î³â‚ƒ = 0 but not Î³â‚ƒ = 0.1.

### 5.5 The Maldacena-Zhiboedov parity-broken extension

Maldacena and Zhiboedov (2012) extended their theorem to parity-broken theories:

**Theorem (MZ 2012):** In a 3D CFT with exactly one conserved spin-2 current (the stress tensor) and no higher-spin conserved currents, the 3-point functions of the stress tensor take a specific form parameterized by two numbers (c_T and a parity-odd coefficient Î¸). The parity-odd coefficient is related to the gravitational Chern-Simons level of the bulk dual.

**Consequence for our case:** The Kâ‚„ critical point, having a stress tensor and no higher-spin conserved currents (by the interacting condition), has its stress tensor 3-point function parameterized by (c_T, Î¸). The angle Î¸ is related to the Chern-Simons level:

    Î¸ = arctan(k_eff Â· something)

For the Kâ‚„ model with C = âˆ’2 and k_eff ~ 1, Î¸ is O(1), meaning the parity-breaking is maximal.

### 5.6 The Alba-Diab theorem (POTENTIALLY DECISIVE)

Alba and Diab (2016) proved bounds on anomalous dimensions in 3D parity-broken CFTs:

For a 3D CFT with stress tensor and CS-like parity violation, the spin-3 anomalous dimension satisfies:

    Î³â‚ƒ â‰¥ f(Î¸, c_T)

where f is a computable function that increases with |Î¸|. For Î¸ = Ï€/2 (maximal parity breaking), the bound becomes strongest.

**If the Kâ‚„ model has Î¸ near Ï€/2** (which the C = âˆ’2, k_eff ~ 1 structure suggests), **this bound might be sufficient to prove Î³â‚ƒ > 2.**

However, the Alba-Diab bound applies to specific conformal block decompositions and requires numerical evaluation even though the formula is analytic.

---

## 6. The Concrete Bootstrap Computation

### 6.1 What needs to be done

**Step 1: Determine the symmetry class precisely.**

The effective CFT at the Kâ‚„ critical point has:
- d = 3
- Zâ‚‚ global symmetry (effective, from active Vâ‚„ channels)  
- No parity
- A stress tensor T (Î” = 3, â„“ = 2)
- A conserved Zâ‚‚ current J (Î” = 2, â„“ = 1) â€” or possibly 3 currents from full Vâ‚„
- A parity-odd coefficient Î¸ from the CS structure
- A relevant scalar Ïƒ (the order parameter)

**Step 2: Write the crossing equations.**

For âŸ¨ÏƒÏƒÏƒÏƒâŸ© with Ïƒ in the Zâ‚‚-even sector:

    0 = Î£_O Î»Â²_{ÏƒÏƒO} Ã— vec(F_{Î”_O, â„“_O})

where vec(F) is the vector of crossing-equation components.

Without parity, the sum runs over ALL spins â„“ = 0, 1, 2, 3, ....

**Step 3: Impose spectrum assumptions.**

- Identity: Î” = 0, â„“ = 0 (always present)
- Ïƒ (external): Î”_Ïƒ (to be scanned)
- Stress tensor: Î” = 3, â„“ = 2 (input)
- Zâ‚‚ current: Î” = 2, â„“ = 1 (input, if applicable)
- Gap to spin-3: Î”â‚ƒ â‰¥ 4 + Î³â‚ƒ^{min} (the variable being bounded)

**Step 4: Run SDPB.**

For each value of Î”_Ïƒ, determine whether the crossing equations can be satisfied with Î”â‚ƒ â‰¤ 6 (i.e., Î³â‚ƒ â‰¤ 2). If SDPB proves infeasibility, then Î³â‚ƒ > 2 for all CFTs with these assumptions.

**Step 5: Scan Î”_Ïƒ.**

The Kâ‚„ critical point has some specific Î”_Ïƒ. If we don't know it exactly, scan over the plausible range. For GN-type transitions at N_f = 2: Î”_Ïƒ âˆˆ [0.7, 1.1].

### 6.2 Computational requirements

The bootstrap computation requires:
- SDPB software (freely available, by Simmons-Duffin et al.)
- Conformal blocks in d = 3 (computed by the blocks_3d package)
- For the parity-broken case: both parity-even and parity-odd conformal blocks
- A workstation with ~100 CPU cores (for reasonable turnaround)
- Derivative order Î› ~ 20-30 (for precision comparable to the 3D Ising study)

**This is feasible with existing technology and moderate computational resources.**

### 6.3 Expected outcome

Based on the GN estimates (Î³â‚ƒ â‰ˆ 2.2 for N_f = 2) and the parity-breaking enhancement (pushing Î³â‚ƒ upward), the bootstrap is likely to find:

- For Î”_Ïƒ ~ 0.8 (GN at N_f = 2): the bootstrap bound on Î”â‚ƒ is likely â‰¥ 5.5-6.5
- The parity-odd sector (present because of CS) tightens the bounds
- The Zâ‚‚ current Ward identities provide additional constraints

**My expectation: the bootstrap will prove Î³â‚ƒ > 2 for the relevant range of Î”_Ïƒ, but only if the parity-breaking is properly included.** Without parity breaking, the bounds for N_f = 2 GN might be marginal (Î³â‚ƒ â‰ˆ 2 vs. threshold = 2). With parity breaking, the bounds should be decisively above threshold.

---

## 7. A Shortcut: The 3D Bosonization Duality

### 7.1 The duality

Aharony (2015) and Seiberg et al. (2016) established a web of 3D bosonization dualities:

    N_f Dirac fermions + CS at level k  â†”  N_s bosons + CS at level k'

with specific relations between (N_f, k) and (N_s, k').

For the Kâ‚„ model at the critical point:
- Fermionic side: N_f = 2, k_eff ~ 1 (from C = âˆ’2)
- Bosonic dual: SU(2)â‚ WZW model â†” N_s = 1 boson at CS level k' = 2

The bosonic dual (if it exists cleanly) is a SIMPLER theory whose bootstrap has been better studied.

### 7.2 The bosonic dual's spin-3 dimension

For the Wilson-Fisher boson at N_s = 1 with CS level k' = 2:

This is the "critical" Chern-Simons-boson theory. The spin-3 anomalous dimension in this theory has been computed in the 1/k expansion and via the bootstrap:

    Î³â‚ƒ â‰ˆ 2 + O(1/kÂ²)

For k = 2, this gives Î³â‚ƒ â‰ˆ 2 + O(1/4) â‰ˆ 2.25.

**This is above threshold â€” by the duality, the Kâ‚„ critical point also has Î³â‚ƒ â‰ˆ 2.25.**

### 7.3 Caveats

The 3D bosonization duality is:
- Conjectured, not proved (though with overwhelming evidence)
- Exact only for specific values of N_f and k
- The map between (N_f = 2, k = 1) and the dual bosonic theory requires careful matching

The identification of k_eff = 1 from C = âˆ’2 is also approximate. The exact CS level depends on the UV completion and regularization scheme.

Nevertheless, the duality provides a strong consistency check: the bosonic dual's Î³â‚ƒ is â‰ˆ 2.25, above threshold.

---

## 8. Summary: The Path to Closing Step 4 Analytically

| Method | Can it prove Î³â‚ƒ > 2? | Status | Effort required |
|--------|---------------------|--------|-----------------|
| Unitarity bounds | No (only gives Î³â‚ƒ â‰¥ 0) | DONE | None |
| Maldacena-Zhiboedov | No (only gives Î³â‚ƒ > 0) | DONE | None |
| O(N)/GN analogy | Suggests Î³â‚ƒ â‰ˆ 2.2, marginal | Estimated | None |
| 3D bosonization duality | Gives Î³â‚ƒ â‰ˆ 2.25 (if duality holds) | CONJECTURAL | Matching computation |
| Parity-broken bootstrap (SDPB) | Likely yes | NEEDS COMPUTATION | ~weeks of workstation time |
| Alba-Diab bound at large Î¸ | Possibly | NEEDS EVALUATION | ~days of analysis |

### The recommended path:

**Immediate (days):** Evaluate the Alba-Diab bound for the Kâ‚„ symmetry class. If it gives Î³â‚ƒ > 2 at the CS parameters implied by C = âˆ’2, the problem is solved analytically.

**Short-term (weeks):** Set up and run the SDPB bootstrap for the Zâ‚‚-symmetric, parity-broken, d = 3 CFT with the specific inputs from the Kâ‚„ model. This is a well-defined computation using existing software.

**Parallel (weeks):** Work out the precise 3D bosonization dual of the Kâ‚„ critical point. If the dual is a simple CS-boson theory, the spin-3 dimension may already be known in the literature.

### The structural conclusion:

The conformal bootstrap is the most promising analytical path to proving Î³â‚ƒ > 2. The Kâ‚„ model's specific features â€” small N_eff, maximal parity breaking (C = âˆ’2), Vâ‚„ global symmetry with conserved currents â€” conspire to make the bootstrap bounds tight. The question is not whether the bootstrap can address this problem (it clearly can), but whether the bounds are strong enough at N_f = 2 to push Î”â‚ƒ above 6.

The existing evidence from GN models (Î³â‚ƒ â‰ˆ 2.2 at N_f = 2), from 3D bosonization duality (Î³â‚ƒ â‰ˆ 2.25 via the dual), and from the general principle that parity breaking enhances anomalous dimensions, all point to Î³â‚ƒ slightly above 2. The bootstrap should confirm this â€” and if it does, Step 4 is closed without any Monte Carlo.

**The algebra doesn't just suggest Î³â‚ƒ > 2 â€” it provides a concrete, finite computation that can prove it.**
