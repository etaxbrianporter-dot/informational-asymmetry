---
title: "The Riemann Hypothesis and the Polynomial/Series Wall"
section: "Mathematical Structure"
status: "terminated"
---

# The Riemann Hypothesis and the Polynomial/Series Wall

## Paper V of the Informational Asymmetry Program

**Brian Porter — February 2026**

---

## Abstract

We investigate two independent paths from the Informational Asymmetry framework toward the Riemann Hypothesis. The first repackages the Connes semilocal trace formula as an informational asymmetry measure δ_NC, establishing a dictionary between bounded asymmetry and the location of zeta zeros. The second exploits the structural identity between the matching chain K₂ₙ and the function-field/number-field transition: three of five structural features (Euler product factorization, Galois disjointness, and the polynomial/series wall itself) are mathematically identical. Both paths terminate at the same obstruction — infinite-dimensional positivity on the adèle class space — which we identify as the polynomial/series wall. The wall is mapped, its boundary precisely located, and its impassability demonstrated from three independent frameworks (motivic cohomology, Connes' trace formula, the Langlands program). The matching chain's polynomial side, where the Standard Model lives, is fully computable; the series side, where the Riemann Hypothesis lives, is not accessible by the methods of this program. We document this termination as a structural result: the same wall that separates finite algebraic computability from infinite analytic truth governs both number theory and the physics of the continuum limit.

---

## 1. Two Paths to the Same Wall

The Informational Asymmetry program produces, as a byproduct of its matching combinatorics, structures that are eerily parallel to the analytic number theory surrounding the Riemann Hypothesis. We pursued two paths to see where they led. Both terminated at the same point.

**Path 1 (Connes-Asymmetry):** Repackage the semilocal trace formula discrepancy as a third-order asymmetry measure δ_NC, using the information-geometric language natural to the program.

**Path 2 (Matching Chain):** Observe that the matching chain K₄ → K₆ → K₈ → ··· has a secular factorization structurally identical to an Euler product, with proved Galois disjointness playing the role of ℚ-independence of log primes.

Both arrive at the polynomial/series wall: the distinction between finite products (polynomials, proved) and infinite products (series, open).

---

## 2. Path 1: The Asymmetry Measure δ_NC

### 2.1 Definition

Let (A, H, D) be the spectral triple for the adèle class space X = ℚ× \ A_ℚ / Ẑ× as in Connes' construction. Define the third-order asymmetry gap:

$$\delta_{NC}(h) = \frac{\left| \text{Tr}(R_\Lambda U(h)) - \left(2h(1)\log'\Lambda + \sum_v \int'_{k_v^*} h(u^{-1})|1-u| \, d^*u \right)\right|}{(\log \Lambda)^{3/2}}$$

where R_Λ = P̂_Λ P_Λ is the cutoff operator, U(h) is the unitary representation for test function h, and the integral is the principal value in the sense of Connes (1998). The (log Λ)^{3/2} normalization is imported from the exponential family analogy: |T|/g^{3/2} for Fisher information g.

### 2.2 The conjecture

**Bounded Asymmetry ⟺ RH:** The Riemann Hypothesis holds for all L-functions with Größencharakter if and only if δ_NC(h) remains bounded as Λ → ∞ for all suitable h.

### 2.3 What survives

**The δ_NC definition is valid.** It is a mathematically clean translation of Connes' trace formula discrepancy into information-geometric language. The dictionary is exact:

| Exponential family | Connes trace formula |
|--------------------|---------------------|
| Sample space | Adèle class space |
| Sufficient statistic T | Tr(R_Λ U(h)) |
| Expected value | Semilocal terms |
| Fisher information g | (log Λ)^{3/2} normalization |

**The off-line growth lemma is correct.** If ρ is a zero with Re(ρ) = 1/2 + ε > 1/2, the von Mangoldt explicit formula gives ψ(x) − x ~ x^{1/2+ε}/ρ, so |δ| diverges logarithmically with x. This is not new — it is a classical consequence of the explicit formula. But the repackaging as "asymmetry divergence" is clean and correct.

### 2.4 What collapses

**The π/3 bound collapses.** The proposed universal saturation bound δ_NC ≤ π/3 was imported from finite-dimensional exponential families, where the cumulant generating function provides a natural scale. But the Euler product is not an exponential family. The adèle class space is infinite-dimensional. There is no cumulant generating function, no natural Fisher information, and no basis for the specific value π/3. The bound is an analogy, not a derivation.

**The conjecture reduces to Connes' open problem.** The assertion "bounded δ_NC ⟺ RH" is, when unpacked, equivalent to the statement that the Weil positivity condition holds on the adèle class space. This is precisely the open problem identified by Connes in 1998 — proving positivity in infinite dimensions. Calling the discrepancy "asymmetry" and normalizing by (log Λ)^{3/2} changes the packaging, not the content.

### 2.5 Termination

Path 1 terminates at infinite-dimensional positivity. The asymmetry framework provides a valid reinterpretation of Connes' trace formula in information-geometric language, but this reinterpretation is orthogonal to the analytic obstruction that prevents proving RH.

---

## 3. Path 2: The Matching Chain as Geometric Euler Product

### 3.1 Secular factorization

The Gelfand pair branching rules show that the physical eigenspace V_{(2n−2,2)} decomposes under hub restriction, with Δd(k) = 4k − 5 new dimensions at level k. The bordered eigenproblem produces a secular factor S_k(λ) — a polynomial of degree 4k − 5 — whose roots are the eigenvalues introduced at level k.

The full characteristic polynomial factorizes as a finite product:

$$\text{char}(G_{\text{phys}}(n), \lambda) = \prod_{k=2}^{n} S_k(\lambda)$$

This is structurally identical to the Euler product factorization of L-functions.

### 3.2 Five structural identities

| Feature | Function field L-function | Matching chain | Status |
|---------|--------------------------|----------------|--------|
| Euler product factorization | ∏_P det(1 − χ(P)u^{deg P})^{−1} | ∏_{k=2}^{n} S_k(λ) | **Identical** |
| Frobenius / bordered update | Finite-rank perturbation | Finite-rank perturbation | Partial |
| Galois disjointness | log primes are ℚ-independent | ℚ(λ₆) ∩ ℚ(λ₈) = ℚ | **Identical** |
| Functional equation / budget | Constrains zero locations | Budget identity constrains spectral weight | Partial |
| The Wall | Polynomial ↔ series | Finite product ↔ infinite product | **Identical** |

Three of five features are structurally identical. The matching chain hits the same wall as the Riemann Hypothesis.

### 3.3 Galois disjointness as ℚ-independence

The Galois disjointness proved for ℚ(λ_vac(3)) and ℚ(λ_vac(4)) is the matching chain's version of the ℚ-independence of log primes:

| Feature | RH direction | Matching chain |
|---------|-------------|----------------|
| Independent objects | log p for each prime p | λ_vac(n) for each level n |
| Independence type | ℚ-linear independence | Galois disjointness |
| Consequence | Adelic torus T^∞ is infinite-dimensional | Compositum has ∏ deg_k growth |
| Prevents | Critical line from closing | Algebraic closure of the limit |
| Status | Proved (unique factorization) | Proved for n = 3, 4 |

The compositum F_n = ℚ(λ_vac(3), ..., λ_vac(n)) has degree [F_n : ℚ] = ∏_{k=3}^{n} [ℚ(λ_vac(k)) : ℚ], growing super-exponentially. The known data:

| Level | deg(λ_vac) | Discriminant primes | Number field |
|-------|-----------|-------------------|-------------|
| K₆ (n=3) | 2 | {5} | ℚ(√5) |
| K₈ (n=4) | 6 | {2, 7, 43, 421} | Galois group C₂ ≀ C₃ |
| K₁₀ (n=5) | ≤ 8 | Unknown | — |

### 3.4 The budget identity as proved spectral constraint

At each finite level, the budget identity (proved from the Johnson scheme) constrains the spectral weight distribution:

$$\lambda_{\text{mid}} \cdot d_{\text{phys}} = 2(n-1) \cdot \lambda_{\text{max}}$$

The self-focusing parameter η(n) = 2(n−1)/(2n−1) → 1 as 1/n², meaning physical observables converge polynomially even as algebraic complexity grows super-exponentially. This is the matching chain's analog of the functional equation — a constraint that spectral data must satisfy, proved at every finite level.

### 3.5 The function field proof of concept

Over function fields F_q(t), the V₄ channel decomposition recovers integer winding numbers (= deg(L_i)). The critical circle is compact (period 2π/log q) because deg(P) ∈ ℤ. The Weil conjectures prove zero-pinning. The polynomial/series wall is the load-bearing distinction: the bridge works over function fields (where everything is polynomial) and breaks over number fields (where everything is series).

---

## 4. The Wall

### 4.1 The precise obstruction

The polynomial/series wall is, at root, the distinction between finite-dimensional and infinite-dimensional algebraic objects:

- **Polynomial side:** Each matching chain level n has char(G_phys(n)) as a polynomial of known degree n(2n−3). Budget identity proved. All eigenvalues algebraic. Galois groups computable. Spectral questions decidable. This is the function field side.

- **Series side:** The limit n → ∞ produces an infinite product. Properties of the limit are open. Algebraic closure is blocked by Galois disjointness. Convergence of λ_vac(n) is numerical but unproved. Whether the limit is a period (in the Kontsevich-Zagier sense) depends on motivic structure that has not been established.

The matching chain and number theory hit the same wall because they face the same mathematical phenomenon: attempting to take an inverse limit of finite-dimensional algebraic objects where the connecting maps introduce genuinely new algebraic content at each step, and the limit is not algebraic.

### 4.2 Three frameworks, one wall

Three known frameworks attempt to cross the polynomial/series wall for the Riemann Hypothesis:

**Motivic cohomology.** The correct framework in principle — L-functions as det(1 − u·Frob|H¹) over function fields. But the motivic Frobenius over ℚ is not known to exist. The Langlands-Rapoport conjecture and the standard conjectures on algebraic cycles remain open.

**Connes' trace formula.** Identifies the correct positivity condition on adèlic space, reducing RH to positivity of a specific functional in infinite dimensions. But positivity is unproved, and the infinite-dimensional analysis required appears to lie beyond current techniques.

**The Langlands program.** Unifies automorphic L-functions across number fields and function fields. Even where proved (Wiles for GL₂/ℚ), does not give RH. The functoriality conjectures, if established in full generality, would imply far more than RH — but they are themselves open.

The V₄ framework from the Informational Asymmetry program is orthogonal to the wall. V₄ works identically on both sides (polynomial and series). The wall concerns the analytic character of L-functions, not their algebraic decomposition.

### 4.3 The wall as Gödelian boundary

The universe may be Gödelian in a precise sense: the fundamental structure (matching algebra at finite levels) is decidable and complete, but the emergent physics (continuum gauge theory, the mass gap, the nature of λ_∞) contains truths that transcend any finite algebraic description.

The polynomial/series wall isn't a technical obstacle we haven't figured out how to cross. It may be an incompleteness barrier — a genuine limit on what finite algebraic systems can say about their own continuum limits.

The critical observation: the physics doesn't care. The self-focusing theorem says physical observables converge polynomially (Δη ~ 1/2n²) even as algebraic complexity diverges super-exponentially. The Standard Model lives entirely on the polynomial side. You don't need to cross the wall to build the universe.

---

## 5. The Riemann Hypothesis from the Matching Chain Perspective

### 5.1 What we see

The matching chain generates a sequence of number fields ℚ(λ_vac(n)) that are pairwise Galois-disjoint (proved for n = 3, 4), with discriminant primes reflecting both structural symmetry (5 from Z₅ at K₆, 7 from Z₇ at K₈) and spectral arithmetic (43 and 421 at K₈, with no obvious combinatorial origin).

The eigenvalues converge numerically. The budget identity holds at every finite level. The spectral weight self-focuses: η(n) → 1. All of this is proved, computable, and on the polynomial side of the wall.

### 5.2 What we conjecture

**Matching Chain Disjointness:** For all n ≠ m, ℚ(λ_vac(n)) ∩ ℚ(λ_vac(m)) = ℚ. The number fields at successive levels are pairwise Galois-disjoint.

If true, this implies the limit λ_∞ = lim λ_vac(n), if it exists, is transcendental over ℚ — no polynomial with integer coefficients captures it. The polynomial/series wall is genuine: the limit lies strictly on the series side.

**Period structure:** If the matching chain has motivic structure, the Kontsevich-Zagier period conjecture would force λ_∞ to be a period (expressible as an integral of an algebraic form over a semi-algebraic domain). This would place λ_∞ in the same class as π, log 2, and ζ(3), but exclude numbers like e, 1/π, and the Euler-Mascheroni constant γ.

### 5.3 What we do not see

A path to proving the Riemann Hypothesis. The asymmetry measure δ_NC is a valid repackaging. The matching chain exhibits the correct structural features. But neither provides new analytic content toward the positivity proof. The wall is identified, mapped, and impassable by these methods.

---

## 6. A New Question: Is the Yang-Mills Mass Gap Algebraic?

The matching chain perspective suggests a question about the Yang-Mills mass gap that appears to be unasked in the literature:

At lattice spacing a, the mass gap Δ(a) is (in principle) an algebraic number in a number field K(a) determined by the gauge group and lattice structure. As a → 0, the number fields proliferate and Δ(a) → Δ_YM. The question: is Δ_YM algebraic or transcendental?

If the Galois disjointness pattern from the matching chain extends to lattice gauge theory, then Δ_YM is transcendental — the continuum limit crosses the polynomial/series wall. The arithmetic structure of the lattice mass gaps carries physical information (via the pro-finite Galois group of the compositum), even though the individual number fields diverge.

This connects the two millennium problems: the Yang-Mills mass gap lives on the series side of the same wall that the Riemann Hypothesis asks about. The nature of that wall — whether it is a Gödelian boundary or merely a reflection of current mathematical limitations — is the deeper question.

---

## 7. Terminal Assessment

### 7.1 Comparison table

| Component | Status | Load-bearing? |
|-----------|--------|---------------|
| δ_NC definition | **Survives** (as translation) | No — reparametrization of Connes |
| Bounded asymmetry ⟺ RH | **Terminal** | Reduces to Connes' open problem |
| π/3 saturation bound | **Collapsed** | No derivation in infinite dimensions |
| Off-line growth lemma | **Proved** | No — classical (von Mangoldt) |
| Numerical validation | **Valid** | No — confirms known asymptotics |
| Euler product structure | **Identical** | Yes — structural identity, not proof |
| Galois disjointness | **Proved** (K₆, K₈) | Yes — matching chain version of wall |
| Budget identity as RH analog | **Partial** | Additive vs. multiplicative |
| Path to RH | **Terminal** | Hits polynomial/series wall |

### 7.2 What remains

The arithmetic direction terminates. Two independent paths (Connes-asymmetry and matching-chain Euler product) arrive at the same wall. What remains is exact (the dictionary, the structural identity, the Galois disjointness). What fails is identified (infinite-dimensional positivity). The boundary between them is precisely mapped.

The matching chain is not a route to the Riemann Hypothesis. It is a combinatorial object that exhibits the same structural architecture as the analytic number theory that surrounds the Hypothesis. The identity of this architecture is itself a result — it says something about the mathematical universe that the same wall governs both matching combinatorics and the distribution of primes.

Whether this is a coincidence, a consequence of some deeper categorical structure, or a genuine connection waiting for new mathematics to exploit, we cannot say. We document the termination and the precise boundary.

---

## References

1. Connes, A. "Trace formula in noncommutative geometry and the zeros of the Riemann zeta function." Selecta Math. 5, 29–106 (1999).
2. Weil, A. "Sur les courbes algébriques et les variétés qui s'en déduisent." Actualités Sci. Ind. 1041, Hermann (1948).
3. Deligne, P. "La conjecture de Weil: I." Publ. Math. IHÉS 43, 273–307 (1974).
4. Davenport, H. *Multiplicative Number Theory.* 3rd ed., Springer (2000).
5. Kontsevich, M. & Zagier, D. "Periods." In *Mathematics Unlimited — 2001 and Beyond*, Springer, 771–808 (2001).
6. Porter, B. "The Polynomial Wall: Matching Chain as Geometric Euler Product." Working document (2026).
7. Porter, B. "Galois Disjointness between Number Fields at Each Level." Working document (2026).
8. Porter, B. "Discriminant Atlas: Arithmetic Structure K₆ through K₁₆." Working document (2026).
