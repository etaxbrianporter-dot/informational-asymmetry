# Informational Asymmetry Program

**From three axioms to the Standard Model and gravity.**

Matching combinatorics on complete graphs K<sub>2n</sub> uniquely determine Lorentzian signature, the gauge group, three fermion generations, the Higgs mass, and the fermion mass hierarchy — with zero free parameters for particle physics and one for gravity.

[![Verification](https://img.shields.io/endpoint?url=https://etaxbrianporter-dot.github.io/informational-asymmetry/badge.json)](https://etaxbrianporter-dot.github.io/informational-asymmetry/)

**Brian Porter · 2026 · Working in public**

---

## The Claim

Three informational axioms — asymmetric transitions, relational comparison, finite capacity — force the saturation equation (2n−1)!! = I(n), which has exactly two solutions. These solutions generate:

| Level | Graph | Physical Content | Free Parameters |
|-------|-------|-----------------|-----------------|
| K₄ | 4 vertices, 3 matchings | Lorentzian signature, arrow of time, Einstein equations | 0 |
| K₆ | 6 vertices, 15 matchings | Higgs mass (126.1 GeV), gauge group, three generations | 0 |
| K₈ | 8 vertices, 105 matchings | Yukawa hierarchy 415:135:1, ρ₂ selection rule | 0 |
| Gravity | — | M_Pl, Λ_CC | 1 (Λ) |

## Verify Everything

```bash
cd verification_suite
pip install numpy scipy
python run_all.py
```

117 independently verified claims. Every numerical result in every paper has a test. Uses only standard Python libraries.

## Repository Structure

```
├── index.html                    # Main website
├── papers/                       # LaTeX sources + compiled PDFs
│   ├── paper1.tex                # K₄ → Lorentzian signature
│   ├── paper2.tex                # Cosmic birefringence predictions
│   └── paper3.tex                # K₈ fermion mass hierarchy
├── docs/                         # Working documents (markdown)
│   ├── executive_summary.md
│   ├── spectral_budget_identity.md
│   ├── galois_disjointness.md
│   ├── ...
│   └── c_derivation_termination.md   # Killed approaches included
├── verification_suite/           # Reproducibility test suite
│   ├── run_all.py
│   ├── lib/                      # Core computational libraries
│   └── tests/                    # 22 test scripts
├── code/                         # Research computation scripts
│   ├── k6_tools.py
│   ├── gz_basis.py
│   └── ...
├── scripts/                      # Build pipeline
│   ├── build_site.py             # Markdown → HTML converter
│   ├── dev.py                    # Local dev server with auto-rebuild
│   └── add_frontmatter.py       # Metadata helper
└── .github/workflows/
    ├── build-deploy.yml          # CI: verify → build → deploy to Pages
    └── verify.yml                # PR gate: verification suite must pass
```

## How It Works

Every push to `main`:
1. **Verification suite runs** — if any claim fails, the build breaks
2. **Markdown docs are converted to HTML** — with git timestamps, section grouping, consistent styling
3. **Site deploys to GitHub Pages** — including papers, docs, verification suite, and code

Every pull request:
- Verification suite runs as a gate — can't merge if claims fail

## Local Development

```bash
pip install markdown pygments pyyaml numpy scipy
python scripts/dev.py
# → Serves at localhost:8000, auto-rebuilds on file changes
```

## Working in Public

This repository contains everything: successful results, failed derivations, bug discoveries, and honest assessments. The killed approaches section documents paths that were computationally falsified and terminated. Test 22 demonstrates two bugs we found in our own code.

This is how science should work.

## Status

See the [full status table](https://bporter.github.io/informational-asymmetry/#status) on the website, or `docs/progress_addendum_feb18.md` in this repo.

## License

All working documents: CC-BY-4.0. Code: MIT.
