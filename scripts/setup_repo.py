#!/usr/bin/env python3
"""
setup_repo.py — Initialize the repository from existing project files.
=======================================================================
Run this once to populate docs/ and papers/ from your working files.
It copies and renames files into the expected structure.

Usage:
    python scripts/setup_repo.py /path/to/project/files
    python scripts/setup_repo.py   # Interactive: asks for path
"""
import sys
import shutil
from pathlib import Path

ROOT = Path(__file__).parent.parent

# Map project files to repo locations
FILE_MAP = {
    # Papers
    'paper1.tex': 'papers/paper1.tex',
    'paperII.tex': 'papers/paper2.tex',
    'paperIII.tex': 'papers/paper3.tex',
    
    # Code
    'k6_tools.py': 'code/k6_tools.py',
    'gz_basis.py': 'code/gz_basis.py',
    'computation_D1.py': 'code/computation_D1.py',
    'computation_D1_manybody.py': 'code/computation_D1_manybody.py',
    'aperture_analysis.py': 'code/aperture_analysis.py',
    'f_epsilon_arrow_v2.py': 'code/f_epsilon_arrow_v2.py',
    
    # Docs — working markdown files
    'K6_normalization_analysis.md': 'docs/K6_normalization_analysis.md',
    'normalization_derivation_v9.md': 'docs/normalization_derivation_v9.md',
    'normalization_derivation_v10.md': 'docs/normalization_derivation_v10.md',
    'K5_normalization_analysis.md': 'docs/K5_normalization_analysis.md',
    'higgs_from_graph_chain.md': 'docs/higgs_from_graph_chain.md',
    'binary_coloring_summary.md': 'docs/binary_coloring_summary.md',
    'doublet_selection_answer.md': 'docs/doublet_selection_answer.md',
    'k8_fermion_results.md': 'docs/k8_fermion_results.md',
    'k8_fermion_results_v2.md': 'docs/k8_fermion_results_v2.md',
    'k8_higgs_results.md': 'docs/k8_higgs_results.md',
    'R8_algebraic_analysis.md': 'docs/R8_algebraic_analysis.md',
    'c_derivation_termination_analysis.md': 'docs/c_derivation_termination_analysis.md',
    'c_derivation_terminated.md': 'docs/c_derivation_terminated.md',
    'light_emergence_theorem.md': 'docs/light_emergence_theorem.md',
    'D1MB_bug_analysis.md': 'docs/D1MB_bug_analysis.md',
    'spectral_budget_identity.md': 'docs/spectral_budget_identity.md',
    'polynomial_reduction_theorem.md': 'docs/polynomial_reduction_theorem.md',
    'polynomial_wall_analysis.md': 'docs/polynomial_wall_analysis.md',
    'K10_computation_results.md': 'docs/K10_computation_results.md',
    'K10_addendum_discriminant_primes.md': 'docs/K10_addendum_discriminant_primes.md',
    'K12_galois_wreath_product.md': 'docs/K12_galois_wreath_product.md',
    'spectral_motive_motivic_status.md': 'docs/spectral_motive_motivic_status.md',
    'M_E_Z_triangle.md': 'docs/M_E_Z_triangle.md',
    'K14_computation_results.md': 'docs/K14_computation_results.md',
    'cfsg_spectral_extractors.md': 'docs/cfsg_spectral_extractors.md',
    'representation_theory_landscape.md': 'docs/representation_theory_landscape.md',
    'continuum_limit_analysis.md': 'docs/continuum_limit_analysis.md',
    'K16_computation_results.md': 'docs/K16_computation_results.md',
    'embedding_dependence_results.md': 'docs/embedding_dependence_results.md',
    'congruence_law_proof.md': 'docs/congruence_law_proof.md',
    'R_oscillation_analysis.md': 'docs/R_oscillation_analysis.md',
    'six_level_chain_update.md': 'docs/six_level_chain_update.md',
    'spectral_function_verdict.md': 'docs/spectral_function_verdict.md',
    'period_structure_terminal.md': 'docs/period_structure_terminal.md',
    'order_one_reduction_chain.md': 'docs/order_one_reduction_chain.md',
    'time_emergence_t_derivation.md': 'docs/time_emergence_t_derivation.md',
    'discriminant_atlas.md': 'docs/discriminant_atlas.md',
    'natural_machine_framework.md': 'docs/natural_machine_framework.md',
    'fixing_lambda_results.md': 'docs/fixing_lambda_results.md',
    'entanglement_ramifications.md': 'docs/entanglement_ramifications.md',
    'c1_investigation_summary.md': 'docs/c1_investigation_summary.md',
    'rho2_hypercharge_results.md': 'docs/rho2_hypercharge_results.md',
    'two_thirds_verdict.md': 'docs/two_thirds_verdict.md',
    'algebraic_analysis_higher_spin.md': 'docs/algebraic_analysis_higher_spin.md',
    'conformal_bootstrap_K4.md': 'docs/conformal_bootstrap_K4.md',
    'gamma3_analytic_evaluation.md': 'docs/gamma3_analytic_evaluation.md',
    'pade_borel_results.md': 'docs/pade_borel_results.md',
    'cT_measurement_spec.md': 'docs/cT_measurement_spec.md',
    'k4d_cosmological_cycling_v1.md': 'docs/k4d_cosmological_cycling_v1.md',
    'D1_manybody_results.md': 'docs/D1_manybody_results.md',
    'K4D_Part2_Formal_Outline.md': 'docs/K4D_Part2_Formal_Outline.md',
    'k5_shadow_graph.md': 'docs/k5_shadow_graph.md',
}


def main():
    if len(sys.argv) > 1:
        src_dir = Path(sys.argv[1])
    else:
        src_path = input("Path to project files directory: ").strip()
        src_dir = Path(src_path)
    
    if not src_dir.exists():
        print(f"ERROR: {src_dir} does not exist")
        sys.exit(1)
    
    # Create directories
    (ROOT / "papers").mkdir(exist_ok=True)
    (ROOT / "docs").mkdir(exist_ok=True)
    (ROOT / "code").mkdir(exist_ok=True)
    
    copied = 0
    missing = 0
    
    for src_name, dest_rel in FILE_MAP.items():
        src = src_dir / src_name
        dest = ROOT / dest_rel
        
        if src.exists():
            dest.parent.mkdir(exist_ok=True)
            shutil.copy2(src, dest)
            print(f"  ✓ {src_name} → {dest_rel}")
            copied += 1
        else:
            # Also check for .docx versions
            docx_name = src_name.replace('.md', '.docx')
            docx_src = src_dir / docx_name
            if docx_src.exists():
                docx_dest = ROOT / dest_rel.replace('.md', '.docx')
                shutil.copy2(docx_src, docx_dest)
                print(f"  ✓ {docx_name} → {dest_rel.replace('.md', '.docx')} (docx)")
                copied += 1
            else:
                print(f"  — {src_name} not found")
                missing += 1
    
    print(f"\nCopied {copied} files, {missing} not found.")
    print(f"\nNext steps:")
    print(f"  1. python scripts/add_frontmatter.py --apply")
    print(f"  2. python scripts/build_site.py")
    print(f"  3. python scripts/dev.py  (preview at localhost:8000)")
    print(f"  4. git init && git add -A && git commit -m 'Initial commit'")
    print(f"  5. Push to GitHub, enable Pages in Settings")


if __name__ == "__main__":
    main()
