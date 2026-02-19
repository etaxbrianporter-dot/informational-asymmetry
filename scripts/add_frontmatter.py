#!/usr/bin/env python3
"""
add_frontmatter.py — Add YAML frontmatter to markdown docs.
=============================================================
Scans docs/*.md, adds frontmatter to files that don't have it.
Uses filename patterns to auto-assign sections.

Usage:
    python scripts/add_frontmatter.py              # Preview changes
    python scripts/add_frontmatter.py --apply      # Write changes
"""
import sys
import re
from pathlib import Path

ROOT = Path(__file__).parent.parent
DOCS_DIR = ROOT / "docs"

# Section assignment rules (pattern → section)
SECTION_RULES = [
    (r'k4|K4|pfaffian|lorentzian|hessian|signature', 'K₄ Spacetime'),
    (r'k6|K6|higgs|normalization|gram|a2|a4', 'K₆ Higgs Sector'),
    (r'k8|K8|fermion|yukawa|generation|heawood', 'K₈ Fermion Sector'),
    (r'k10|K10|k12|K12|k14|K14|k16|K16', 'Higher Levels'),
    (r'birefringence|CMB|torsion|litebird', 'Cosmic Birefringence'),
    (r'johnson|budget|spectral_budget|polynomial_reduction|kneser', 'Mathematical Structure'),
    (r'galois|discriminant|number_field|cyclotomic', 'Arithmetic Structure'),
    (r'gravity|einstein|gamma3|QMC|observer|dimension', 'Gravitational Sector'),
    (r'termination|killed|wall|verdict|bug', 'Killed Approaches'),
    (r'executive|progress|classification|unified|summary', 'Program Overview'),
    (r'entangle|emergence|light|time|arrow', 'Emergence'),
]

def guess_section(filename, content):
    """Guess the section from filename and content."""
    text = filename + " " + content[:500]
    for pattern, section in SECTION_RULES:
        if re.search(pattern, text, re.IGNORECASE):
            return section
    return ""

def extract_title(content):
    """Extract title from first # heading."""
    match = re.match(r'^#\s+(.+)', content.strip(), re.MULTILINE)
    if match:
        return match.group(1).strip()
    return None

def has_frontmatter(content):
    """Check if file already has YAML frontmatter."""
    return content.strip().startswith("---")

def add_frontmatter(filepath, apply=False):
    """Add frontmatter to a single file."""
    content = filepath.read_text(encoding='utf-8', errors='replace')
    
    if has_frontmatter(content):
        return None, "already has frontmatter"
    
    title = extract_title(content) or filepath.stem.replace('_', ' ').title()
    section = guess_section(filepath.stem, content)
    
    frontmatter = f"""---
title: {title}
section: {section}
status: active
---

"""
    new_content = frontmatter + content
    
    if apply:
        filepath.write_text(new_content, encoding='utf-8')
        return title, f"added (section: {section})"
    else:
        return title, f"would add (section: {section})"


def main():
    apply = "--apply" in sys.argv
    
    if not DOCS_DIR.exists():
        print(f"No docs/ directory found. Create it and add .md files.")
        return
    
    md_files = sorted(DOCS_DIR.glob("*.md"))
    print(f"Found {len(md_files)} markdown files in docs/\n")
    
    for f in md_files:
        title, status = add_frontmatter(f, apply=apply)
        if title:
            print(f"  {'✓' if apply else '·'} {f.name}: {status}")
        else:
            print(f"  — {f.name}: {status}")
    
    if not apply:
        print(f"\nDry run. Use --apply to write changes.")


if __name__ == "__main__":
    main()
