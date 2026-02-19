#!/usr/bin/env python3
"""
build_docs_manifest.py — Auto-generate docs.json from docs/*.html

Each HTML doc should include these <meta> tags in its <head>:

    <meta name="doc-category" content="Core Framework">
    <meta name="doc-description" content="One-paragraph claim through five-axiom chain">
    <meta name="doc-date" content="Feb 2026">
    <meta name="doc-sort" content="10">          <!-- optional: controls order within category -->

The <title> tag becomes the display title.

Category sort order is defined in CATEGORY_ORDER below.
Files without meta tags go into "Uncategorized" at the bottom.

Usage:
    python build_docs_manifest.py              # scans ./docs/, writes ./docs.json
    python build_docs_manifest.py path/to/docs  # custom docs dir
"""

import os
import re
import sys
import json
from pathlib import Path

# ── Category display order ──
# Categories not in this list appear at the end, alphabetically.
CATEGORY_ORDER = [
    "Core Framework",
    "K\u2086 Higgs Sector",
    "K\u2088 Fermion Sector",
    "Gauge Couplings",
    "Mathematical Structure",
    "Gravitational Sector",
    "Killed Approaches",
]


def extract_meta(html_text, name):
    """Extract content from <meta name="..." content="...">"""
    pattern = rf'<meta\s+name="{re.escape(name)}"\s+content="([^"]*)"'
    m = re.search(pattern, html_text, re.IGNORECASE)
    if m:
        return m.group(1).strip()
    # Also try reversed attribute order: content before name
    pattern2 = rf'<meta\s+content="([^"]*)"\s+name="{re.escape(name)}"'
    m2 = re.search(pattern2, html_text, re.IGNORECASE)
    return m2.group(1).strip() if m2 else None


def extract_title(html_text):
    """Extract text from <title>...</title>"""
    m = re.search(r'<title>([^<]+)</title>', html_text, re.IGNORECASE)
    return m.group(1).strip() if m else None


def scan_docs(docs_dir):
    """Scan docs directory and return list of doc entries."""
    docs_path = Path(docs_dir)
    if not docs_path.is_dir():
        print(f"Error: {docs_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    entries = []
    for html_file in sorted(docs_path.glob("*.html")):
        if html_file.name.startswith("_"):
            continue  # skip templates and hidden files
        text = html_file.read_text(encoding="utf-8", errors="replace")

        title = extract_meta(text, "doc-title") or extract_title(text)
        category = extract_meta(text, "doc-category") or "Uncategorized"
        desc = extract_meta(text, "doc-description") or ""
        date = extract_meta(text, "doc-date") or ""
        sort_key = extract_meta(text, "doc-sort") or "50"

        if not title:
            # Fallback: humanize filename
            title = html_file.stem.replace("_", " ").replace("-", " ").title()

        try:
            sort_val = int(sort_key)
        except ValueError:
            sort_val = 50

        entries.append({
            "file": html_file.name,
            "href": f"docs/{html_file.name}",
            "title": title,
            "category": category,
            "desc": desc,
            "date": date,
            "sort": sort_val,
        })

    return entries


def group_by_category(entries):
    """Group entries by category, respecting CATEGORY_ORDER."""
    groups = {}
    for e in entries:
        cat = e["category"]
        groups.setdefault(cat, []).append(e)

    # Sort within each category
    for cat in groups:
        groups[cat].sort(key=lambda x: (x["sort"], x["title"]))

    # Order categories
    ordered = []
    seen = set()
    for cat in CATEGORY_ORDER:
        if cat in groups:
            ordered.append({
                "category": cat,
                "docs": [
                    {"date": d["date"], "href": d["href"], "title": d["title"], "desc": d["desc"]}
                    for d in groups[cat]
                ]
            })
            seen.add(cat)

    # Remaining categories alphabetically
    for cat in sorted(groups.keys()):
        if cat not in seen:
            ordered.append({
                "category": cat,
                "docs": [
                    {"date": d["date"], "href": d["href"], "title": d["title"], "desc": d["desc"]}
                    for d in groups[cat]
                ]
            })

    return ordered


def main():
    docs_dir = sys.argv[1] if len(sys.argv) > 1 else "docs"
    output = sys.argv[2] if len(sys.argv) > 2 else "docs.json"

    entries = scan_docs(docs_dir)
    manifest = group_by_category(entries)

    total = sum(len(g["docs"]) for g in manifest)
    cats = len(manifest)

    with open(output, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)

    print(f"docs.json: {total} documents in {cats} categories")
    for g in manifest:
        print(f"  {g['category']}: {len(g['docs'])} docs")


if __name__ == "__main__":
    main()
