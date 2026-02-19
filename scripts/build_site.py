#!/usr/bin/env python3
"""
build_site.py — Convert markdown working documents to styled HTML pages.
========================================================================
Reads docs/*.md, applies template, outputs to site_build/docs/*.html.
Generates docs index with git-derived timestamps.

Usage:
    python scripts/build_site.py              # Full build
    python scripts/build_site.py --watch      # Rebuild on file changes (local dev)
"""
import os
import sys
import re
import json
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

try:
    import markdown
    from markdown.extensions.toc import TocExtension
    from markdown.extensions.tables import TableExtension
    from markdown.extensions.fenced_code import FencedCodeExtension
    from markdown.extensions.codehilite import CodeHiliteExtension
    HAS_MARKDOWN = True
except ImportError:
    HAS_MARKDOWN = False
    print("WARNING: 'markdown' package not installed. Using basic converter.")

# ── Paths ──
ROOT = Path(__file__).parent.parent
DOCS_DIR = ROOT / "docs"
TEMPLATE_DIR = ROOT / "templates"
BUILD_DIR = ROOT / "site_build"
DOCS_BUILD = BUILD_DIR / "docs"


def git_last_modified(filepath):
    """Get the last git commit date for a file."""
    try:
        result = subprocess.run(
            ["git", "log", "-1", "--format=%aI", "--", str(filepath)],
            capture_output=True, text=True, cwd=ROOT
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()[:10]  # YYYY-MM-DD
    except FileNotFoundError:
        pass
    # Fallback to file mtime
    mtime = os.path.getmtime(filepath)
    return datetime.fromtimestamp(mtime).strftime("%Y-%m-%d")


def git_first_added(filepath):
    """Get the first git commit date for a file."""
    try:
        result = subprocess.run(
            ["git", "log", "--diff-filter=A", "--format=%aI", "--", str(filepath)],
            capture_output=True, text=True, cwd=ROOT
        )
        if result.returncode == 0 and result.stdout.strip():
            dates = result.stdout.strip().split('\n')
            return dates[-1][:10]
    except FileNotFoundError:
        pass
    return git_last_modified(filepath)


def parse_frontmatter(content):
    """
    Extract YAML-like frontmatter from markdown.
    ---
    title: Something
    section: K6 Higgs Sector
    status: active
    ---
    """
    meta = {}
    body = content
    
    if content.startswith("---"):
        parts = content.split("---", 2)
        if len(parts) >= 3:
            for line in parts[1].strip().split("\n"):
                if ":" in line:
                    key, val = line.split(":", 1)
                    meta[key.strip()] = val.strip()
            body = parts[2]
    
    # If no frontmatter title, extract from first # heading
    if "title" not in meta:
        h1_match = re.match(r'^#\s+(.+)', body.strip())
        if h1_match:
            meta["title"] = h1_match.group(1)
    
    if "title" not in meta:
        meta["title"] = "Untitled"
    
    return meta, body


def convert_markdown(text):
    """Convert markdown to HTML."""
    if HAS_MARKDOWN:
        md = markdown.Markdown(extensions=[
            TocExtension(permalink=False, toc_depth=3),
            TableExtension(),
            FencedCodeExtension(),
            CodeHiliteExtension(css_class='highlight', guess_lang=False),
            'markdown.extensions.nl2br',
        ])
        return md.convert(text)
    else:
        # Basic fallback: escape HTML, handle headers, code blocks, paragraphs
        import html as html_mod
        text = html_mod.escape(text)
        # Headers
        text = re.sub(r'^### (.+)$', r'<h3>\1</h3>', text, flags=re.MULTILINE)
        text = re.sub(r'^## (.+)$', r'<h2>\1</h2>', text, flags=re.MULTILINE)
        text = re.sub(r'^# (.+)$', r'<h1>\1</h1>', text, flags=re.MULTILINE)
        # Bold/italic
        text = re.sub(r'\*\*(.+?)\*\*', r'<strong>\1</strong>', text)
        text = re.sub(r'\*(.+?)\*', r'<em>\1</em>', text)
        # Code
        text = re.sub(r'`(.+?)`', r'<code>\1</code>', text)
        # Paragraphs
        text = re.sub(r'\n\n', '</p><p>', text)
        text = f'<p>{text}</p>'
        return text


def load_template():
    """Load the HTML template for doc pages."""
    template_path = TEMPLATE_DIR / "doc.html"
    if template_path.exists():
        return template_path.read_text()
    
    # Inline default template
    return '''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{{title}} — Informational Asymmetry Program</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=EB+Garamond:ital,wght@0,400;0,500;0,600;0,700;1,400;1,500&family=JetBrains+Mono:wght@400;500;600&family=Source+Sans+3:wght@300;400;500;600;700&display=swap" rel="stylesheet">
<style>
:root {
  --ink: #1a1a1a;
  --paper: #fafaf7;
  --accent: #8b2500;
  --mid: #6b6b6b;
  --rule: #d4d0c8;
  --code-bg: #f0ede6;
  --serif: 'EB Garamond', Georgia, serif;
  --sans: 'Source Sans 3', system-ui, sans-serif;
  --mono: 'JetBrains Mono', 'Consolas', monospace;
}
* { margin: 0; padding: 0; box-sizing: border-box; }
body {
  font-family: var(--sans);
  color: var(--ink);
  background: var(--paper);
  line-height: 1.7;
  -webkit-font-smoothing: antialiased;
}
.page-header {
  border-bottom: 1px solid var(--rule);
  padding: 1.5rem 0;
}
.page-header .container {
  display: flex;
  justify-content: space-between;
  align-items: baseline;
  flex-wrap: wrap;
  gap: 1rem;
}
.site-title {
  font-family: var(--serif);
  font-size: 1.2rem;
  font-weight: 600;
  color: var(--ink);
  text-decoration: none;
}
.site-title:hover { color: var(--accent); }
nav a {
  font-size: 0.8rem;
  font-weight: 500;
  color: var(--mid);
  text-decoration: none;
  letter-spacing: 0.02em;
  text-transform: uppercase;
  margin-left: 1.5rem;
}
nav a:hover { color: var(--ink); }

.container { max-width: 960px; margin: 0 auto; padding: 0 2rem; }

.doc-header {
  padding: 3rem 0 2rem;
  border-bottom: 1px solid var(--rule);
}
.doc-header h1 {
  font-family: var(--serif);
  font-size: 2.2rem;
  font-weight: 600;
  line-height: 1.2;
  letter-spacing: -0.02em;
  margin-bottom: 0.75rem;
}
.doc-meta {
  font-size: 0.85rem;
  color: var(--mid);
}
.doc-meta .badge {
  display: inline-block;
  font-size: 0.7rem;
  font-weight: 600;
  text-transform: uppercase;
  letter-spacing: 0.04em;
  padding: 0.2em 0.6em;
  border-radius: 3px;
  background: #e8f5e9;
  color: #2d6a2e;
  margin-left: 0.5rem;
}

article {
  padding: 2rem 0 4rem;
  max-width: 720px;
}
article h1 { font-family: var(--serif); font-size: 1.8rem; font-weight: 600; margin: 2rem 0 1rem; }
article h2 { font-family: var(--serif); font-size: 1.5rem; font-weight: 600; margin: 2rem 0 0.75rem; border-bottom: 1px solid var(--rule); padding-bottom: 0.3rem; }
article h3 { font-family: var(--sans); font-size: 1.05rem; font-weight: 600; margin: 1.5rem 0 0.5rem; }
article p { margin-bottom: 1rem; }
article ul, article ol { margin: 0 0 1rem 1.5rem; }
article li { margin-bottom: 0.3rem; }
article blockquote {
  border-left: 3px solid var(--accent);
  padding: 0.5rem 1.5rem;
  margin: 1rem 0;
  color: var(--mid);
  background: #f5f3ee;
}
article code {
  font-family: var(--mono);
  font-size: 0.85em;
  background: var(--code-bg);
  padding: 0.15em 0.4em;
  border-radius: 3px;
}
article pre {
  background: #2b2b2b;
  color: #e0e0e0;
  font-family: var(--mono);
  font-size: 0.8rem;
  line-height: 1.6;
  padding: 1.5rem 2rem;
  border-radius: 6px;
  overflow-x: auto;
  margin: 1.5rem 0;
}
article pre code {
  background: none;
  padding: 0;
  color: inherit;
}
article table {
  border-collapse: collapse;
  width: 100%;
  margin: 1rem 0;
  font-size: 0.9rem;
}
article th {
  text-align: left;
  font-weight: 600;
  padding: 0.5rem 0.75rem;
  border-bottom: 2px solid var(--rule);
  font-size: 0.8rem;
  text-transform: uppercase;
  letter-spacing: 0.04em;
  color: var(--mid);
}
article td {
  padding: 0.5rem 0.75rem;
  border-bottom: 1px solid #eae7e0;
}
article strong { font-weight: 600; }
article em { font-style: italic; }
article hr {
  border: none;
  border-top: 1px solid var(--rule);
  margin: 2rem 0;
}
article img {
  max-width: 100%;
  border-radius: 4px;
  margin: 1rem 0;
}

footer {
  border-top: 1px solid var(--rule);
  padding: 1.5rem 0;
  font-size: 0.8rem;
  color: var(--mid);
}
footer a { color: var(--mid); }
footer a:hover { color: var(--accent); }
</style>
</head>
<body>

<div class="page-header">
  <div class="container">
    <a href="../index.html" class="site-title">Informational Asymmetry Program</a>
    <nav>
      <a href="../index.html#framework">Framework</a>
      <a href="../index.html#papers">Papers</a>
      <a href="../index.html#docs">Docs</a>
      <a href="../index.html#verification">Verification</a>
    </nav>
  </div>
</div>

<div class="container">
  <div class="doc-header">
    <h1>{{title}}</h1>
    <div class="doc-meta">
      {{meta_line}}
    </div>
  </div>
  <article>
    {{content}}
  </article>
</div>

<footer>
  <div class="container">
    Informational Asymmetry Program · Brian Porter · 2026 · 
    <a href="../index.html">Home</a> · 
    <a href="https://github.com/bporter/informational-asymmetry">Source</a>
  </div>
</footer>

</body>
</html>'''


def build_doc(md_path, template):
    """Convert a single markdown file to HTML."""
    content = md_path.read_text(encoding='utf-8', errors='replace')
    meta, body = parse_frontmatter(content)
    
    html_body = convert_markdown(body)
    
    # Build metadata line
    modified = git_last_modified(md_path)
    created = git_first_added(md_path)
    section = meta.get("section", "")
    status = meta.get("status", "")
    
    meta_parts = [f"Brian Porter"]
    if created == modified:
        meta_parts.append(f"· {modified}")
    else:
        meta_parts.append(f"· Created {created} · Updated {modified}")
    if section:
        meta_parts.append(f"· {section}")
    if status:
        meta_parts.append(f'<span class="badge">{status}</span>')
    
    meta_line = " ".join(meta_parts)
    
    # Apply template
    page = template.replace("{{title}}", meta.get("title", "Untitled"))
    page = page.replace("{{content}}", html_body)
    page = page.replace("{{meta_line}}", meta_line)
    
    return page, meta, modified


def build_docs_index(doc_entries):
    """Generate a JSON manifest of all docs for the main site to consume."""
    return json.dumps(doc_entries, indent=2)


def copy_docx_and_pdf(src_dir, dest_dir):
    """Copy non-markdown docs (docx, pdf, py) to build."""
    for ext in ['*.docx', '*.pdf', '*.py', '*.tex']:
        for f in src_dir.glob(ext):
            shutil.copy2(f, dest_dir / f.name)


def main():
    print("=" * 60)
    print("  Building Informational Asymmetry Program site")
    print("=" * 60)
    
    # Create build directories
    BUILD_DIR.mkdir(exist_ok=True)
    DOCS_BUILD.mkdir(exist_ok=True)
    
    template = load_template()
    
    # Build all markdown docs
    doc_entries = []
    md_files = sorted(DOCS_DIR.glob("*.md")) if DOCS_DIR.exists() else []
    
    print(f"\n  Found {len(md_files)} markdown documents in docs/")
    
    for md_path in md_files:
        try:
            html, meta, modified = build_doc(md_path, template)
            out_name = md_path.stem + ".html"
            out_path = DOCS_BUILD / out_name
            out_path.write_text(html, encoding='utf-8')
            
            doc_entries.append({
                "file": out_name,
                "title": meta.get("title", md_path.stem),
                "section": meta.get("section", ""),
                "status": meta.get("status", ""),
                "modified": modified,
                "source": md_path.name,
            })
            print(f"  ✓ {md_path.name} → {out_name}")
        except Exception as e:
            print(f"  ✗ {md_path.name}: {e}")
    
    # Also convert any .md in root (like README)
    root_md = ROOT / "README.md"
    if root_md.exists():
        try:
            html, meta, modified = build_doc(root_md, template)
            (BUILD_DIR / "about.html").write_text(html, encoding='utf-8')
            print(f"  ✓ README.md → about.html")
        except Exception as e:
            print(f"  ✗ README.md: {e}")
    
    # Write docs manifest
    manifest_path = BUILD_DIR / "docs_manifest.json"
    manifest_path.write_text(build_docs_index(doc_entries), encoding='utf-8')
    print(f"\n  Wrote manifest: {len(doc_entries)} entries → docs_manifest.json")
    
    # Copy non-markdown assets
    if DOCS_DIR.exists():
        copy_docx_and_pdf(DOCS_DIR, DOCS_BUILD)
    
    # Copy papers
    papers_dir = ROOT / "papers"
    papers_build = BUILD_DIR / "papers"
    if papers_dir.exists():
        papers_build.mkdir(exist_ok=True)
        for f in papers_dir.iterdir():
            shutil.copy2(f, papers_build / f.name)
        print(f"  Copied papers/")
    
    # Copy code directory
    code_dir = ROOT / "code"
    code_build = BUILD_DIR / "code"
    if code_dir.exists():
        if code_build.exists():
            shutil.rmtree(code_build)
        shutil.copytree(code_dir, code_build)
        print(f"  Copied code/")
    
    print(f"\n  Build complete → {BUILD_DIR}/")
    print("=" * 60)


if __name__ == "__main__":
    main()
