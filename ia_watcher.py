#!/usr/bin/env python3
"""
ia_watcher.py — File watcher for Informational Asymmetry pipeline
=================================================================
Watches ~/Downloads for new files and routes them into the repo:

  *.py                              → code/
  paper1*, paperII*, paperIII*      → papers/
  test_*, verification_*, lib files → verification_suite/
  everything else                   → docs/

After each routed file: rebuilds docs manifest, commits, and pushes.

Usage:
    python ia_watcher.py                 # watch ~/Downloads
    python ia_watcher.py --dry-run       # log routing without moving/pushing
    python ia_watcher.py --downloads /path/to/dir   # custom watch dir

Requires: pip install watchdog
"""
from __future__ import annotations

import os
import re
import sys
import time
import shutil
import logging
import argparse
import subprocess
from pathlib import Path
from datetime import datetime

try:
    from watchdog.observers import Observer
    from watchdog.events import FileSystemEventHandler
except ImportError:
    print("ERROR: watchdog not installed. Run: pip install watchdog")
    sys.exit(1)

REPO_DIR = Path.home() / "scripts" / "informational-asymmetry"
DEFAULT_DOWNLOADS = Path("/mnt/c/Users/brian/Downloads")

IGNORE_PATTERNS = [
    r".*:Zone\.Identifier$",
    r"^\..*",
    r".*\.tmp$",
    r".*\.crdownload$",
    r".*\.part$",
]

PAPER_PATTERNS = [
    r"^paper1[\._\-\s]",
    r"^paper1\.(pdf|tex|docx|md)$",
    r"^paperII[\._\-\s]",
    r"^paperII\.(pdf|tex|docx|md)$",
    r"^paperIII[\._\-\s]",
    r"^paperIII\.(pdf|tex|docx|md)$",
]

VERIFICATION_PATTERNS = [
    r"^test_\d+",
    r"^verification_bridge",
    r"^run_all\.py$",
    r"^(gram|johnson|k4|matching)_tools\.py$",
    r"^PLAN\.md$",
    r"^requirements\.txt$",
]

LOG_FILE = REPO_DIR / "watcher.log"

def setup_logging():
    fmt = "%(asctime)s  %(levelname)-7s  %(message)s"
    datefmt = "%Y-%m-%d %H:%M:%S"
    logging.basicConfig(
        level=logging.INFO, format=fmt, datefmt=datefmt,
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(LOG_FILE, encoding="utf-8"),
        ],
    )

log = logging.getLogger("ia_watcher")

def should_ignore(filename: str) -> bool:
    for pat in IGNORE_PATTERNS:
        if re.match(pat, filename):
            return True
    return False

def classify_file(filename: str) -> tuple[str, Path]:
    name = filename
    lower = filename.lower()
    for pat in PAPER_PATTERNS:
        if re.match(pat, name, re.IGNORECASE):
            return "papers", REPO_DIR / "papers"
    for pat in VERIFICATION_PATTERNS:
        if re.match(pat, name, re.IGNORECASE):
            if name.endswith("_tools.py"):
                return "verification_suite/lib", REPO_DIR / "verification_suite" / "lib"
            if name.startswith("test_"):
                return "verification_suite/tests", REPO_DIR / "verification_suite" / "tests"
            return "verification_suite", REPO_DIR / "verification_suite"
    if lower.endswith(".py"):
        return "code", REPO_DIR / "code"
    return "docs", REPO_DIR / "docs"

def strip_zone_identifier(dest_dir: Path):
    for zf in dest_dir.glob("*:Zone.Identifier"):
        try:
            zf.unlink()
        except OSError:
            pass

def run_cmd(cmd, cwd=None, check=True):
    result = subprocess.run(
        cmd, capture_output=True, text=True,
        cwd=cwd or REPO_DIR, timeout=120
    )
    if check and result.returncode != 0:
        log.error(f"Command failed: {' '.join(cmd)}")
        log.error(f"  stderr: {result.stderr.strip()}")
        return None
    return result.stdout.strip()

def rebuild_manifest():
    log.info("Rebuilding docs manifest...")
    result = run_cmd(
        [sys.executable, str(REPO_DIR / "build_docs_manifest.py")],
        check=False
    )
    if result is not None:
        log.info(f"  {result}")
    else:
        log.warning("  Manifest rebuild had issues (non-fatal)")

def rebuild_site():
    build_script = REPO_DIR / "scripts" / "build_site.py"
    if build_script.exists():
        log.info("Running site build...")
        result = run_cmd(
            [sys.executable, str(build_script)],
            check=False
        )
        if result:
            for line in result.split("\n"):
                if "Build complete" in line or "✓" in line:
                    log.info(f"  {line.strip()}")

def run_verification():
    run_all = REPO_DIR / "verification_suite" / "run_all.py"
    if not run_all.exists():
        log.warning("Verification suite not found, skipping")
        return True
    log.info("Running verification suite...")
    result = subprocess.run(
        [sys.executable, str(run_all)],
        capture_output=True, text=True,
        cwd=REPO_DIR / "verification_suite",
        timeout=300
    )
    output = result.stdout + result.stderr
    if "ALL TESTS PASSED" in output:
        log.info("  ✓ All tests passing")
        return True
    else:
        log.error("  ✗ TESTS FAILED — will commit but flag in message")
        return False

def git_commit_and_push(filename: str, category: str, tests_passed: bool):
    status = run_cmd(["git", "status", "--porcelain"], check=False)
    if not status:
        log.info("  No git changes to commit")
        return
    if tests_passed:
        msg = f"auto: {filename} → {category}"
    else:
        msg = f"auto: {filename} → {category} [TESTS FAILING]"
    log.info(f"  Committing: {msg}")
    run_cmd(["git", "add", "-A"])
    result = run_cmd(["git", "commit", "-m", msg], check=False)
    if result is None:
        log.error("  Commit failed")
        return
    log.info("  Pushing to remote...")
    push_result = run_cmd(["git", "push"], check=False)
    if push_result is not None:
        log.info("  ✓ Pushed successfully")
    else:
        log.error("  ✗ Push failed — will retry on next file event")

class DownloadHandler(FileSystemEventHandler):
    def __init__(self, dry_run=False):
        super().__init__()
        self.dry_run = dry_run
        self._debounce = {}

    def on_created(self, event):
        if event.is_directory:
            return
        self._handle(event.src_path)

    def on_moved(self, event):
        if event.is_directory:
            return
        self._handle(event.dest_path)

    def _handle(self, filepath):
        path = Path(filepath)
        filename = path.name
        if should_ignore(filename):
            return
        now = time.time()
        if filename in self._debounce and (now - self._debounce[filename]) < 5:
            return
        self._debounce[filename] = now
        time.sleep(1)
        if not path.exists():
            return
        category, dest_dir = classify_file(filename)
        dest_path = dest_dir / filename
        log.info(f"")
        log.info(f"{'─' * 50}")
        log.info(f"  New file: {filename}")
        log.info(f"  Route:    → {category}/")
        if self.dry_run:
            log.info(f"  [DRY RUN] Would copy to {dest_path}")
            return
        try:
            dest_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(str(path), str(dest_path))
            log.info(f"  Copied to {dest_path.relative_to(REPO_DIR)}")
        except Exception as e:
            log.error(f"  Copy failed: {e}")
            return
        strip_zone_identifier(dest_dir)
        if category == "docs":
            rebuild_manifest()
        tests_passed = True
        if category.startswith("verification_suite") or category == "code":
            tests_passed = run_verification()
        git_commit_and_push(filename, category, tests_passed)
        log.info(f"  Done.")

def main():
    global REPO_DIR
    parser = argparse.ArgumentParser(
        description="Watch Downloads and route files to informational-asymmetry repo"
    )
    parser.add_argument("--downloads", type=Path, default=DEFAULT_DOWNLOADS,
        help=f"Directory to watch (default: {DEFAULT_DOWNLOADS})")
    parser.add_argument("--dry-run", action="store_true",
        help="Log routing decisions without copying or pushing")
    parser.add_argument("--repo", type=Path, default=REPO_DIR,
        help=f"Path to repo (default: {REPO_DIR})")
    args = parser.parse_args()
    REPO_DIR = args.repo
    setup_logging()
    if not args.downloads.is_dir():
        log.error(f"Watch directory does not exist: {args.downloads}")
        sys.exit(1)
    if not (REPO_DIR / ".git").is_dir():
        log.error(f"Repo not found or not a git repo: {REPO_DIR}")
        sys.exit(1)
    log.info("═" * 50)
    log.info("  Informational Asymmetry — File Watcher")
    log.info("═" * 50)
    log.info(f"  Watching:  {args.downloads}")
    log.info(f"  Repo:      {REPO_DIR}")
    log.info(f"  Dry run:   {args.dry_run}")
    log.info(f"  Routing:")
    log.info(f"    *.py              → code/")
    log.info(f"    paper1/II/III.*   → papers/")
    log.info(f"    test_*/verify*    → verification_suite/")
    log.info(f"    everything else   → docs/")
    log.info("")
    handler = DownloadHandler(dry_run=args.dry_run)
    observer = Observer()
    observer.schedule(handler, str(args.downloads), recursive=False)
    observer.start()
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        log.info("\nShutting down watcher...")
        observer.stop()
    observer.join()

if __name__ == "__main__":
    main()
