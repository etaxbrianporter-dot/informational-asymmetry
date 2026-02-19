#!/usr/bin/env python3
"""
dev.py — Local development server with auto-rebuild.
=====================================================
Watches docs/ for changes, rebuilds, serves at localhost:8000.

Usage:
    python scripts/dev.py
"""
import os
import sys
import time
import http.server
import socketserver
import threading
import subprocess
from pathlib import Path

ROOT = Path(__file__).parent.parent
BUILD_DIR = ROOT / "site_build"
PORT = 8000


def build():
    """Run the build script."""
    print("\n  Rebuilding...", flush=True)
    result = subprocess.run(
        [sys.executable, "scripts/build_site.py"],
        cwd=ROOT, capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"  Build error:\n{result.stderr}")
    else:
        print("  Build complete.", flush=True)


def watch():
    """Watch for file changes and rebuild."""
    last_build = time.time()
    watch_dirs = [ROOT / "docs", ROOT / "templates", ROOT / "papers"]
    watch_files = [ROOT / "index.html"]
    
    def get_mtimes():
        mtimes = {}
        for d in watch_dirs:
            if d.exists():
                for f in d.rglob("*"):
                    if f.is_file():
                        mtimes[str(f)] = f.stat().st_mtime
        for f in watch_files:
            if f.exists():
                mtimes[str(f)] = f.stat().st_mtime
        return mtimes
    
    prev_mtimes = get_mtimes()
    
    while True:
        time.sleep(1)
        curr_mtimes = get_mtimes()
        if curr_mtimes != prev_mtimes:
            changed = set(curr_mtimes.keys()) ^ set(prev_mtimes.keys())
            for k in set(curr_mtimes.keys()) & set(prev_mtimes.keys()):
                if curr_mtimes[k] != prev_mtimes[k]:
                    changed.add(k)
            if changed:
                names = [os.path.basename(f) for f in list(changed)[:3]]
                print(f"\n  Changed: {', '.join(names)}")
                build()
                prev_mtimes = get_mtimes()


def serve():
    """Serve the build directory."""
    os.chdir(BUILD_DIR)
    handler = http.server.SimpleHTTPRequestHandler
    handler.extensions_map.update({'.html': 'text/html', '.json': 'application/json'})
    
    with socketserver.TCPServer(("", PORT), handler) as httpd:
        print(f"\n  Serving at http://localhost:{PORT}")
        print(f"  Watching for changes... (Ctrl+C to stop)\n")
        httpd.serve_forever()


def main():
    # Initial build
    build()
    
    # Copy index.html to build dir
    import shutil
    src = ROOT / "index.html"
    if src.exists():
        BUILD_DIR.mkdir(exist_ok=True)
        shutil.copy2(src, BUILD_DIR / "index.html")
    
    # Start watcher in background
    watcher = threading.Thread(target=watch, daemon=True)
    watcher.start()
    
    # Serve (blocks)
    try:
        serve()
    except KeyboardInterrupt:
        print("\n  Stopped.")


if __name__ == "__main__":
    main()
