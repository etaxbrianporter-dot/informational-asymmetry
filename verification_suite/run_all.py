#!/usr/bin/env python3
"""
run_all.py — Master Verification Suite Runner
==============================================
Informational Asymmetry Program
Brian Porter, February 2026

Executes all test scripts and produces a summary report.
Each test is independent and can also be run standalone.

Usage:
    python run_all.py          # Run all tests
    python run_all.py --quick  # Run only tests 01-09 (fastest subset)
"""
import os
import sys
import subprocess
import time
import re
from datetime import datetime

def run_test(test_file, test_dir):
    """Run a single test script and capture results."""
    start = time.time()
    try:
        result = subprocess.run(
            [sys.executable, test_file],
            cwd=test_dir,
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout per test
        )
        elapsed = time.time() - start
        output = result.stdout + result.stderr
        
        # Parse PASS/FAIL counts from output
        pass_match = re.search(r'(\d+) passed', output)
        fail_match = re.search(r'(\d+) failed', output)
        
        n_pass = int(pass_match.group(1)) if pass_match else 0
        n_fail = int(fail_match.group(1)) if fail_match else 0
        
        return {
            'file': test_file,
            'passed': n_pass,
            'failed': n_fail,
            'time': elapsed,
            'output': output,
            'returncode': result.returncode,
            'error': None
        }
    except subprocess.TimeoutExpired:
        return {
            'file': test_file,
            'passed': 0,
            'failed': 1,
            'time': 300,
            'output': '',
            'returncode': -1,
            'error': 'TIMEOUT (300s)'
        }
    except Exception as e:
        return {
            'file': test_file,
            'passed': 0,
            'failed': 1,
            'time': 0,
            'output': '',
            'returncode': -1,
            'error': str(e)
        }


def main():
    quick = '--quick' in sys.argv
    verbose = '--verbose' in sys.argv or '-v' in sys.argv
    
    # Find test directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    test_dir = os.path.join(script_dir, 'tests')
    
    if not os.path.isdir(test_dir):
        print(f"ERROR: Test directory not found: {test_dir}")
        sys.exit(1)
    
    # Discover test files
    test_files = sorted([
        f for f in os.listdir(test_dir) 
        if f.startswith('test_') and f.endswith('.py')
    ])
    
    if quick:
        test_files = [f for f in test_files if int(f.split('_')[1]) <= 9]
    
    # Header
    print("=" * 72)
    print("  INFORMATIONAL ASYMMETRY PROGRAM: VERIFICATION REPORT")
    print("=" * 72)
    print(f"  Date:   {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Python: {sys.version.split()[0]}")
    try:
        import numpy; print(f"  NumPy:  {numpy.__version__}")
    except: print("  NumPy:  NOT FOUND")
    try:
        import scipy; print(f"  SciPy:  {scipy.__version__}")
    except: print("  SciPy:  NOT FOUND")
    print(f"  Tests:  {len(test_files)} {'(quick mode)' if quick else ''}")
    print("=" * 72)
    print()
    
    # Run tests
    results = []
    total_pass = 0
    total_fail = 0
    total_time = 0
    
    for tf in test_files:
        # Extract test name from filename
        parts = tf.replace('.py', '').split('_', 2)
        test_num = parts[1]
        test_name = parts[2].replace('_', ' ').title() if len(parts) > 2 else ''
        
        sys.stdout.write(f"  Test {test_num}: {test_name:<45s} ")
        sys.stdout.flush()
        
        r = run_test(tf, test_dir)
        results.append(r)
        
        total_pass += r['passed']
        total_fail += r['failed']
        total_time += r['time']
        
        if r['error']:
            status = f"ERROR ({r['error']})"
        elif r['failed'] > 0:
            status = f"FAIL ({r['passed']}/{r['passed']+r['failed']})"
        else:
            status = f"PASS ({r['passed']} claims)"
        
        print(f"{status:>25s}  [{r['time']:.1f}s]")
        
        if verbose and (r['failed'] > 0 or r['error']):
            print()
            for line in r['output'].split('\n'):
                if '✗' in line or 'Error' in line or 'Traceback' in line:
                    print(f"    {line}")
            print()
    
    # Summary
    print()
    print("=" * 72)
    if total_fail == 0:
        print(f"  ALL TESTS PASSED: {total_pass} claims verified, 0 failures")
    else:
        print(f"  RESULT: {total_pass} passed, {total_fail} FAILED")
        print()
        print("  Failed tests:")
        for r in results:
            if r['failed'] > 0 or r['error']:
                fail_count = r['failed']
                msg = r['error'] if r['error'] else f"{fail_count} failures"
                print(f"    {r['file']}: {msg}")
    
    print(f"  Total time: {total_time:.1f}s")
    print("=" * 72)
    
    return 1 if total_fail > 0 else 0


if __name__ == '__main__':
    sys.exit(main())
