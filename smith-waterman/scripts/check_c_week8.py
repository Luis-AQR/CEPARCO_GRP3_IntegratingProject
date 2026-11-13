"""
check_c_week8.py - Week 8: compare the sequential C program with the independent
pure-Python reference (scripts/check_expected.py) on score, best cell, matrix
checksum and, for the small tests, the whole score matrix.

Works on Windows, macOS and Linux. Build the C program first (from the repo root):
    build.bat                                   (Windows)

Usage (from the repo root):
    python scripts\\check_c_week8.py             tests + bench_1000
    python scripts\\check_c_week8.py --all       also bench_5000 (slower: pure Python)
    python scripts\\check_c_week8.py --exe bin\\sw_seq.exe    use a specific executable
"""
import argparse
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_expected import read_fasta_pair, smith_waterman  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DATA = os.path.join(ROOT, "data")


def find_exe(explicit=None):
    """Locate the compiled program: bin/sw_seq.exe (Windows) or bin/sw_seq."""
    candidates = [explicit] if explicit else [
        os.path.join(ROOT, "bin", "sw_seq.exe"),
        os.path.join(ROOT, "bin", "sw_seq"),
    ]
    for c in candidates:
        if c and os.path.isfile(c):
            return os.path.abspath(c)
    sys.exit("error: compiled program not found (looked for bin\\sw_seq.exe). "
             "Run build.bat first, or pass --exe PATH.")


def run_c(exe, path):
    out = subprocess.run([exe, path], capture_output=True, text=True, check=True).stdout
    score = int(re.search(r"Max score:\s*(-?\d+)", out).group(1))
    m = re.search(r"Max cell\s*:\s*H\((\d+), (\d+)\)", out)
    cell = (int(m.group(1)), int(m.group(2))) if m else (0, 0)
    checksum = int(re.search(r"Checksum\s*:\s*(\d+)", out).group(1))
    rows = []
    lines = out.splitlines()
    for k, line in enumerate(lines):
        if line.startswith("Score matrix"):
            for row in lines[k + 4:]:
                if not row.strip():
                    break
                toks = row.split()
                vals = [t for t in toks[1:] if re.fullmatch(r"-?\d+", t)]
                rows.append([int(v) for v in vals])
            break
    return score, cell, checksum, rows


def main():
    ap = argparse.ArgumentParser(description="Compare sw_seq with the Python reference.")
    ap.add_argument("--all", action="store_true", help="also check bench_5000")
    ap.add_argument("--exe", help="path to the compiled sw_seq executable")
    args = ap.parse_args()

    exe = find_exe(args.exe)
    names = ["test_paper", "test_identical", "test_nomatch", "test_gap", "bench_1000"]
    if args.all:
        names.append("bench_5000")

    print(f"Using {exe}\n")
    bad = 0
    for n in names:
        p = os.path.join(DATA, n + ".fasta")
        a, b = read_fasta_pair(p)
        H, score, cell, *_ = smith_waterman(a, b)
        checksum = sum(sum(r) for r in H)
        c_score, c_cell, c_sum, c_rows = run_c(exe, p)
        ok = (c_score == score and c_cell == cell and c_sum == checksum)
        if c_rows and len(a) <= 40 and len(b) <= 40:
            ok = ok and c_rows == H
        print(f"  {'ok  ' if ok else 'FAIL'} {n:16s} C: {c_score} at {c_cell} sum={c_sum}"
              f"   Python: {score} at {cell} sum={checksum}"
              + ("   (matrix compared)" if c_rows else ""))
        bad += (not ok)
    print("\nAll C results match the reference." if not bad else f"\n{bad} mismatch(es).")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
