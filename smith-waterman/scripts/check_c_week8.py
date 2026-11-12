#!/usr/bin/env python3
"""
check_c_week8.py - Week 8: compare the sequential C program with the independent
pure-Python reference (scripts/check_expected.py) on score, best cell and,
for the small tests, the whole score matrix.

Usage (from smith-waterman/):
    gcc -O2 -Wall -Wextra -o c/sw_seq c/sw_seq.c
    python3 scripts/check_c_week8.py            # tests + bench_1000
    python3 scripts/check_c_week8.py --all      # also bench_5000 (slower, pure Python)
"""
import os, re, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_expected import read_fasta_pair, smith_waterman

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
EXE = os.path.join(ROOT, "c", "sw_seq")
DATA = os.path.join(ROOT, "data")


def run_c(path):
    out = subprocess.run([EXE, path], capture_output=True, text=True, check=True).stdout
    score = int(re.search(r"Max score:\s*(-?\d+)", out).group(1))
    m = re.search(r"Max cell\s*:\s*H\((\d+), (\d+)\)", out)
    cell = (int(m.group(1)), int(m.group(2))) if m else (0, 0)
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
    return score, cell, rows


def main():
    names = ["test_paper", "test_identical", "test_nomatch", "test_gap", "bench_1000"]
    if "--all" in sys.argv:
        names.append("bench_5000")
    bad = 0
    for n in names:
        p = os.path.join(DATA, n + ".fasta")
        a, b = read_fasta_pair(p)
        H, score, cell, *_ = smith_waterman(a, b)
        c_score, c_cell, c_rows = run_c(p)
        ok = (c_score == score and c_cell == cell)
        if c_rows and len(a) <= 40 and len(b) <= 40:
            ok = ok and c_rows == H
        print(f"  {'ok  ' if ok else 'FAIL'} {n:16s} C: {c_score} at {c_cell}   Python: {score} at {cell}"
              + ("   (matrix compared)" if c_rows else ""))
        bad += (not ok)
    print("\nAll C results match the reference." if not bad else f"\n{bad} mismatch(es).")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
