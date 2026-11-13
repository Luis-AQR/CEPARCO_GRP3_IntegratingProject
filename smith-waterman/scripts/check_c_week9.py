"""
check_c_week9.py - Week 9: compare the traceback of the sequential C program with the
independent pure-Python reference (scripts/check_expected.py).

For every input it checks:
  * the C program exits with 0 and prints "Re-score : <score> (OK)";
  * score and best cell equal the Python reference;
  * the aligned strings and the coordinates equal the Python traceback (both use the
    same priority: diagonal, then up, then left);
  * the aligned strings re-score (in Python, from the strings only) to the max score
    and the un-gapped strings equal the slices of the input sequences.
test_nomatch must report "No local alignment" and print no alignment.

Usage (from the repo root, after building bin\\sw_seq.exe):
    python scripts\\check_c_week9.py             tests + bench_1000
    python scripts\\check_c_week9.py --all       also bench_5000 (slower: pure Python)
    python scripts\\check_c_week9.py --exe PATH  use a specific executable
"""
import argparse
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_expected  # noqa: E402
from check_expected import read_fasta_pair, smith_waterman  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DATA = os.path.join(ROOT, "data")


def find_exe(explicit=None):
    candidates = [explicit] if explicit else [
        os.path.join(ROOT, "bin", "sw_seq.exe"),
        os.path.join(ROOT, "bin", "sw_seq"),
    ]
    for c in candidates:
        if c and os.path.isfile(c):
            return os.path.abspath(c)
    sys.exit("error: compiled program not found (looked for bin\\sw_seq.exe). "
             "Build it first, or pass --exe PATH.")


def rescore(al_a, al_b):
    """Score an alignment from the two aligned strings only."""
    if len(al_a) != len(al_b):
        return None
    total = 0
    for x, y in zip(al_a, al_b):
        if x == "-" and y == "-":
            return None
        if x == "-" or y == "-":
            total -= check_expected.D
        else:
            total += check_expected.M_SCORE if x == y else -check_expected.MU
    return total


def run_c(exe, path):
    """Run sw_seq and parse score, cell, alignment strings, range and re-score line."""
    proc = subprocess.run([exe, path, "--no-matrix"], capture_output=True, text=True)
    out = proc.stdout
    info = {"exit": proc.returncode}
    info["score"] = int(re.search(r"Max score:\s*(-?\d+)", out).group(1))
    m = re.search(r"Max cell\s*:\s*H\((\d+), (\d+)\)", out)
    info["cell"] = (int(m.group(1)), int(m.group(2))) if m else (0, 0)
    info["no_alignment"] = "No local alignment" in out
    m = re.search(r"Range\s*:\s*a\[(\d+)\.\.(\d+)\], b\[(\d+)\.\.(\d+)\]", out)
    info["range"] = tuple(int(g) for g in m.groups()) if m else None
    a_parts, b_parts = [], []
    for line in out.splitlines():
        m = re.match(r"^\s+([ab])\s+\d+\s+(\S+)\s*$", line)
        if m:
            (a_parts if m.group(1) == "a" else b_parts).append(m.group(2))
    info["a_aln"] = "".join(a_parts)
    info["b_aln"] = "".join(b_parts)
    m = re.search(r"Re-score\s*:\s*(-?\d+) \((OK|MISMATCH)", out)
    info["rescore_line"] = (int(m.group(1)), m.group(2)) if m else None
    return info


def check_one(exe, name):
    path = os.path.join(DATA, name + ".fasta")
    a, b = read_fasta_pair(path)
    _, score, cell, start, ref_a, ref_b = smith_waterman(a, b)
    c = run_c(exe, path)
    problems = []

    if c["exit"] != 0:
        problems.append(f"exit code {c['exit']}")
    if c["score"] != score:
        problems.append(f"score {c['score']} != {score}")
    if c["cell"] != cell:
        problems.append(f"cell {c['cell']} != {cell}")

    if score == 0:
        if not c["no_alignment"] or c["a_aln"] or c["b_aln"]:
            problems.append("expected 'No local alignment' and no alignment")
        return c, score, cell, problems

    # alignment strings and coordinates vs the Python traceback
    ref_range = (start[0], cell[0], start[1], cell[1])
    if c["a_aln"] != ref_a or c["b_aln"] != ref_b:
        problems.append("aligned strings differ from the Python traceback")
    if c["range"] != ref_range:
        problems.append(f"range {c['range']} != {ref_range}")

    # validity: independent re-score and slices of the input sequences
    rs = rescore(c["a_aln"], c["b_aln"])
    if rs != score:
        problems.append(f"Python re-score {rs} != {score}")
    if c["range"]:
        i1, i2, j1, j2 = c["range"]
        if c["a_aln"].replace("-", "") != a[i1 - 1:i2]:
            problems.append("a alignment is not a[i1..i2]")
        if c["b_aln"].replace("-", "") != b[j1 - 1:j2]:
            problems.append("b alignment is not b[j1..j2]")
    if c["rescore_line"] != (score, "OK"):
        problems.append(f"C re-score line {c['rescore_line']}")
    return c, score, cell, problems


def main():
    ap = argparse.ArgumentParser(description="Check the C traceback against the Python reference.")
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
        c, score, cell, problems = check_one(exe, n)
        detail = "no alignment" if c["no_alignment"] else f"len {len(c['a_aln'])} range {c['range']}"
        print(f"  {'ok  ' if not problems else 'FAIL'} {n:16s} score {score} at {cell}  {detail}"
              + (f"   problems: {'; '.join(problems)}" if problems else ""))
        bad += bool(problems)
    print("\nAll C results match the reference." if not bad else f"\n{bad} mismatch(es).")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
