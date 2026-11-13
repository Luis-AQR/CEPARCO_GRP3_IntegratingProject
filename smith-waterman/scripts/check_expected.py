"""
check_expected.py - Week 7: confirm the expected results of the test sequences.

A tiny, independent pure-Python Smith-Waterman written straight from the
recurrence in the proposal (Eq. 1 / Eq. 2).  Its job is to prove, before any C or CUDA
exists, that the test sequences in data/ really give the results documented in
the README.  Later the C and CUDA outputs are compared against the same values.

Scoring (finalised): match +1, mismatch -1, gap -2 (linear).
Traceback priority: diagonal > up > left (as in the proposal's pseudocode).
Maximum cell: first one in row-major order (strict '>').

Usage:  python scripts\\check_expected.py          (run from the repo root)
"""
import os
import sys

M_SCORE, MU, D = 1, 1, 2


def read_fasta_pair(path):
    recs, cur = [], None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                cur = []
                recs.append(cur)
            elif line:
                if cur is None:
                    raise ValueError("sequence data before first '>' header")
                cur.append(line.upper())
    seqs = ["".join(r) for r in recs]
    if len(seqs) != 2 or not all(seqs):
        raise ValueError(f"expected exactly two non-empty sequences, found {len(seqs)}")
    for s in seqs:
        if not s.isalpha():
            raise ValueError("sequence contains non-letter characters")
    return seqs[0], seqs[1]


def smith_waterman(a, b):
    M, N = len(a), len(b)
    H = [[0] * (N + 1) for _ in range(M + 1)]
    best, bi, bj = 0, 0, 0
    for i in range(1, M + 1):
        for j in range(1, N + 1):
            s = M_SCORE if a[i - 1] == b[j - 1] else -MU
            h = max(0, H[i - 1][j - 1] + s, H[i - 1][j] - D, H[i][j - 1] - D)
            H[i][j] = h
            if h > best:
                best, bi, bj = h, i, j
    ra, rb = [], []
    i, j = bi, bj
    while best > 0 and i > 0 and j > 0 and H[i][j] > 0:
        s = M_SCORE if a[i - 1] == b[j - 1] else -MU
        if H[i][j] == H[i - 1][j - 1] + s:
            ra.append(a[i - 1]); rb.append(b[j - 1]); i -= 1; j -= 1
        elif H[i][j] == H[i - 1][j] - D:
            ra.append(a[i - 1]); rb.append("-"); i -= 1
        else:
            ra.append("-"); rb.append(b[j - 1]); j -= 1
    al_a, al_b = "".join(reversed(ra)), "".join(reversed(rb))
    start = (i + 1, j + 1) if al_a else None
    return H, best, (bi, bj), start, al_a, al_b


# Values derived by hand from the recurrence (see README, "Test sequences").
EXPECTED = {
    "test_paper.fasta": dict(
        score=3, end=(4, 5), start=(2, 3), a="GCT", b="GCT",
        matrix=[[0, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 2, 0],
                [0, 0, 1, 0, 0, 3]]),
    "test_identical.fasta": dict(score=8, end=(8, 8), start=(1, 1), a="ACGTACGT", b="ACGTACGT"),
    "test_nomatch.fasta":   dict(score=0, end=(0, 0), start=None, a="", b=""),
    "test_gap.fasta":       dict(score=10, end=(12, 13), start=(1, 1),
                                 a="ACGTCA-TGCATG", b="ACGTCACTGCATG"),
}

BENCH_LENGTHS = {"bench_1000.fasta": 1000, "bench_5000.fasta": 5000, "bench_10000.fasta": 10000}


def main():
    data = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
    failures = 0

    print("Test sequences (pure-Python reference, match +1 / mismatch -1 / gap -2):")
    for name, exp in EXPECTED.items():
        a, b = read_fasta_pair(os.path.join(data, name))
        H, score, end, start, al_a, al_b = smith_waterman(a, b)
        got = dict(score=score, end=end, start=start, a=al_a, b=al_b)
        problems = [k for k in ("score", "end", "start", "a", "b") if got[k] != exp[k]]
        if "matrix" in exp and H != exp["matrix"]:
            problems.append("matrix")
        status = "ok  " if not problems else "FAIL"
        failures += bool(problems)
        print(f"  {status} {name:22s} {len(a):>3} x {len(b):<3} score={score:<3} end={end} "
              f"{al_a or '(no alignment)'}" + (f"  mismatch in: {problems}" if problems else ""))

    print("\nBenchmark sequences (format check only; scored by the C baseline in Week 9):")
    for name, n in BENCH_LENGTHS.items():
        try:
            a, b = read_fasta_pair(os.path.join(data, name))
            ok = len(a) == n and len(b) == n and set(a + b) <= set("ACGT")
            detail = f"{len(a)} x {len(b)}, alphabet {''.join(sorted(set(a + b)))}"
        except ValueError as e:
            ok, detail = False, str(e)
        failures += (not ok)
        print(f"  {'ok  ' if ok else 'FAIL'} {name:22s} {detail}")

    print("\nAll checks passed." if not failures else f"\n{failures} check(s) FAILED.")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
