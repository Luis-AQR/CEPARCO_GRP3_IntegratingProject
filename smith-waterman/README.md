# Smith-Waterman: Full Alignment Reconstruction, C vs CUDA

CEPARCO Integrating Project, Group 3 (Alcasid, Bilanes, Corpuz, Roa)

This project implements the Smith-Waterman local sequence alignment algorithm twice,
a sequential **C** baseline and a **CUDA** version that fills the scoring matrix one
anti-diagonal at a time, and compares them on identical inputs
(`Speedup = T_sequentialC / T_CUDA`). See the project proposal for the algorithm and the
parallelisation approach.

**Status: Week 7 - project setup.** The scoring parameters and test sequences are
finalised and the Colab notebook and repository are set up. No algorithm code yet
(C in Weeks 8-9, CUDA from Week 10).

## Contents

```
data/                      test_*.fasta  hand-checkable inputs
                           bench_*.fasta 1,000 / 5,000 / 10,000-base benchmark inputs
scripts/gen_sequences.py   reproducible generator for the benchmark inputs
scripts/check_expected.py  independent Python check of the expected results below
colab/SW_Week7_Setup.ipynb Colab notebook: GPU / nvcc check, clone, validate inputs
results/                   benchmark output goes here (Weeks 9 and 12)
```

## Scoring parameters (final)

Linear gap penalty, exactly Eq. 1 and Eq. 2 of the proposal:

| Parameter | Value |
|-----------|-------|
| match, *m* | +1 |
| mismatch, *mu* | -1 |
| gap, *d* | -2 (once per gap column) |

These are the values of the proposal's worked example (Section III), so the demo can
reproduce that example exactly. Run time does not depend on the parameter values, so they
are not a performance variable. Both implementations must use the same values on every run.

Conventions both implementations will follow (so their outputs can be compared):

* **Matrix:** `(M+1) x (N+1)`, sequence *a* (length M) down the rows, sequence *b*
  (length N) across the columns; row 0 and column 0 are 0.
* **Best cell:** the *first* cell in row-major order holding the maximum score.
* **Traceback:** start at the best cell, stop at a 0; at each cell try diagonal, then up
  (gap in *b*), then left (gap in *a*) - the order of the proposal's pseudocode.
* **Several optimal alignments:** the strings may legitimately differ, so a result is
  judged by score, coordinates, and by re-scoring the reconstructed alignment (Section IV).
* **Integers:** 32-bit scores; the full matrix is kept for traceback, so a 10,000 x 10,000
  input needs about 380 MB.

## Input format

A FASTA file with **exactly two** records. Sequence 1 runs down the rows, sequence 2 across
the columns. Letters only (case-insensitive); multi-line records and CRLF line endings are fine.

## Test sequences

Expected values below were derived by hand from the recurrence and are re-checked by
`python3 scripts/check_expected.py`.

| File | Purpose | Expected (match +1, mismatch -1, gap -2) |
|------|---------|------------------------------------------|
| `test_paper.fasta` | proposal's worked example: `AGCT` (rows) vs `ATGCT` (columns) | max score **3** at H(4,5); alignment `GCT` / `GCT`; a[2..4], b[3..5] |
| `test_identical.fasta` | identical 8-mers | score **8** at H(8,8); whole sequences align |
| `test_nomatch.fasta` | `AAAAAAAA` vs `CCCCCCCC` | score **0**; no local alignment |
| `test_gap.fasta` | `ACGTCATGCATG` vs the same with one extra `C` inserted | score **10** at H(12,13); `ACGTCA-TGCATG` / `ACGTCACTGCATG` (12 matches, 1 gap) |
| `bench_1000.fasta`, `bench_5000.fasta`, `bench_10000.fasta` | random DNA with a planted, mutated shared region (20% of the length, 10% mutations); sizes from the proposal | scored by the C baseline in Week 9 |

Expected score matrix for `test_paper.fasta` (rows A, G, C, T; columns A, T, G, C, T):

```
        0  1  2  3  4  5
           A  T  G  C  T
  0     0  0  0  0  0  0
  1 A   0  1  0  0  0  0
  2 G   0  0  0  1  0  0
  3 C   0  0  0  0  2  0
  4 T   0  0  1  0  0  3
```

**Correction to the proposal's Figure 7:** the figure shows 0 at H(4,2), but it is **1**
(a[4] = T matches b[2] = T and H(3,1) = 0). The maximum score and the `GCT` alignment are
unchanged. Figures 7 and 8 should be updated.

The benchmark inputs are reproducible:

```bash
python3 scripts/gen_sequences.py --length 10000 --seed 1 --out data/bench_10000.fasta
```

## Setup

**Google Colab** (all final timings come from here, in one session, so C and CUDA are compared
on the same machine): open `colab/SW_Week7_Setup.ipynb`, choose Runtime -> Change runtime type
-> GPU (usually a T4; the exact GPU depends on availability, so record which one you got),
and run the cells. It checks the GPU and `nvcc`, compiles a tiny CUDA test, clones this
repository, and validates the inputs.

**GitHub** (once, by one member):

```bash
git init
git add .
git commit -m "Week 7: project setup, scoring parameters, test sequences, Colab notebook"
git branch -M main
git remote add origin https://github.com/<user>/<repo>.git
git push -u origin main
```

Then add the other members under Settings -> Collaborators, and put the repository URL
into the `REPO_URL` cell of the notebook. Members with an NVIDIA GPU may also build locally
(Visual Studio 2019/2022 + CUDA Toolkit), but final timings come from Colab only.

## Plan

| Week | Milestone | Deliverable |
|------|-----------|-------------|
| 7 | Finalise parameters and test sequences; Colab + GitHub setup | **this repository** |
| 8 | Sequential C: input reading, matrix filling, score tracking | C program that prints the score matrix and max score |
| 9 | C traceback; baseline execution times | verified C program, baseline timings |
| 10 | Mid-project demo | C baseline + working CUDA matrix filling |
| 11 | CUDA traceback and alignment reconstruction; compare with C | full CUDA alignment with correctness comparison |
| 12 | Larger inputs, `nsys` profiling, speedup; tiled wavefront if launch overhead is high | timing tables, speedup, profiling notes |
| 13 | Documentation and demo | final submission |

## Reference

[1] T. F. Smith and M. S. Waterman, "Identification of common molecular subsequences,"
J. Mol. Biol., vol. 147, no. 1, pp. 195-197, 1981.
