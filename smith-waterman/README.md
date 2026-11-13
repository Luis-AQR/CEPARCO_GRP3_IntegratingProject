# Smith-Waterman: Full Alignment Reconstruction, C vs CUDA

CEPARCO Integrating Project, Group 3 (Alcasid, Bilanes, Corpuz, Roa)

This project implements the Smith-Waterman local sequence alignment algorithm twice,
a sequential **C** baseline and a **CUDA** version that fills the scoring matrix one
anti-diagonal at a time, and compares them on identical inputs
(`Speedup = T_sequentialC / T_CUDA`). See the project proposal for the algorithm and the
parallelisation approach.

**Status: Week 8 - sequential C, matrix filling.** `c/sw_seq.c` reads the FASTA input,
fills the scoring matrix row by row and reports the maximum score and its cell, and matches
the independent Python reference on every input. Traceback and baseline timings come in
Week 9, CUDA from Week 10.

## Contents

```
data/                      test_*.fasta  hand-checkable inputs
                           bench_*.fasta 1,000 / 5,000 / 10,000-base benchmark inputs
scripts/gen_sequences.py   reproducible generator for the benchmark inputs
scripts/check_expected.py  independent Python check of the expected results below
scripts/check_c_week8.py   compares the C program with the Python reference (score, cell, matrix)
c/sw_seq.c                 sequential C baseline (Week 8: input, matrix filling, score tracking)
bin/                       compiled sw_seq.exe (created by the build command below; do not commit)
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
The C program also accepts `--match N`, `--mismatch N` and `--gap N` (positive integers) to change
them at run time, since the proposal lists the scoring parameters as inputs; the defaults above
are used for every documented result and every check.

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
| `bench_1000.fasta`, `bench_5000.fasta`, `bench_10000.fasta` | random DNA with a planted, mutated shared region (20% of the length, 10% mutations); sizes from the proposal | scored by the C baseline in Week 8, see *Week 8 results* |

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
repository, validates the inputs, and (Week 8) builds and checks the C program.
Environment recorded on 12 Nov 2026: Tesla T4, driver 580.82.07, CUDA 13.0, Python 3.13.

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

## Build and run the C baseline (Week 8)

`sw_seq.c` is plain C99 (no POSIX calls), so it builds with MSVC or MinGW-w64 on Windows.
From `smith-waterman/`, with one compiler on `PATH` (for MSVC use the *x64 Native Tools Command Prompt for VS*):

```bat
mkdir bin
gcc -std=c99 -O2 -Wall -Wextra -o bin\sw_seq.exe c\sw_seq.c
```

or with MSVC:

```bat
mkdir bin
cl /nologo /O2 /W4 /Fobin\ /Febin\sw_seq.exe c\sw_seq.c
```

Run and check (use `python3` instead of `python` if that is what your machine has):

```bat
bin\sw_seq.exe data\test_paper.fasta                 :: prints the matrix (inputs up to 40 x 40), max score and checksum
bin\sw_seq.exe data\bench_1000.fasta --no-matrix     :: larger inputs: score, cell and checksum only
python scripts\check_expected.py                     :: Week 7 check of the expected results
python scripts\check_c_week8.py                      :: C vs Python reference (add --all for bench_5000)
```

`check_c_week8.py` looks for `bin\sw_seq.exe`; use `--exe PATH` if the program is somewhere else.
It should end with `All C results match the reference.` Options: `--no-matrix`, `--match N`,
`--mismatch N`, `--gap N`. Use the optimised build (`-O2` / `/O2`) for anything timed, and take
final timings on Colab only.

Expected for `test_paper.fasta`: max score 3 at H(4, 5), and the matrix shown above.

### What the program does

* reads a FASTA file with exactly two records and rejects non-letters, missing or extra records and empty sequences;
* fills the `(M+1) x (N+1)` matrix row by row with `H(i,j) = max(0, diag, up, left)` and tracks the maximum score and its cell (the *first* maximum in row-major order; "no local alignment" if every cell is 0);
* prints the matrix for inputs up to 40 x 40, then `Max score`, `Max cell` and a **checksum** (the sum of all cells), which lets the Week 10 CUDA matrix be compared with this one on the large inputs without printing them;
* rejects sequences longer than 46,000 (about 8.5 GB of matrix) and reports an allocation failure instead of crashing.

### Week 8 results

Scoring +1 / -1 / -2. Every result equals the independent Python reference (`scripts/check_expected.py`).

| File | M x N | C: max score at cell | Checksum | Python agrees |
|------|-------|----------------------|----------|---------------|
| `test_paper` | 4 x 5 | **3** at H(4,5) | 8 | yes (whole matrix) |
| `test_identical` | 8 x 8 | **8** at H(8,8) | 100 | yes (whole matrix) |
| `test_nomatch` | 8 x 8 | **0**, no cell | 0 | yes (whole matrix) |
| `test_gap` | 12 x 13 | **10** at H(12,13) | 222 | yes (whole matrix) |
| `bench_1000` | 1,000 x 1,000 | **160** at H(784,995) | 1,563,453 | yes |
| `bench_5000` | 5,000 x 5,000 | **805** at H(2313,1448) | 1,445,075,095 | yes (`--all`) |
| `bench_10000` | 10,000 x 10,000 | **1536** at H(5403,9807) | 1,009,197,635 | not run (pure Python is too slow) |

These are the values the CUDA version must reproduce in Weeks 10-11. Not done yet: traceback
(no aligned strings or coordinates) and timing code, both planned for Week 9.

## Plan

| Week | Milestone | Deliverable |
|------|-----------|-------------|
| 7 | Finalise parameters and test sequences; Colab + GitHub setup | repository with README and sample inputs (done) |
| 8 | Sequential C: input reading, matrix filling, score tracking | **this repository**: C program that prints the score matrix and max score |
| 9 | C traceback; baseline execution times | verified C program, baseline timings |
| 10 | Mid-project demo | C baseline + working CUDA matrix filling |
| 11 | CUDA traceback and alignment reconstruction; compare with C | full CUDA alignment with correctness comparison |
| 12 | Larger inputs, `nsys` profiling, speedup; tiled wavefront if launch overhead is high | timing tables, speedup, profiling notes |
| 13 | Documentation and demo | final submission |

## Reference

[1] T. F. Smith and M. S. Waterman, "Identification of common molecular subsequences,"
J. Mol. Biol., vol. 147, no. 1, pp. 195-197, 1981.