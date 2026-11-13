# Smith-Waterman: Full Alignment Reconstruction, C vs CUDA
## Weekly Progress Report (through Week 9, traceback)

**Course:** Advanced Computer Architecture 1 (CEPARCO), Integrating Project
**Group 3:** Alcasid, Aizy; Bilanes, Kendrick; Corpuz, Gerald; Roa, Luis

The project implements the Smith-Waterman local alignment algorithm twice, a sequential **C**
baseline and a **CUDA** version that fills the scoring matrix one anti-diagonal at a time, and
compares them on identical inputs (`Speedup = T_sequentialC / T_CUDA`). Algorithm and
parallelisation approach: see the project proposal.

---

## Background: the Smith-Waterman algorithm

### What it is

Smith-Waterman, published by Temple Smith and Michael Waterman in 1981 [1], is a
dynamic-programming algorithm for **local sequence alignment**: given two sequences, it
finds the pair of subsequences, one from each, that are most similar to each other, even
when the sequences as a whole are not similar. The matching regions do not have to be
identical; mismatches and gaps (insertions and deletions) are allowed but penalised. It was
first developed for finding common molecular subsequences in DNA, RNA and protein [1], and it
is guaranteed to return the best local alignment under the chosen scoring scheme [4].

A *global* alignment, by contrast, forces the whole of both sequences to line up. That is the
wrong question when only a part of each sequence is related, which is why a local method is
used here.

### How it works

The algorithm builds a scoring matrix `H` of `(M+1) x (N+1)` cells, where `H(i,j)` is the best
score of any alignment that ends at the i-th symbol of sequence *a* and the j-th symbol of
sequence *b*. It has three steps:

1. **Initialisation.** Row 0 and column 0 are set to 0.
2. **Matrix filling.** Each remaining cell is the maximum of four candidates:

   ```
   H(i,j) = max( 0,
                 H(i-1,j-1) + s(a_i, b_j),     diagonal: match (+m) or mismatch (-mu)
                 H(i-1,j)   - d,               up:       gap in b (deletion)
                 H(i,j-1)   - d )              left:     gap in a (insertion)
   ```

   The `0` is the "reset". A badly scoring prefix is dropped and a fresh alignment can start
   at any cell, so scores never go negative. This one term is what makes the alignment local
   rather than global.
3. **Traceback.** Start at the highest-scoring cell and follow the choices that produced each
   cell back to a 0. The path read in order is the aligned pair of subsequences, with `-` where
   a gap was inserted.

Example (this project's test case `test_paper`, with +1 / -1 / -2): `AGCT` against `ATGCT`
gives a maximum score of 3 at H(4,5), and the traceback returns `GCT` / `GCT`. The leading `A`
of `AGCT` and the leading `AT` of `ATGCT` are not part of the best local alignment and are left out.
The complete matrix for this example is shown under *Test sequences and expected results*.

### Why it needs parallelising

Every cell of the matrix is computed, so both running time and memory grow with `M x N`
(the full matrix is kept for the traceback). A 10,000 x 10,000 comparison is 100 million cells.
Heuristic search tools (BLAST is the best known) skip most of the matrix to run faster, but
they can miss the optimal alignment. To keep maximum sensitivity, this project processes the
full matrix and speeds it up with parallelism instead [2][3][5].

Parallelism is possible because `H(i,j)` depends only on its left, upper and upper-left
neighbours. All cells on the same anti-diagonal (`i + j` constant) therefore do not depend on
each other and can be computed at the same time by GPU threads. The CUDA version in this
project does exactly that, one anti-diagonal per kernel launch.

### Use cases

Any problem that asks "which parts of these two sequences are alike, even if the rest is not?"
fits, whenever the sequences can contain small differences:

| Field | What is compared | What a local alignment finds |
|-------|------------------|------------------------------|
| **Bioinformatics / genetics** | DNA, RNA or protein sequences; a query against a large sequence database [2][5] | shared genes, conserved regions or motifs between otherwise different sequences |
| **Plagiarism detection** | words in sentences or documents; notes in songs | copied or lightly reworded passages and melodies, even with words or notes inserted, removed or changed |
| **Cybersecurity** | malware binaries | shared code fragments and attack patterns between samples, even when the surrounding code differs |

Exactness matters in these settings: a missed match can be a missed gene, an undetected copied
passage or an unrecognised malware variant.

---

## Progress at a glance

| Week | Milestone (from the proposal) | Status |
|------|-------------------------------|--------|
| 7 | Finalise scoring parameters and test sequences; set up Colab notebook and GitHub repository | **Done** - Colab confirmed; merged into `main` |
| **8** | Sequential C: input reading, matrix filling, score tracking | **Done** (see *Week 8 report*) - C output matches the Python reference on all 6 inputs checked |
| 9 | C traceback; record baseline execution times | **Done** - C traceback verified against the Python reference; baseline times recorded on Colab (T4 session, 13 Nov 2026), and the Colab correctness and sanitizer checks passed (see *Week 9 report*) |
| 10 | Mid-project demo: C baseline + working CUDA matrix filling | Not started |
| 11 | CUDA traceback and alignment reconstruction; compare with C | Not started |
| 12 | Larger inputs, `nsys` profiling, speedup; tiled wavefront if launch overhead is high | Not started |
| 13 | Documentation and demo | Not started |

---

## Week 9 report (completed)

**Milestone:** C traceback; record baseline execution times.
**Deliverable:** verified C program with traceback, and baseline execution times (`results/week9_baseline.csv`, `results/week9_baseline.md`).

### What was done

| # | Task | Result |
|---|------|--------|
| 1 | Traceback | `traceback()` starts at the best cell and follows diagonal, then up, then left, until a 0. Moves are recomputed from `H` and the sequences, so no direction matrix is stored. The aligned strings are built backwards in buffers of `M + N + 1`. |
| 2 | Alignment output | Prints the length (matches / mismatches / gaps), the 1-based ranges `a[i1..i2]`, `b[j1..j2]`, and the aligned pair in 60-column blocks with match bars. |
| 3 | Independent re-score | `rescore_alignment()` scores the aligned strings only (it does not read `H`), and checks that the un-gapped strings equal `a[i1..i2]` and `b[j1..j2]`. The program prints `Re-score : <n> (OK)` and exits with 1 if it differs from the max score. |
| 4 | Flags | `--no-traceback` gives exactly the Week 8 output (checked against the Week 8 build on `test_paper`, `test_nomatch`, `bench_1000`). |
| 5 | Verification script | `scripts/check_c_week9.py` runs `sw_seq` and compares score, cell, aligned strings and ranges with the Python traceback in `check_expected.py` (same priority), and re-scores the parsed strings in Python. |
| 6 | Timing code | `sw_seq` times allocation, fill and traceback separately with a monotonic clock; opt-in flags `--time`, `--repeat N`, `--warmup`, `--csv` (the default output is unchanged). |
| 7 | Baseline run | `scripts/run_baseline.py` and `colab/SW_Week9_Timing.ipynb` take the official timings in one Colab session and write the results files. |

### Results

Scoring +1 / -1 / -2. Every alignment re-scores to the max score.

| File | Max score at cell | Range | Alignment | Re-score |
|------|-------------------|-------|-----------|----------|
| `test_paper` | 3 at H(4,5) | a[2..4], b[3..5] | `GCT` / `GCT` | 3 OK |
| `test_identical` | 8 at H(8,8) | a[1..8], b[1..8] | whole 8-mers | 8 OK |
| `test_nomatch` | 0, no cell | n/a | "No local alignment", traceback skipped | n/a |
| `test_gap` | 10 at H(12,13) | a[1..12], b[1..13] | `ACGTCA-TGCATG` / `ACGTCACTGCATG` | 10 OK |
| `bench_1000` | 160 at H(784,995) | a[566..784], b[776..995] | length 225 (198 matches, 16 mismatches, 11 gaps) | 160 OK |
| `bench_5000` | 805 at H(2313,1448) | a[1298..2313], b[439..1448] | length 1043 (954 matches, 29 mismatches, 60 gaps) | 805 OK |
| `bench_10000` | 1536 at H(5403,9807) | a[3387..5403], b[7808..9807] | length 2069 (1863 matches, 85 mismatches, 121 gaps) | 1536 OK (C only) |

`check_c_week9.py --all` ends with `All C results match the reference.` (strings and ranges equal the
Python traceback on all five inputs it covers); `check_c_week8.py --all` still passes, so the fill and
checksums are unchanged. Custom scoring also re-scores correctly (`test_gap --match 2 --mismatch 1
--gap 3`: 21; `bench_1000 --match 3 --mismatch 2 --gap 1`: 1337).

### Technical decisions

* **Priority** diagonal, then up (gap in b), then left (gap in a); the CUDA traceback (Week 11) must use the same order.
* **Validity** is judged by the re-score and the coordinates. String equality with the Python reference is also checked because both sides use the same fixed priority.
* **Checks can fail.** A scratch build whose traceback prefers up over diagonal on a tie was caught by `check_c_week9.py` on `bench_1000` (strings differ, re-score still OK); a wrong gap penalty in the re-score was caught as `MISMATCH`. Swapping only up and left changed nothing on these inputs, because an up/left tie hardly ever occurs once the diagonal is tried first.

### Baseline timing (official, Colab)

C is the timing baseline for `Speedup = T_sequentialC / T_CUDA`; the Python reference decides correctness for both C and CUDA. All official timings come from **one Colab session**, so the whole team uses the same environment and the CUDA timings can be taken on the same kind of machine. Local Windows runs are for development only and are not reported.

**Environment** (recorded in `results/week9_baseline.md`): 13 Nov 2026 12:14; Linux 6.6.159, glibc 2.39; Intel Xeon @ 2.00 GHz, 2 logical CPUs (the C code is single-threaded); Tesla T4, 15360 MiB, driver 580.82.07; gcc 13.3.0, `-std=c99 -O2 -Wall -Wextra`; Python 3.13.16.

**Method**

* **What is timed** (per run): `alloc` = `malloc` + zeroing the whole matrix (so the page faults are charged here, not to the fill; the zeroing goes through a volatile pointer because GCC otherwise turns `malloc` + `memset` into a lazy `calloc`), `fill` = `fill_matrix()` only, `traceback` = `traceback()` only, `total` = the three added. Reading the FASTA file and printing are not timed.
* **Protocol:** warm-up run (run 0, not in the statistics) + 5 timed runs for 1,000 and 5,000, warm-up + 3 for 10,000, each with and without the traceback. The script stops if any run does not reproduce the verified score, cell and checksum; none did.
* **Timer:** `clock_gettime(CLOCK_MONOTONIC)` on Linux/Colab, `QueryPerformanceCounter` on Windows.

**Results** (seconds, median of the timed runs; full = fill + traceback)

| Input | Cells | Runs | Fill median (min - max) | Alloc + zero | Traceback | Total |
|-------|-------|------|-------------------------|--------------|-----------|-------|
| `bench_1000` | 1,002,001 | 5 | **0.0020** (0.0020 - 0.0021) | 0.0001 | 0.000003 | **0.0021** |
| `bench_5000` | 25,010,001 | 5 | **0.0495** (0.0485 - 0.0512) | 0.0464 | 0.000049 | **0.0960** |
| `bench_10000` | 100,020,001 | 3 | **0.2111** (0.2072 - 0.2159) | 0.1925 | 0.000092 | **0.4086** |

Fill-only runs (`--no-traceback`): fill 0.0019 / 0.0508 / 0.1980 s and total 0.0021 / 0.0998 / 0.3859 s for 1,000 / 5,000 / 10,000, within a few percent of the full runs, which confirms the traceback is negligible.

**What the numbers show**

* The fill costs about 2 ns per cell at every size (2.0, 1.98 and 2.11 ns), so the time grows linearly with `M x N`, as expected for a full-matrix method.
* Allocating and zeroing the matrix is not small: at 5,000 and 10,000 it is about as long as the fill (95 MB and 381 MB of page faults). The CUDA timings must therefore report `cudaMalloc`, any memset, copies and kernel time separately, and the speedup must say what it compares: **fill vs kernel** (the part that is parallelised) and **total vs total** (end to end). The denominators are the Fill and Total columns above.
* The spread between runs is within about 5%, so the medians are stable enough to quote to two significant figures.
* These numbers are for this Colab machine only. Do not mix them with the laptop numbers or with a later session; if the baseline is re-run, re-run the CUDA timings in the same session.

### Open items

- [x] Colab run of `SW_Week9_Timing.ipynb` (same session as the timings): `check_expected.py` printed `All checks passed.`; `check_c_week8.py --all` and `check_c_week9.py --all` each ended with `All C results match the reference.` (six inputs, `bench_5000` included, same scores, cells and checksums as the local runs); the AddressSanitizer + UndefinedBehaviorSanitizer build ran `test_paper`, `test_identical`, `test_nomatch`, `test_gap` and `bench_1000` with `--repeat 2 --warmup`, every alignment re-scored OK, and the cell ended with `sanitizer run finished: no reports`.
- [x] Merge the Week 7, Week 8 and Week 9 branches into `main`: the Week 7 and Week 8 branches were already in `main`; the Week 9 branch is merged with this commit.

### Reproducing the Week 9 checks

From `smith-waterman/`, after building `bin\sw_seq.exe` as in the Week 8 section:

```bat
bin\sw_seq.exe data\test_gap.fasta --no-matrix       :: alignment, ranges, Re-score line
python scripts\check_c_week9.py                       :: C traceback vs Python (add --all for bench_5000)
bin\sw_seq.exe data\bench_1000.fasta --no-matrix --repeat 5 --warmup --csv   :: timing (development only)
```

On Colab (official timings), open `colab/SW_Week9_Timing.ipynb`, or from `smith-waterman/`:

```bash
gcc -std=c99 -O2 -Wall -Wextra -o bin/sw_seq c/sw_seq.c
python3 scripts/run_baseline.py --exe bin/sw_seq      # writes results/week9_baseline.csv / .md
```

---

## Week 8 report

**Milestone:** sequential C: input reading, matrix filling, score tracking.
**Deliverable:** a working C program that outputs the score matrix and the maximum score (`c/sw_seq.c`).

### What was done

| # | Task | Result |
|---|------|--------|
| 1 | FASTA reader | `read_fasta_pair()` accepts exactly two records, upper-cases letters, handles multi-line records, blank lines and CRLF, and rejects non-letters, missing/extra records and empty sequences with an error message. |
| 2 | Matrix filling | `fill_matrix()` computes `H(i,j) = max(0, diag, up, left)` row by row with the Week 7 scoring (+1 / -1 / -2). The matrix is `calloc`-ed, so row 0 and column 0 are the required 0s. |
| 3 | Score tracking | Maximum score and its cell are updated with a strict `>`, so the reported cell is the **first** maximum in row-major order (the shared C/CUDA convention). If every cell is 0 the program reports "no local alignment". |
| 4 | Output | Prints file, lengths, scoring, matrix size in MB, the full matrix for inputs up to 40 x 40, then `Max score` and `Max cell`. `--no-matrix` suppresses the matrix. |
| 5 | Verification against the Python reference | `scripts/check_c_week8.py` runs the compiled program and compares score and best cell with `check_expected.py` on `test_paper`, `test_identical`, `test_nomatch`, `test_gap`, `bench_1000` and (with `--all`) `bench_5000`; the full matrix is compared for the four small tests. |

### Results

Scoring +1 / -1 / -2. Every row below was produced by `sw_seq`; "Python" is the independent reference.

| File | Size (M x N) | Matrix | C: max score at cell | Python: score at cell | Match |
|------|--------------|--------|----------------------|-----------------------|-------|
| `test_paper` | 4 x 5 | 30 cells | **3** at H(4,5) | 3 at (4,5) | yes (matrix compared) |
| `test_identical` | 8 x 8 | 81 cells | **8** at H(8,8) | 8 at (8,8) | yes (matrix compared) |
| `test_nomatch` | 8 x 8 | 81 cells | **0**, no cell | 0 at (0,0) | yes (matrix compared) |
| `test_gap` | 12 x 13 | 182 cells | **10** at H(12,13) | 10 at (12,13) | yes (matrix compared) |
| `bench_1000` | 1,000 x 1,000 | 3.8 MB | **160** at H(784,995) | 160 at (784,995) | yes |
| `bench_5000` | 5,000 x 5,000 | 95.4 MB | **805** at H(2313,1448) | 805 at (2313,1448) | yes |
| `bench_10000` | 10,000 x 10,000 | 381.5 MB | **1536** at H(5403,9807) | not run (pure Python is too slow) | C only |

The `test_paper` matrix printed by the C program is identical to the expected matrix in the
*Test sequences* section below, including H(4,2) = 1.

These scores and cells are the values the CUDA version must reproduce in Week 10-11. No timings
were recorded this week (planned for Week 9).

### Technical decisions

* **Row-major, single pass.** Cells are computed in the order the traceback and the CUDA
  comparison need; no optimisation beyond `-O2` since this is the baseline for the speedup.
* **Full matrix in memory** (`int`, 4 bytes per cell), as planned for traceback: 381.5 MB at 10,000 x 10,000.
* **Size guard.** Sequences longer than 46,000 are rejected so `(M+1) x (N+1) x 4` bytes stays under about 8 GB, and an allocation failure is reported instead of crashing.
* **Scoring** lives in a `Scoring` struct filled from `DEFAULT_MATCH` / `DEFAULT_MISMATCH` / `DEFAULT_GAP` `#define`s, so the CUDA code can use the same values; `--match`, `--mismatch`, `--gap` override them at run time.
* **Checksum.** The program also prints the sum of all matrix cells, so the Week 10 CUDA matrix can be compared with the C matrix on the large inputs without printing them: `bench_1000` 1,563,453; `bench_5000` 1,445,075,095; `bench_10000` 1,009,197,635 (`test_paper` 8, `test_identical` 100, `test_nomatch` 0, `test_gap` 222).

### Limitations / not yet done

* No traceback, so no aligned strings or alignment coordinates yet (Week 9).
* No timing code yet; baseline times come in Week 9.
* `bench_10000` was run in C only. The matrix comparison is limited to inputs up to 40 x 40; larger inputs are compared by score and cell.
* `bench_5000` is only checked when `check_c_week8.py` is run with `--all`.

### Open items

- [ ] Run `check_c_week8.py` on Windows with your compiler and confirm it ends with `All C results match the reference.`
- [ ] Compile and run `sw_seq` on Colab (step 7 of the notebook); the Week 8 results above were produced locally
- [x] Merge the Week 7 and Week 8 branches into `main` (done, both are in `main`)
- [ ] Update Figures 7 and 8 in the proposal (H(4,2) = 1)
- [ ] Move the 23 MB pitch video out of the repository (`git rm --cached "CEPARCO Integrating Project Pitch.mp4"`) and link it here: `<video link>`

### Plan for Week 9

**Deliverable:** verified C program with traceback, and baseline execution times.

* add traceback from the best cell to a 0 (diagonal, then up, then left)
* print the aligned pair and its coordinates; re-score the alignment to confirm it equals the max score
* check against the expected alignments (`test_paper`: `GCT` / `GCT`; `test_gap`: `ACGTCA-TGCATG` / `ACGTCACTGCATG`)
* time the matrix fill on `bench_1000/5000/10000` with repeated runs, and record the machine/GPU session

### Reproducing the Week 8 checks

From `smith-waterman/`:

```bat
mkdir bin
gcc -std=c99 -O2 -Wall -Wextra -o bin\sw_seq.exe c\sw_seq.c         :: MinGW-w64
cl /nologo /O2 /W4 /Fobin\ /Febin\sw_seq.exe c\sw_seq.c             :: or MSVC (x64 Native Tools Command Prompt)
bin\sw_seq.exe data\test_paper.fasta                 :: matrix + max score 3 at H(4, 5) + checksum
bin\sw_seq.exe data\bench_10000.fasta --no-matrix    :: score, cell and checksum only
python scripts\check_c_week8.py                      :: C vs Python reference (add --all for bench_5000)
```

`check_c_week8.py` looks for `bin\sw_seq.exe` (use `--exe PATH` otherwise) and should end with
`All C results match the reference.` Options: `--no-matrix`, `--match N`, `--mismatch N`, `--gap N`
(positive integers, defaults 1 / 1 / 2; the checks use the defaults). Do not commit `bin/`.

---

## Week 7 report (completed)

**Milestone:** finalise scoring parameters and test sequences; set up the Colab notebook and GitHub repository.
**Deliverable:** GitHub repository with README and sample input files.

### What was done

| # | Task | Result |
|---|------|--------|
| 1 | Scoring parameters finalised | match **+1**, mismatch **-1**, gap **-2** (linear). These are the proposal's worked-example values (Section III), so the demo can reproduce that example exactly. |
| 2 | Comparison conventions fixed | Same rules for C and CUDA so outputs can be compared (listed below). |
| 3 | Test sequences prepared | 4 small hand-checkable inputs and 3 benchmark inputs (1,000 / 5,000 / 10,000 bases, the sizes in the proposal) in `data/`. |
| 4 | Expected results derived and independently checked | Hand-derived values for the 4 small inputs, confirmed by an independent pure-Python implementation (`scripts/check_expected.py`) written directly from the recurrence. |
| 5 | Benchmark inputs made reproducible | `scripts/gen_sequences.py` regenerates all three files from a fixed seed; output was compared byte-for-byte with the committed files. |
| 6 | Colab notebook set up | `colab/SW_Week7_Setup.ipynb`: GPU / `nvcc` check, CUDA smoke test, clone, validation of the inputs, environment record. |
| 7 | GitHub repository set up | Branch `build/week7-requirements` with README, sample inputs, scripts and notebook. |

### Technical decisions

* **Scoring:** `H(i,j) = max(0, H(i-1,j-1)+s, H(i-1,j)-d, H(i,j-1)-d)`, with `s = +1` on a match and `-1` on a mismatch, `d = 2`, and row 0 / column 0 equal to 0 (Eq. 1 and Eq. 2 of the proposal). Run time does not depend on these values, so they are not a performance variable.
* **Matrix:** `(M+1) x (N+1)`; sequence *a* down the rows, sequence *b* across the columns.
* **Best cell:** the *first* cell in row-major order holding the maximum score, so C and CUDA report the same cell.
* **Traceback:** from the best cell until a 0 is reached; at each cell try diagonal, then up (gap in *b*), then left (gap in *a*), as in the proposal's pseudocode.
* **Validity of an alignment:** several optimal alignments can exist, so results are compared by score, coordinates, and by re-scoring the reconstructed alignment (Section IV of the proposal), not by string equality alone.
* **Data types:** 32-bit scores, full matrix kept for traceback; a 10,000 x 10,000 input needs about 380 MB.
* **Input:** FASTA with exactly two records, letters only; sequence 1 runs down the rows.

### Test sequences and expected results

| File | Purpose | Expected (+1 / -1 / -2) |
|------|---------|-------------------------|
| `test_paper.fasta` | the proposal's worked example, `AGCT` (rows) vs `ATGCT` (columns) | score **3** at H(4,5); `GCT` / `GCT`; a[2..4], b[3..5] |
| `test_identical.fasta` | identical 8-mers | score **8** at H(8,8); whole sequences align |
| `test_nomatch.fasta` | `AAAAAAAA` vs `CCCCCCCC` | score **0**; no local alignment |
| `test_gap.fasta` | `ACGTCATGCATG` vs the same with one extra `C` | score **10** at H(12,13); `ACGTCA-TGCATG` / `ACGTCACTGCATG` (12 matches, 1 gap) |
| `bench_1000/5000/10000.fasta` | random DNA with a planted, mutated shared region (20% of the length, 10% mutations) | scored by the C baseline in Week 9 |

Expected score matrix for `test_paper.fasta`:

```
        0  1  2  3  4  5
           A  T  G  C  T
  0     0  0  0  0  0  0
  1 A   0  1  0  0  0  0
  2 G   0  0  0  1  0  0
  3 C   0  0  0  0  2  0
  4 T   0  0  1  0  0  3
```

### Verification

| Check | Where | Result |
|-------|-------|--------|
| 4 small inputs vs hand-derived scores, coordinates, alignments, and the full worked-example matrix (`check_expected.py`) | local, Python 3 | pass |
| `check_expected.py` fails when a parameter is changed (gap 2 -> 1) | local | fails (as intended) |
| 3 benchmark files regenerate identically from seed 1 | Google Colab | pass (`identical` x 3) |
| `check_expected.py` on Colab | Google Colab | **to confirm** |
| CUDA smoke test prints `-> OK` | Google Colab | **to confirm** |

### Issues found and how they were resolved

1. **Error in the proposal's Figure 7.** The figure shows 0 at H(4,2); the correct value is **1** (a[4] = T matches b[2] = T and H(3,1) = 0). The maximum score and the `GCT` alignment are unaffected. *Action:* correct Figures 7 and 8 in the proposal; the matrix above is the corrected one.
2. **Colab notebook bug.** The benchmark-regeneration cell used a `for` loop with `!` shell commands; Colab did not expand the loop variable, so the generator received an empty `--length` and every file was wrongly reported as `DIFFERENT`. *Fix:* the cell is now a `%%bash` cell; all three files report `identical`.
3. **Repository layout.** The project folder was copied into the repository as `smith-waterman/smith-waterman/`, so `scripts/check_expected.py` was not found from the clone root. *Fix:* the notebook now locates the project root automatically. *Recommended:* flatten the layout so `data/`, `scripts/`, `colab/` sit at the repository root.

### Open items to close Week 7

- [x] Colab environment recorded (12 Nov 2026): Tesla T4 (15360 MiB), driver 580.82.07, CUDA 13.0 (`nvcc` 13.0.88), Python 3.13.16, gcc 13.3.0
- [x] `check_expected.py` prints `All checks passed.` and the CUDA smoke test prints `-> OK` on Colab
- [x] Merge `build/week7-requirements` into `main` (done)
- [ ] Update Figures 7 and 8 in the proposal (H(4,2) = 1)
- [ ] Optional: flatten the repository layout

### Risks and notes

* **Colab GPU varies** between sessions (usually a T4). The GPU model is recorded with every set of timings, and all final timings come from a single session so C and CUDA run on the same machine.
* **Memory:** the full matrix is stored for traceback (about 380 MB at 10,000 x 10,000), so larger inputs are limited by GPU memory.
* **Baseline is correct and timed.** The C score/cell results are verified (Week 8) and the baseline times were recorded on Colab in Week 9 (see the Week 9 report).
* **Speedup is not yet known.** One kernel launch per anti-diagonal adds launch overhead; the tiled wavefront approach is the planned fallback in Week 12 if profiling shows it limits speedup.

## Repository contents (as of Week 8)

```
data/                      test_*.fasta (hand-checkable) and bench_*.fasta (1k / 5k / 10k)
scripts/gen_sequences.py   reproducible generator for the benchmark inputs
scripts/check_expected.py  independent Python check of the expected results
scripts/check_c_week8.py   compares the C program with the Python reference (score, cell, matrix)
scripts/check_c_week9.py   compares the C traceback with the Python reference (alignment, ranges, re-score)
c/sw_seq.c                 sequential C baseline (Week 8: input, matrix filling, score tracking; Week 9: traceback)
scripts/run_baseline.py    records the baseline C times (CSV + markdown summary in results/)
colab/SW_Week7_Setup.ipynb Colab notebook: GPU / nvcc check, clone, validate inputs
colab/SW_Week9_Timing.ipynb Colab notebook: build, correctness checks, sanitizers, baseline timings
results/                   benchmark output goes here (Weeks 9 and 12)
```

## Reproducing the Week 7 checks

```bat
python scripts\check_expected.py                       :: expect: All checks passed.
python scripts\gen_sequences.py --length 10000 --seed 1 --out b10k.fasta
fc b10k.fasta data\bench_10000.fasta
del b10k.fasta
```

Or open `colab/SW_Week7_Setup.ipynb` in Google Colab with a GPU runtime and run all cells.
The Week 8 commands are listed under *Reproducing the Week 8 checks*.

## References

Numbering follows the project proposal.

[1] T. F. Smith and M. S. Waterman, "Identification of common molecular subsequences,"
J. Mol. Biol., vol. 147, no. 1, pp. 195-197, Mar. 1981, doi: 10.1016/0022-2836(81)90087-5.

[2] M. Farrar, "Striped Smith-Waterman speeds database searches six times over other SIMD
implementations," Bioinformatics, vol. 23, no. 2, pp. 156-161, Jan. 2007,
doi: 10.1093/bioinformatics/btl582.

[3] A. Khajeh-Saeed, S. Poole, and J. B. Perot, "Acceleration of the Smith-Waterman algorithm
using single and multiple graphics processors," J. Comput. Phys., vol. 229, no. 11,
pp. 4247-4258, Jun. 2010, doi: 10.1016/j.jcp.2010.02.009.

[4] A. Coghlan, "The Smith-Waterman algorithm," Wellcome Trust Sanger Institute, SlideShare.
[Online]. Available: https://www.slideshare.net/slideshow/the-smith-waterman-algorithm/16572305

[5] B. Schmidt, F. Kallenborn, A. Chacon, and C. Hundt, "CUDASW++4.0: Ultra-fast GPU-based
Smith-Waterman protein sequence database search," BMC Bioinformatics, vol. 25, Art. no. 342,
Nov. 2024, doi: 10.1186/s12859-024-05965-6.