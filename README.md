# Smith-Waterman: Full Alignment Reconstruction, C vs CUDA
## Weekly Progress Report

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
| **7** | Finalise scoring parameters and test sequences; set up Colab notebook and GitHub repository | **Done** (this report) - 3 confirmations still open, see *Open items* |
| 8 | Sequential C: input reading, matrix filling, score tracking | Not started |
| 9 | C traceback; record baseline execution times | Not started |
| 10 | Mid-project demo: C baseline + working CUDA matrix filling | Not started |
| 11 | CUDA traceback and alignment reconstruction; compare with C | Not started |
| 12 | Larger inputs, `nsys` profiling, speedup; tiled wavefront if launch overhead is high | Not started |
| 13 | Documentation and demo | Not started |

---

## Week 7 report

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

- [ ] Run the last notebook cell and record the Colab environment: GPU `<model>`, driver `<version>`, CUDA `<version>`
- [ ] Confirm `check_expected.py` prints `All checks passed.` and the CUDA smoke test prints `-> OK` on Colab
- [ ] Merge `build/week7-requirements` into `main`
- [ ] Update Figures 7 and 8 in the proposal (H(4,2) = 1)
- [ ] Optional: flatten the repository layout; consider moving the 23 MB pitch video out of the repository, since every Colab clone downloads it

### Plan for Week 8

**Deliverable:** a working C program that outputs the score matrix and the maximum score.

* read and validate the FASTA input
* fill the matrix row by row and track the maximum score and its cell
* print the matrix for small inputs
* check against the expected values in this report (`test_paper`: matrix above, score 3 at H(4,5); `test_identical`: 8; `test_nomatch`: 0; `test_gap`: 10 at H(12,13))

### Risks and notes

* **Colab GPU varies** between sessions (usually a T4). The GPU model is recorded with every set of timings, and all final timings come from a single session so C and CUDA run on the same machine.
* **Memory:** the full matrix is stored for traceback (about 380 MB at 10,000 x 10,000), so larger inputs are limited by GPU memory.
* **Speedup is not yet known.** One kernel launch per anti-diagonal adds launch overhead; the tiled wavefront approach is the planned fallback in Week 12 if profiling shows it limits speedup.

### Team contributions (Week 7)

| Member | Contribution |
|--------|--------------|
| Alcasid, Aizy | |
| Bilanes, Kendrick | |
| Corpuz, Gerald | |
| Roa, Luis | |

---

## Repository contents (as of Week 7)

```
data/                      test_*.fasta (hand-checkable) and bench_*.fasta (1k / 5k / 10k)
scripts/gen_sequences.py   reproducible generator for the benchmark inputs
scripts/check_expected.py  independent Python check of the expected results
colab/SW_Week7_Setup.ipynb Colab notebook: GPU / nvcc check, clone, validate inputs
results/                   benchmark output goes here (Weeks 9 and 12)
```

## Reproducing the Week 7 checks

```bash
python3 scripts/check_expected.py                      # expect: All checks passed.
python3 scripts/gen_sequences.py --length 10000 --seed 1 --out /tmp/b10k.fasta
cmp /tmp/b10k.fasta data/bench_10000.fasta && echo identical
```

Or open `colab/SW_Week7_Setup.ipynb` in Google Colab with a GPU runtime and run all cells.

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
