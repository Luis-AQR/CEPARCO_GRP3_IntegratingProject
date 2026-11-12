/*
 * sw_seq.c - Sequential Smith-Waterman (Week 8)
 * CEPARCO Integrating Project, Group 3
 *
 * Week 8 scope: read two sequences from a FASTA file, fill the scoring matrix
 * row by row, track the maximum score and its cell. Traceback comes in Week 9.
 *
 * Scoring (final, same as the proposal's worked example):
 *     match +1, mismatch -1, gap -2 (linear)
 *
 * Conventions (from the Week 7 README, the CUDA version must follow them too):
 *   - Matrix is (M+1) x (N+1); sequence a (length M) runs down the rows,
 *     sequence b (length N) runs across the columns; row 0 and column 0 are 0.
 *   - Best cell = the FIRST cell in row-major order that holds the maximum
 *     score (so the update uses a strict '>').
 *   - 32-bit integer scores.
 *
 * Build:  gcc -O2 -Wall -Wextra -o sw_seq sw_seq.c
 * Run:    ./sw_seq ../data/test_paper.fasta
 *         ./sw_seq ../data/bench_1000.fasta --no-matrix
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/*scoring parameters*/
#define MATCH_SCORE     1    /* m  */
#define MISMATCH_PEN    1    /* mu (subtracted) */
#define GAP_PEN         2    /* d  (subtracted, once per gap column) */

/* Print the whole matrix only when it is small enough to read. */
#define PRINT_LIMIT     40

static int max4(int a, int b, int c, int d)
{
    int m = a;
    if (b > m) m = b;
    if (c > m) m = c;
    if (d > m) m = d;
    return m;
}

/* Matrix cell (i, j) in a flat array of (M+1) x (N+1) ints. */
#define IDX(i, j, ncols) ((size_t)(i) * (size_t)(ncols) + (size_t)(j))

/*input reading*/

/*
 * Reads a FASTA file with exactly two records and returns the two sequences
 * as newly allocated, upper-case, NUL-terminated strings
 * Returns 0 on success, -1 on any error (message printed to stderr)
 * Handles multi-line records, blank lines and CRLF line endings
 */
static int read_fasta_pair(const char *path, char **out_a, char **out_b)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "error: cannot open '%s'\n", path);
        return -1;
    }

    char *seq[2] = { NULL, NULL };
    size_t len[2] = { 0, 0 }, cap[2] = { 0, 0 };
    int nrec = 0;                 /* records seen so far */
    int cur = -1;                 /* index of the record being filled */

    /* portable getline replacement: a line buffer that grows with fgets */
    size_t buf_size = 4096;
    char *line = (char *)malloc(buf_size);
    if (!line) { fclose(f); return -1; }

    for (;;) {
        size_t used = 0;
        int got_any = 0;
        for (;;) {
            if (used + 2 > buf_size) {
                buf_size *= 2;
                char *tmp = (char *)realloc(line, buf_size);
                if (!tmp) { free(line); fclose(f); return -1; }
                line = tmp;
            }
            if (!fgets(line + used, (int)(buf_size - used), f)) break;
            got_any = 1;
            used += strlen(line + used);
            if (used > 0 && line[used - 1] == '\n') break;
        }
        if (!got_any) break;                       /* end of file */

        /* strip trailing CR / LF / spaces */
        while (used > 0 && (line[used - 1] == '\n' || line[used - 1] == '\r' ||
                            line[used - 1] == ' '  || line[used - 1] == '\t'))
            line[--used] = '\0';
        if (used == 0) continue;                   /* blank line */

        if (line[0] == '>') {                      /* header: start a record */
            if (nrec == 2) {
                fprintf(stderr, "error: more than two records in '%s'\n", path);
                goto fail;
            }
            cur = nrec++;
            continue;
        }
        if (cur < 0) {
            fprintf(stderr, "error: sequence data before the first '>' header\n");
            goto fail;
        }
        for (size_t k = 0; k < used; k++) {
            unsigned char ch = (unsigned char)line[k];
            if (!isalpha(ch)) {
                fprintf(stderr, "error: non-letter character '%c' in record %d\n",
                        ch, cur + 1);
                goto fail;
            }
            if (len[cur] + 2 > cap[cur]) {
                cap[cur] = cap[cur] ? cap[cur] * 2 : 1024;
                char *tmp = (char *)realloc(seq[cur], cap[cur]);
                if (!tmp) goto fail;
                seq[cur] = tmp;
            }
            seq[cur][len[cur]++] = (char)toupper(ch);
            seq[cur][len[cur]] = '\0';
        }
    }

    if (nrec != 2 || len[0] == 0 || len[1] == 0) {
        fprintf(stderr, "error: expected exactly two non-empty records, found %d\n", nrec);
        goto fail;
    }

    free(line);
    fclose(f);
    *out_a = seq[0];
    *out_b = seq[1];
    return 0;

fail:
    free(line);
    free(seq[0]);
    free(seq[1]);
    fclose(f);
    return -1;
}

/*matrix filling + score tracking*/

/*
 * Fills H (already zeroed by calloc, so row 0 and column 0 are the 0s required
 * by the initialisation step) and returns the maximum score through
 * best_score / best_i / best_j.
 */
static void fill_matrix(const char *a, int M, const char *b, int N,
                        int *H, int *best_score, int *best_i, int *best_j)
{
    const int ncols = N + 1;
    int maxScore = 0, maxI = 0, maxJ = 0;

    for (int i = 1; i <= M; i++) {
        for (int j = 1; j <= N; j++) {
            int s    = (a[i - 1] == b[j - 1]) ? MATCH_SCORE : -MISMATCH_PEN;
            int diag = H[IDX(i - 1, j - 1, ncols)] + s;          /* match / mismatch */
            int up   = H[IDX(i - 1, j,     ncols)] - GAP_PEN;    /* gap in b */
            int left = H[IDX(i,     j - 1, ncols)] - GAP_PEN;    /* gap in a */

            int h = max4(0, diag, up, left);
            H[IDX(i, j, ncols)] = h;

            if (h > maxScore) {           
                maxScore = h;
                maxI = i;
                maxJ = j;
            }
        }
    }
    *best_score = maxScore;
    *best_i = maxI;
    *best_j = maxJ;
}

/*output*/

static void print_matrix(const char *a, int M, const char *b, int N, const int *H)
{
    const int ncols = N + 1;
    printf("\nScore matrix H (%d x %d):\n\n", M + 1, N + 1);

    printf("      ");
    for (int j = 0; j <= N; j++) printf("%4d", j);
    printf("\n      ");
    printf("%4s", "");
    for (int j = 1; j <= N; j++) printf("%4c", b[j - 1]);
    printf("\n");

    for (int i = 0; i <= M; i++) {
        printf("%3d %c ", i, i == 0 ? ' ' : a[i - 1]);
        for (int j = 0; j <= N; j++) printf("%4d", H[IDX(i, j, ncols)]);
        printf("\n");
    }
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "usage: %s <input.fasta> [--no-matrix]\n", argv[0]);
        return 2;
    }
    int show_matrix = 1;
    for (int k = 2; k < argc; k++)
        if (strcmp(argv[k], "--no-matrix") == 0) show_matrix = 0;

    char *a = NULL, *b = NULL;
    if (read_fasta_pair(argv[1], &a, &b) != 0) return 1;

    size_t Ms = strlen(a), Ns = strlen(b);
    if (Ms > 46000 || Ns > 46000) {   /* keeps (M+1)*(N+1)*4 bytes under 8 GB */
        fprintf(stderr, "error: sequences too long for this program\n");
        free(a); free(b);
        return 1;
    }
    int M = (int)Ms, N = (int)Ns;

    size_t cells = (size_t)(M + 1) * (size_t)(N + 1);
    int *H = (int *)calloc(cells, sizeof(int)); 
    if (!H) {
        fprintf(stderr, "error: cannot allocate %.1f MB for the matrix\n",
                (double)cells * sizeof(int) / (1024.0 * 1024.0));
        free(a); free(b);
        return 1;
    }

    int best_score, best_i, best_j;
    fill_matrix(a, M, b, N, H, &best_score, &best_i, &best_j);

    printf("Smith-Waterman (sequential C)\n");
    printf("File     : %s\n", argv[1]);
    printf("Lengths  : M = %d (rows), N = %d (columns)\n", M, N);
    printf("Scoring  : match +%d, mismatch -%d, gap -%d\n",
           MATCH_SCORE, MISMATCH_PEN, GAP_PEN);
    printf("Matrix   : %zu cells, %.1f MB\n", cells,
           (double)cells * sizeof(int) / (1024.0 * 1024.0));

    if (show_matrix && M <= PRINT_LIMIT && N <= PRINT_LIMIT)
        print_matrix(a, M, b, N, H);
    else if (show_matrix)
        printf("(matrix not printed: larger than %d x %d)\n", PRINT_LIMIT, PRINT_LIMIT);

    printf("\nMax score: %d\n", best_score);
    if (best_score > 0)
        printf("Max cell : H(%d, %d)\n", best_i, best_j);
    else
        printf("Max cell : none (no local alignment, every cell is 0)\n");

    free(H);
    free(a);
    free(b);
    return 0;
}
