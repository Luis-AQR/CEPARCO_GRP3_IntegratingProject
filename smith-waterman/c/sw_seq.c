/*
 * sw_seq.c - Sequential Smith-Waterman (Week 8 fill + Week 9 traceback)
 * CEPARCO Integrating Project, Group 3
 *
 * Week 8 scope: read two sequences from a FASTA file, fill the scoring matrix
 * row by row, track the maximum score and its cell, print the matrix (small
 * inputs) and the result.
 * Week 9 scope: traceback from the best cell, print the aligned pair with its
 * coordinates, and re-score the alignment independently of the matrix.
 * Traceback priority (shared with the CUDA version): diagonal, then up
 * (gap in b), then left (gap in a); stops at the first 0.
 *
 * Default scoring (final, same as the proposal's worked example):
 *     match +1, mismatch -1, gap -2 (linear)
 * The defaults can be overridden at run time (--match / --mismatch / --gap),
 * because the proposal lists the scoring parameters as inputs. The Python
 * checks and all documented results use the defaults.
 *
 * Conventions (from the Week 7 README, the CUDA version must follow them too):
 *   - Matrix is (M+1) x (N+1); sequence a (length M) runs down the rows,
 *     sequence b (length N) runs across the columns; row 0 and column 0 are 0.
 *   - Best cell = the FIRST cell in row-major order that holds the maximum
 *     score (so the update uses a strict '>').
 *   - 32-bit integer scores.
 *
 * Plain standard C (C99): builds unchanged with MSVC (cl), MinGW-w64 gcc, or
 * any other C compiler. No POSIX functions are used.
 *
 * Build (Windows, from smith-waterman/, one of):
 *     gcc -std=c99 -O2 -Wall -Wextra -o bin\sw_seq.exe c\sw_seq.c      (MinGW-w64)
 *     cl /nologo /O2 /W4 /Fobin\ /Febin\sw_seq.exe c\sw_seq.c          (MSVC)
 * Run:
 *     bin\sw_seq.exe data\test_paper.fasta
 *     bin\sw_seq.exe data\bench_1000.fasta --no-matrix
 *     bin\sw_seq.exe data\test_paper.fasta --match 2 --mismatch 1 --gap 3
 */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS /* MSVC: allow fopen/strcpy without warnings */
#endif

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Default scoring parameters */
#define DEFAULT_MATCH 1    /* m  */
#define DEFAULT_MISMATCH 1 /* mu (subtracted) */
#define DEFAULT_GAP 2      /* d  (subtracted, once per gap column) */

/* Print the whole matrix only when it is small enough to read. */
#define PRINT_LIMIT 40

/* Longest sequence accepted: (46001)^2 cells * 4 bytes is about 8.5 GB. */
#define MAX_SEQ_LEN 46000

typedef struct {
    int match;    /* added when a[i] == b[j]        */
    int mismatch; /* subtracted when a[i] != b[j]   */
    int gap;      /* subtracted for each gap column */
} Scoring;

typedef struct {
    int score;                   /* maximum value in the matrix            */
    int i, j;                    /* first cell (row-major) holding it      */
    unsigned long long checksum; /* sum of every cell, for C/CUDA compares */
} FillResult;

typedef struct {
    char *buf_a;          /* owned buffers (free with free_alignment)          */
    char *buf_b;
    char *a_aln;          /* aligned a, '-' = gap, NUL-terminated (into buf_a) */
    char *b_aln;          /* aligned b, same length as a_aln (into buf_b)      */
    int len;              /* number of alignment columns                       */
    int i_start, i_end;   /* 1-based inclusive range in a                      */
    int j_start, j_end;   /* 1-based inclusive range in b                      */
    int matches, mismatches, gaps;
} Alignment;

/* Width of one block of the printed alignment. */
#define ALIGN_LINE 60

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

/* ------------------------------------------------------------------ */
/* Input reading                                                      */
/* ------------------------------------------------------------------ */

/*
 * Reads a FASTA file with exactly two records and returns the two sequences
 * as newly allocated, upper-case, NUL-terminated strings.
 * Returns 0 on success, -1 on any error (message printed to stderr).
 * Handles multi-line records, blank lines and CRLF line endings.
 */
static int read_fasta_pair(const char *path, char **out_a, char **out_b)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "error: cannot open '%s'\n", path);
        return -1;
    }

    char *seq[2] = {NULL, NULL};
    size_t len[2] = {0, 0}, cap[2] = {0, 0};
    int nrec = 0; /* records seen so far */
    int cur = -1; /* index of the record being filled */

    /* portable getline replacement: a line buffer that grows with fgets */
    size_t buf_size = 4096;
    char *line = (char *)malloc(buf_size);
    if (!line) {
        fclose(f);
        return -1;
    }

    for (;;) {
        size_t used = 0;
        int got_any = 0;
        for (;;) {
            if (used + 2 > buf_size) {
                buf_size *= 2;
                char *tmp = (char *)realloc(line, buf_size);
                if (!tmp) {
                    free(line);
                    free(seq[0]);
                    free(seq[1]);
                    fclose(f);
                    return -1;
                }
                line = tmp;
            }
            if (!fgets(line + used, (int)(buf_size - used), f)) break;
            got_any = 1;
            used += strlen(line + used);
            if (used > 0 && line[used - 1] == '\n') break;
        }
        if (!got_any) break; /* end of file */

        /* strip trailing CR / LF / spaces */
        while (used > 0 && (line[used - 1] == '\n' || line[used - 1] == '\r' ||
                            line[used - 1] == ' ' || line[used - 1] == '\t'))
            line[--used] = '\0';
        if (used == 0) continue; /* blank line */

        if (line[0] == '>') { /* header: start a record */
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

/* ------------------------------------------------------------------ */
/* Matrix filling + score tracking                                    */
/* ------------------------------------------------------------------ */

/*
 * Fills H (already zeroed by calloc, so row 0 and column 0 are the 0s required
 * by the initialisation step) row by row, left to right, and reports the
 * maximum score, its cell and a checksum of the whole matrix.
 */
static FillResult fill_matrix(const char *a, int M, const char *b, int N,
                              int *H, Scoring sc)
{
    const int ncols = N + 1;
    FillResult r;
    r.score = 0;
    r.i = 0;
    r.j = 0;
    r.checksum = 0;

    for (int i = 1; i <= M; i++) {
        for (int j = 1; j <= N; j++) {
            int s = (a[i - 1] == b[j - 1]) ? sc.match : -sc.mismatch;
            int diag = H[IDX(i - 1, j - 1, ncols)] + s;      /* match / mismatch */
            int up = H[IDX(i - 1, j, ncols)] - sc.gap;       /* gap in b */
            int left = H[IDX(i, j - 1, ncols)] - sc.gap;     /* gap in a */

            int h = max4(0, diag, up, left);
            H[IDX(i, j, ncols)] = h;
            r.checksum += (unsigned long long)h;

            if (h > r.score) { /* strict '>': keeps the FIRST maximum */
                r.score = h;
                r.i = i;
                r.j = j;
            }
        }
    }
    return r;
}

/* ------------------------------------------------------------------ */
/* Traceback + independent re-score                                   */
/* ------------------------------------------------------------------ */

static void free_alignment(Alignment *al)
{
    free(al->buf_a);
    free(al->buf_b);
    al->buf_a = al->buf_b = al->a_aln = al->b_aln = NULL;
}

/*
 * Traceback from the best cell (bi, bj) until a 0 is reached. At each cell the
 * moves are tried in a fixed order: diagonal, up (gap in b), left (gap in a).
 * The moves are recomputed from H and the sequences (no direction matrix).
 * The aligned strings are built backwards, so no reverse pass is needed.
 * Returns 0 on success, -1 on allocation failure or an inconsistent matrix.
 */
static int traceback(const char *a, const char *b, int M, int N, const int *H,
                     Scoring sc, int bi, int bj, Alignment *out)
{
    const int ncols = N + 1;
    const size_t cap = (size_t)M + (size_t)N + 1; /* longest possible alignment + 1 */

    memset(out, 0, sizeof *out);
    out->buf_a = (char *)malloc(cap + 1);
    out->buf_b = (char *)malloc(cap + 1);
    if (!out->buf_a || !out->buf_b) {
        fprintf(stderr, "error: cannot allocate the alignment buffers\n");
        free_alignment(out);
        return -1;
    }

    size_t pos = cap;
    out->buf_a[pos] = '\0';
    out->buf_b[pos] = '\0';

    int i = bi, j = bj;
    out->i_end = bi;
    out->j_end = bj;

    while (H[IDX(i, j, ncols)] > 0) {
        const int h = H[IDX(i, j, ncols)];

        if (i > 0 && j > 0) {
            int s = (a[i - 1] == b[j - 1]) ? sc.match : -sc.mismatch;
            if (h == H[IDX(i - 1, j - 1, ncols)] + s) { /* 1. diagonal */
                --pos;
                out->buf_a[pos] = a[i - 1];
                out->buf_b[pos] = b[j - 1];
                if (a[i - 1] == b[j - 1]) out->matches++;
                else out->mismatches++;
                i--;
                j--;
                continue;
            }
        }
        if (i > 0 && h == H[IDX(i - 1, j, ncols)] - sc.gap) { /* 2. up: gap in b */
            --pos;
            out->buf_a[pos] = a[i - 1];
            out->buf_b[pos] = '-';
            out->gaps++;
            i--;
            continue;
        }
        if (j > 0 && h == H[IDX(i, j - 1, ncols)] - sc.gap) { /* 3. left: gap in a */
            --pos;
            out->buf_a[pos] = '-';
            out->buf_b[pos] = b[j - 1];
            out->gaps++;
            j--;
            continue;
        }

        fprintf(stderr, "error: traceback stuck at H(%d, %d) = %d: no move explains it\n",
                i, j, h);
        free_alignment(out);
        return -1;
    }

    out->i_start = i + 1; /* the cell where the walk stopped (a 0) is not aligned */
    out->j_start = j + 1;
    out->len = (int)(cap - pos);
    out->a_aln = out->buf_a + pos;
    out->b_aln = out->buf_b + pos;
    return 0;
}

/*
 * Scores the alignment from the two aligned strings only (H is not used) and
 * checks that the strings are well formed: same length, no gap/gap column, and
 * the un-gapped strings equal a[i_start..i_end] and b[j_start..j_end].
 * Returns the score, or INT_MIN if the alignment is malformed.
 */
static int rescore_alignment(const Alignment *al, const char *a, const char *b, Scoring sc)
{
    if (!al->a_aln || !al->b_aln) return INT_MIN;
    if ((int)strlen(al->a_aln) != al->len || (int)strlen(al->b_aln) != al->len) return INT_MIN;

    int score = 0;
    int ia = al->i_start - 1; /* 0-based index of the next residue of a */
    int jb = al->j_start - 1;
    for (int k = 0; k < al->len; k++) {
        char x = al->a_aln[k], y = al->b_aln[k];
        if (x == '-' && y == '-') return INT_MIN;
        if (x == '-') {
            if (b[jb] != y) return INT_MIN;
            jb++;
            score -= sc.gap;
        } else if (y == '-') {
            if (a[ia] != x) return INT_MIN;
            ia++;
            score -= sc.gap;
        } else {
            if (a[ia] != x || b[jb] != y) return INT_MIN;
            ia++;
            jb++;
            score += (x == y) ? sc.match : -sc.mismatch;
        }
    }
    /* the walk must end exactly at the recorded end cell */
    if (ia != al->i_end || jb != al->j_end) return INT_MIN;
    return score;
}

/* ------------------------------------------------------------------ */
/* Output                                                             */
/* ------------------------------------------------------------------ */

/*
 * Prints the alignment in blocks of ALIGN_LINE columns. Each block has an
 * "a" line, a match-bar line and a "b" line; the number is the 1-based index
 * of the next residue of that sequence at the start of the block.
 */
static void print_alignment(const Alignment *al)
{
    printf("Alignment: length %d (%d matches, %d mismatches, %d gaps)\n", al->len,
           al->matches, al->mismatches, al->gaps);
    printf("Range    : a[%d..%d], b[%d..%d]\n", al->i_start, al->i_end, al->j_start,
           al->j_end);

    int ia = al->i_start, jb = al->j_start;
    for (int k = 0; k < al->len; k += ALIGN_LINE) {
        int n = al->len - k < ALIGN_LINE ? al->len - k : ALIGN_LINE;
        printf("\n  a %6d  %.*s\n", ia, n, al->a_aln + k);
        printf("            ");
        for (int c = 0; c < n; c++) {
            char x = al->a_aln[k + c], y = al->b_aln[k + c];
            putchar((x == y && x != '-') ? '|' : ' ');
        }
        printf("\n  b %6d  %.*s\n", jb, n, al->b_aln + k);
        for (int c = 0; c < n; c++) {
            if (al->a_aln[k + c] != '-') ia++;
            if (al->b_aln[k + c] != '-') jb++;
        }
    }
    printf("\n");
}

static void print_matrix(const char *b, int N, const char *a, int M, const int *H)
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

static void usage(const char *prog)
{
    fprintf(stderr,
            "usage: %s <input.fasta> [--no-matrix] [--no-traceback] [--match N] [--mismatch N] [--gap N]\n"
            "  FASTA with exactly two records; sequence 1 runs down the rows.\n"
            "  Defaults: match +%d, mismatch -%d, gap -%d (N must be a positive integer).\n",
            prog, DEFAULT_MATCH, DEFAULT_MISMATCH, DEFAULT_GAP);
}

/* Parses a positive int; returns 0 on success, -1 otherwise. */
static int parse_positive(const char *text, int *out)
{
    char *end = NULL;
    errno = 0;
    long v = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || v < 1 || v > 1000000) return -1;
    *out = (int)v;
    return 0;
}

int main(int argc, char **argv)
{
    if (argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) {
        usage(argv[0]);
        return 2;
    }

    Scoring sc = {DEFAULT_MATCH, DEFAULT_MISMATCH, DEFAULT_GAP};
    int show_matrix = 1;
    int do_traceback = 1;
    for (int k = 2; k < argc; k++) {
        int *target = NULL;
        if (strcmp(argv[k], "--no-matrix") == 0) {
            show_matrix = 0;
            continue;
        } else if (strcmp(argv[k], "--no-traceback") == 0) {
            do_traceback = 0;
            continue;
        } else if (strcmp(argv[k], "--match") == 0) {
            target = &sc.match;
        } else if (strcmp(argv[k], "--mismatch") == 0) {
            target = &sc.mismatch;
        } else if (strcmp(argv[k], "--gap") == 0) {
            target = &sc.gap;
        } else {
            fprintf(stderr, "error: unknown option '%s'\n", argv[k]);
            usage(argv[0]);
            return 2;
        }
        if (k + 1 >= argc || parse_positive(argv[k + 1], target) != 0) {
            fprintf(stderr, "error: %s needs a positive integer\n", argv[k]);
            return 2;
        }
        k++;
    }

    char *a = NULL, *b = NULL;
    if (read_fasta_pair(argv[1], &a, &b) != 0) return 1;

    size_t Ms = strlen(a), Ns = strlen(b);
    if (Ms > MAX_SEQ_LEN || Ns > MAX_SEQ_LEN) {
        fprintf(stderr, "error: sequences longer than %d are not supported\n", MAX_SEQ_LEN);
        free(a);
        free(b);
        return 1;
    }
    int M = (int)Ms, N = (int)Ns;

    size_t cells = (size_t)(M + 1) * (size_t)(N + 1);
    double mb = (double)cells * sizeof(int) / (1024.0 * 1024.0);
    int *H = (int *)calloc(cells, sizeof(int));
    if (!H) {
        fprintf(stderr, "error: cannot allocate %.1f MB for the matrix\n", mb);
        free(a);
        free(b);
        return 1;
    }

    FillResult res = fill_matrix(a, M, b, N, H, sc);

    printf("Smith-Waterman (sequential C)\n");
    printf("File     : %s\n", argv[1]);
    printf("Lengths  : M = %d (rows), N = %d (columns)\n", M, N);
    printf("Scoring  : match +%d, mismatch -%d, gap -%d\n", sc.match, sc.mismatch, sc.gap);
    printf("Matrix   : %llu cells, %.1f MB\n", (unsigned long long)cells, mb);

    if (show_matrix && M <= PRINT_LIMIT && N <= PRINT_LIMIT)
        print_matrix(b, N, a, M, H);
    else if (show_matrix)
        printf("(matrix not printed: larger than %d x %d)\n", PRINT_LIMIT, PRINT_LIMIT);

    printf("\nMax score: %d\n", res.score);
    if (res.score > 0)
        printf("Max cell : H(%d, %d)\n", res.i, res.j);
    else
        printf("Max cell : none (no local alignment, every cell is 0)\n");

    int status = 0;
    if (res.score > 0 && do_traceback) {
        Alignment al;
        if (traceback(a, b, M, N, H, sc, res.i, res.j, &al) != 0) {
            status = 1;
        } else {
            printf("\n");
            print_alignment(&al);
            int rs = rescore_alignment(&al, a, b, sc);
            if (rs == res.score) {
                printf("Re-score : %d (OK)\n", rs);
            } else {
                if (rs == INT_MIN)
                    printf("Re-score : malformed alignment (MISMATCH, expected %d)\n", res.score);
                else
                    printf("Re-score : %d (MISMATCH, expected %d)\n", rs, res.score);
                status = 1;
            }
            free_alignment(&al);
        }
    } else if (res.score == 0 && do_traceback) {
        printf("\nNo local alignment (every cell is 0), traceback skipped.\n");
    }

    printf("Checksum : %llu (sum of all matrix cells)\n", res.checksum);

    free(H);
    free(a);
    free(b);
    return status;
}