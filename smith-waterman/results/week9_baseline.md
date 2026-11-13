# Week 9 baseline: sequential C execution times

Scoring +1 / -1 / -2. Times in seconds. Statistics exclude the warm-up run (run 0). `alloc` = malloc + zeroing the matrix; reading the file is not timed.

## Environment

- **date:** 2026-11-13 12:14
- **platform:** Linux-6.6.159+-x86_64-with-glibc2.39
- **cpu:** Intel(R) Xeon(R) CPU @ 2.00GHz
- **logical_cpus:** 2
- **gpu:** Tesla T4, 15360 MiB, 580.82.07
- **compiler:** gcc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0
- **flags:** -std=c99 -O2 -Wall -Wextra
- **python:** 3.13.16

## Summary (full = fill + traceback)

| Input | Runs | Fill median | Fill min-max | Traceback median | Alloc median | Total median |
|---|---|---|---|---|---|---|
| bench_1000 | 5 | 0.0020 | 0.0020 - 0.0021 | 0.000003 | 0.0001 | 0.0021 |
| bench_5000 | 5 | 0.0495 | 0.0485 - 0.0512 | 0.000049 | 0.0464 | 0.0960 |
| bench_10000 | 3 | 0.2111 | 0.2072 - 0.2159 | 0.000092 | 0.1925 | 0.4086 |

## Fill-only runs (`--no-traceback`), to confirm the traceback is negligible

| Input | Runs | Fill median | Total median |
|---|---|---|---|
| bench_1000 | 5 | 0.0019 | 0.0021 |
| bench_5000 | 5 | 0.0508 | 0.0998 |
| bench_10000 | 3 | 0.1980 | 0.3859 |

Every run reproduced the verified score, best cell and checksum.
