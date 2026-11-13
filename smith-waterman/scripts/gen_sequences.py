"""
gen_sequences.py - generate a reproducible FASTA pair for benchmarking.

Both sequences are random DNA of the requested length.  A "shared" region
(default 20% of the length) is planted in both: sequence b gets a mutated copy
(substitutions, insertions, deletions) of the region found in sequence a.  This
gives the Smith-Waterman run a clear, long local alignment to find and to
reconstruct, while the rest of the sequences stay unrelated.

Usage:
    python scripts\\gen_sequences.py --length 1000 --seed 1 --out data/bench_1000.fasta
"""
import argparse
import random

BASES = "ACGT"


def rand_dna(rng, n):
    return "".join(rng.choice(BASES) for _ in range(n))


def mutate(rng, seq, rate):
    """Substitute / delete / insert with a total per-base probability `rate`."""
    out = []
    for ch in seq:
        r = rng.random()
        if r < rate / 3:                      # substitution
            out.append(rng.choice([b for b in BASES if b != ch]))
        elif r < 2 * rate / 3:                # deletion
            continue
        elif r < rate:                        # insertion before this base
            out.append(rng.choice(BASES))
            out.append(ch)
        else:
            out.append(ch)
    return "".join(out)


def wrap(seq, width=70):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--length", type=int, required=True, help="length of each sequence")
    ap.add_argument("--seed", type=int, default=1, help="random seed (default 1)")
    ap.add_argument("--shared", type=float, default=0.2,
                    help="fraction of the length that is a shared region (default 0.2)")
    ap.add_argument("--mutation", type=float, default=0.10,
                    help="per-base mutation rate inside the shared region (default 0.10)")
    ap.add_argument("--out", required=True, help="output FASTA path")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    n = args.length
    core_len = max(1, int(n * args.shared))
    core_a = rand_dna(rng, core_len)
    core_b = mutate(rng, core_a, args.mutation)

    # sequence a = random + core_a + random
    left_a = rng.randint(0, n - core_len)
    a = rand_dna(rng, left_a) + core_a + rand_dna(rng, n - core_len - left_a)

    # sequence b = random + core_b + random, trimmed/padded to exactly n bases
    core_b = core_b[:n]
    room = n - len(core_b)
    left_b = rng.randint(0, room)
    b = rand_dna(rng, left_b) + core_b + rand_dna(rng, room - left_b)
    assert len(a) == n and len(b) == n

    with open(args.out, "w", newline="\n") as f:   # LF on every OS
        f.write(f">seq_a_len{n}_seed{args.seed}\n{wrap(a)}\n")
        f.write(f">seq_b_len{n}_seed{args.seed}_shared_region_mutated\n{wrap(b)}\n")
    print(f"wrote {args.out}: 2 x {n} bases, shared region ~{core_len} bases, "
          f"{args.mutation:.0%} mutation")


if __name__ == "__main__":
    main()