"""
run_baseline.py - Week 9: record the baseline (sequential C) execution times.

Official timings are taken on Colab, in ONE session, so that everyone on the team uses the
same environment and the CUDA timings can later be taken on the same machine.

For each benchmark input it runs sw_seq twice, each with one warm-up run plus N timed runs:
  * mode "full":      fill + traceback (the normal program)
  * mode "fill_only": --no-traceback (shows the traceback cost is negligible)
and writes, in --out-dir (default: results/):
  * week9_baseline.csv   one row per run, with the environment columns
  * week9_baseline.md    environment record + summary table (min / median / max)

Usage (from smith-waterman/):
    python3 scripts/run_baseline.py --build              # gcc -O2 build, then time
    python3 scripts/run_baseline.py --exe bin/sw_seq     # use an existing build
    python3 scripts/run_baseline.py --repeat 5 --repeat-large 3
    python3 scripts/run_baseline.py --out-dir /tmp/test  # do not touch results/ (testing)
"""
import argparse
import csv
import datetime
import os
import platform
import re
import statistics
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DATA = os.path.join(ROOT, "data")
INPUTS = [("bench_1000", 1000), ("bench_5000", 5000), ("bench_10000", 10000)]
CFLAGS = ["-std=c99", "-O2", "-Wall", "-Wextra"]


def sh(cmd):
    """Run a command, return stripped stdout or '' on any failure."""
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=30).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return ""


def cpu_model():
    out = sh(["lscpu"])
    m = re.search(r"Model name:\s*(.+)", out)
    if m:
        return m.group(1).strip()
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def environment(compiler, flags):
    gpu = sh(["nvidia-smi", "--query-gpu=name,memory.total,driver_version",
              "--format=csv,noheader"]) or "none detected"
    cc_ver = sh([compiler, "--version"]).splitlines()
    return {
        "date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "platform": platform.platform(),
        "cpu": cpu_model(),
        "logical_cpus": str(os.cpu_count()),
        "gpu": gpu.replace("\n", "; "),
        "compiler": cc_ver[0] if cc_ver else "unknown",
        "flags": " ".join(flags),
        "python": platform.python_version(),
    }


def build(compiler, out_path):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    src = os.path.join(ROOT, "c", "sw_seq.c")
    cmd = [compiler] + CFLAGS + ["-o", out_path, src]
    print("building:", " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit("build failed:\n" + r.stdout + r.stderr)
    if r.stderr.strip():
        print(r.stderr)


def run_once(exe, name, repeat, traceback_on):
    args = [exe, os.path.join(DATA, name + ".fasta"), "--no-matrix",
            "--repeat", str(repeat), "--warmup", "--csv"]
    if not traceback_on:
        args.append("--no-traceback")
    p = subprocess.run(args, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"{name}: sw_seq failed (exit {p.returncode})\n{p.stdout[-500:]}\n{p.stderr}")
    rows = []
    for line in p.stdout.splitlines():
        if line.startswith("CSV,") and not line.startswith("CSV,input,"):
            f = line.split(",")
            rows.append({
                "input": f[1], "M": int(f[2]), "N": int(f[3]), "run": int(f[4]),
                "alloc_s": float(f[5]), "fill_s": float(f[6]),
                "traceback_s": float(f[7]), "total_s": float(f[8]),
                "score": int(f[9]), "cell_i": int(f[10]), "cell_j": int(f[11]),
                "checksum": int(f[12]),
            })
    if len(rows) != repeat + 1:
        sys.exit(f"{name}: expected {repeat + 1} CSV rows, got {len(rows)}")
    return rows


# Reference values from the Week 8 / Week 9 reports (C baseline, scoring +1/-1/-2).
EXPECTED = {
    "bench_1000": (160, 784, 995, 1563453),
    "bench_5000": (805, 2313, 1448, 1445075095),
    "bench_10000": (1536, 5403, 9807, 1009197635),
}


def main():
    ap = argparse.ArgumentParser(description="Record baseline C execution times.")
    ap.add_argument("--exe", help="path to the compiled sw_seq (default: bin/sw_seq[.exe])")
    ap.add_argument("--build", action="store_true", help="build bin/sw_seq first")
    ap.add_argument("--cc", default="gcc", help="compiler for --build (default: gcc)")
    ap.add_argument("--repeat", type=int, default=5, help="timed runs for 1k and 5k (default 5)")
    ap.add_argument("--repeat-large", type=int, default=3, help="timed runs for 10k (default 3)")
    ap.add_argument("--out-dir", default=os.path.join(ROOT, "results"))
    args = ap.parse_args()

    exe = args.exe
    if args.build:
        exe = exe or os.path.join(ROOT, "bin", "sw_seq")
        build(args.cc, exe)
    if not exe:
        for c in ("sw_seq", "sw_seq.exe"):
            p = os.path.join(ROOT, "bin", c)
            if os.path.isfile(p):
                exe = p
                break
    if not exe or not os.path.isfile(exe):
        sys.exit("error: compiled program not found. Use --build or --exe PATH.")
    exe = os.path.abspath(exe)

    env = environment(args.cc, CFLAGS)
    print("Environment:")
    for k, v in env.items():
        print(f"  {k:13s}: {v}")
    if "Linux" not in env["platform"]:
        print("\nNOTE: this is not Linux/Colab. These numbers are for development only; "
              "official timings must come from Colab.\n")

    all_rows = []
    for name, _ in INPUTS:
        reps = args.repeat_large if name == "bench_10000" else args.repeat
        for mode, tb in (("full", True), ("fill_only", False)):
            print(f"running {name} ({mode}), warm-up + {reps} runs ...", flush=True)
            for r in run_once(exe, name, reps, tb):
                exp = EXPECTED[name]
                if (r["score"], r["cell_i"], r["cell_j"], r["checksum"]) != exp:
                    sys.exit(f"{name}: result {r['score']} at ({r['cell_i']},{r['cell_j']}) "
                             f"checksum {r['checksum']} differs from the verified baseline {exp}")
                r["mode"] = mode
                r["warmup"] = 1 if r["run"] == 0 else 0
                r.update(env)
                all_rows.append(r)

    os.makedirs(args.out_dir, exist_ok=True)
    cols = ["date", "platform", "cpu", "logical_cpus", "gpu", "compiler", "flags", "input", "M",
            "N", "mode", "run", "warmup", "alloc_s", "fill_s", "traceback_s", "total_s", "score",
            "cell_i", "cell_j", "checksum"]
    csv_path = os.path.join(args.out_dir, "week9_baseline.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(all_rows)

    md_path = os.path.join(args.out_dir, "week9_baseline.md")
    with open(md_path, "w") as f:
        f.write("# Week 9 baseline: sequential C execution times\n\n")
        f.write("Scoring +1 / -1 / -2. Times in seconds. Statistics exclude the warm-up run "
                "(run 0). `alloc` = malloc + zeroing the matrix; reading the file is not timed.\n\n")
        f.write("## Environment\n\n")
        for k, v in env.items():
            f.write(f"- **{k}:** {v}\n")
        f.write("\n## Summary (full = fill + traceback)\n\n")
        f.write("| Input | Runs | Fill median | Fill min-max | Traceback median | Alloc median "
                "| Total median |\n|---|---|---|---|---|---|---|\n")
        for name, _ in INPUTS:
            rs = [r for r in all_rows if r["input"] == name and r["mode"] == "full" and r["run"] > 0]
            fl = [r["fill_s"] for r in rs]
            f.write(f"| {name} | {len(rs)} | {statistics.median(fl):.4f} | "
                    f"{min(fl):.4f} - {max(fl):.4f} | "
                    f"{statistics.median(r['traceback_s'] for r in rs):.6f} | "
                    f"{statistics.median(r['alloc_s'] for r in rs):.4f} | "
                    f"{statistics.median(r['total_s'] for r in rs):.4f} |\n")
        f.write("\n## Fill-only runs (`--no-traceback`), to confirm the traceback is negligible\n\n")
        f.write("| Input | Runs | Fill median | Total median |\n|---|---|---|---|\n")
        for name, _ in INPUTS:
            rs = [r for r in all_rows if r["input"] == name and r["mode"] == "fill_only" and r["run"] > 0]
            f.write(f"| {name} | {len(rs)} | {statistics.median(r['fill_s'] for r in rs):.4f} | "
                    f"{statistics.median(r['total_s'] for r in rs):.4f} |\n")
        f.write("\nEvery run reproduced the verified score, best cell and checksum.\n")

    print(f"\nwrote {csv_path}\nwrote {md_path}")
    with open(md_path) as f:
        print("\n" + f.read())
    return 0


if __name__ == "__main__":
    sys.exit(main())
