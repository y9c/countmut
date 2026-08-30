#!/usr/bin/env python3
"""Monitor the RSS (memory) of a running countmut process in real time.

Usage:
    scripts/countmut_monitor.py [--interval 0.2] [--verbose-pass] -- CMD [ARGS...]

Samples the resident set size of the command's whole process tree every
<interval> seconds, reports the peak, and optionally passes the child's stderr
through (so `countmut --verbose` progress lines -- which already carry a
"rss=NNNMB" field -- show up live).  On exit prints the peak and a compact
timeline (MB vs seconds).

Example:
    countmut_monitor.py --interval 0.2 -- \
        countmut_core --bam big.bam --fa ref.fa --mode mutation --threads 8
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import threading
import time


def _children(pid: int):
    out = [pid]
    try:
        with open(f"/proc/{pid}/task/{pid}/children") as f:
            for tok in f.read().split():
                out.extend(_children(int(tok)))
    except (FileNotFoundError, ProcessLookupError, ValueError):
        pass
    return out


def _rss_mb(pid: int) -> int:
    try:
        with open(f"/proc/{pid}/statm") as f:
            _size, resident = map(int, f.read().split()[:2])
        return resident * os.sysconf("SC_PAGE_SIZE") // (1024 * 1024)
    except (FileNotFoundError, ValueError, IndexError, OSError):
        return 0


def _total_rss_mb(pid: int) -> int:
    return sum(_rss_mb(p) for p in _children(pid))


def _forward(stream) -> None:
    try:
        for line in stream:
            sys.stderr.write(line)
            sys.stderr.flush()
    except (ValueError, OSError):
        pass


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--interval", type=float, default=0.2, help="sampling interval (s)")
    ap.add_argument("--top", type=int, default=60, help="max timeline samples printed")
    ap.add_argument("--no-pass", action="store_true", help="do not forward child stderr")
    ap.add_argument("cmd", nargs=argparse.REMAINDER)
    args = ap.parse_args()
    cmd = args.cmd[1:] if args.cmd and args.cmd[0] == "--" else args.cmd
    if not cmd:
        ap.error("no command to run (use `-- CMD ...`)")

    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL,
                            stderr=None if args.no_pass else subprocess.PIPE, text=True)
    if not args.no_pass:
        threading.Thread(target=_forward, args=(proc.stderr,), daemon=True).start()

    t0 = time.time()
    samples: list[tuple[float, int]] = []
    peak = 0
    while proc.poll() is None:
        rss = _total_rss_mb(proc.pid)
        t = time.time() - t0
        peak = max(peak, rss)
        samples.append((t, rss))
        sys.stderr.write(f"[monitor] t={t:6.1f}s rss={rss:6}MB peak={peak:6}MB\n")
        sys.stderr.flush()
        time.sleep(max(0.01, args.interval))
    proc.wait()

    # final flush (child may have few/no samples)
    rss = _total_rss_mb(proc.pid) if proc.returncode is None else 0
    peak = max(peak, rss)
    if proc.returncode not in (None, 0):
        sys.stderr.write(f"[monitor] WARNING: command exited {proc.returncode}\n")

    sys.stderr.write(f"[monitor] peak RSS = {peak} MB ({len(samples)} samples)\n")
    if samples:
        n = min(len(samples), args.top)
        step = max(1, len(samples) // n)
        sys.stderr.write("[monitor] timeline (s:MB)  ")
        sys.stderr.write(" ".join(f"{samples[i][0]:.0f}:{samples[i][1]}" for i in range(0, len(samples), step)))
        sys.stderr.write("\n")
    sys.exit(proc.returncode or 0)


if __name__ == "__main__":
    main()
