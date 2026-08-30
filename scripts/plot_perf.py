#!/usr/bin/env python3
"""Render the countmut performance chart (see docs/perf-scaling.png)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# median-of-3 seconds on the bimodal benchmark
# (232k reads / 23.2M read-bases / 3002 contigs / 2 ~3k-deep rRNA hotspots)
threads = [1, 2, 4, 8, 16]
rw_mut = [2.14, 1.19, 1.01, 0.95, 0.92]
pl_base = [1.96, 1.10, 0.83, 0.70, 0.69]

# -e overhead at threads=8 (read-walk mutation)
expr_labels = ["baseline\n(no -e)", "read-const\nmapq >= 20", "per-base\nbq >= 20 and\ndist5 >= 2"]
expr_times = [0.95, 0.96, 1.23]

plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.2))

ax1.plot(threads, rw_mut, "o-", color="#d62728", label="read-walk · mutation", lw=2)
ax1.plot(threads, pl_base, "s-", color="#1f77b4", label="pileup · base", lw=2)
ax1.set_xlabel("threads")
ax1.set_ylabel("wall-clock (s)")
ax1.set_title("Thread scaling (median of 3)")
ax1.set_ylim(0, 2.4)
ax1.legend(frameon=False)

ax2.bar(expr_labels, expr_times, color=["#7f7f7f", "#2ca02c", "#d62728"], width=0.6)
ax2.set_ylabel("wall-clock (s) @ 8 threads")
ax2.set_title("Filter overhead (read-walk · mutation)")
ax2.set_ylim(0, 1.5)
for x, y in zip(range(3), expr_times):
    ax2.text(x, y + 0.03, f"{y:.2f}", ha="center", fontsize=10)

fig.suptitle("countmut — bimodal benchmark (232k reads, deep rRNA hotspots)", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig("docs/perf-scaling.png", dpi=150)
print("saved docs/perf-scaling.png")
