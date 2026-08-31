#!/usr/bin/env python3
"""
Simple benchmark to test the read-based approach performance.
"""

import subprocess
import time


def run_benchmark(region_size, iterations=3):
    """Run benchmark with different region sizes."""

    base_cmd = [
        "uv",
        "run",
        "countmut",
        "-i",
        "/home/yec/Desktop/S818.genes.bam",
        "-r",
        "/home/yec/Desktop/genes.fa",
        "--strand",
        "forward",
        "--region",
        f"rRNA-Hsa-nucleus_locus:1-{region_size}",
    ]

    print(f"🧬 Benchmarking region size: {region_size} bp")
    print("=" * 50)

    times = []
    for i in range(iterations):
        print(f"Run {i + 1}/{iterations}...")
        start_time = time.time()

        result = subprocess.run(base_cmd, capture_output=True, text=True)

        end_time = time.time()
        elapsed = end_time - start_time
        times.append(elapsed)

        if result.returncode != 0:
            print(f"❌ Error in run {i + 1}: {result.stderr}")
            return None

        print(f"   Time: {elapsed:.2f}s")

    avg_time = sum(times) / len(times)
    print(f"\n📊 Results for {region_size} bp:")
    print(f"   Average time: {avg_time:.2f}s")
    print(f"   Min time: {min(times):.2f}s")
    print(f"   Max time: {max(times):.2f}s")
    print(f"   Processing rate: {region_size / avg_time:.0f} bp/s")

    return avg_time


if __name__ == "__main__":
    print("🚀 CountMut Read-Based Processing Benchmark")
    print("=" * 60)

    # Test different region sizes
    region_sizes = [500, 1000, 2000, 5000]
    results = {}

    for size in region_sizes:
        results[size] = run_benchmark(size, iterations=2)
        print()

    print("📈 Performance Summary:")
    print("=" * 30)
    for size, time_taken in results.items():
        if time_taken:
            rate = size / time_taken
            print(f"{size:5d} bp: {time_taken:6.2f}s ({rate:6.0f} bp/s)")

    print("\n✅ Benchmark completed!")
