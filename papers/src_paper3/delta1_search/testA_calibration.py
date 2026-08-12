import numpy as np

TOL = 6.15e-5

vals_sorted = np.load("testA_vals_sorted.npy")

rng = np.random.RandomState(42)
targets = rng.uniform(11.5, 14.0, size=200)

counts = np.empty(200, dtype=np.int64)
for i, t in enumerate(targets):
    lo = t * (1 - TOL)
    hi = t * (1 + TOL)
    l = np.searchsorted(vals_sorted, lo, side="left")
    r = np.searchsorted(vals_sorted, hi, side="right")
    counts[i] = r - l

print("mean hits:", counts.mean())
print("min hits:", counts.min())
print("median hits:", np.median(counts))
print("max hits:", counts.max())

with open("testA_calibration.txt", "w") as f:
    f.write(f"200 targets uniform in [11.5,14.0], seed=42, tol={TOL}\n")
    f.write(f"mean={counts.mean():.3f} min={counts.min()} median={np.median(counts):.1f} max={counts.max()}\n")
    f.write("Observed hits for D1: 40\n")
    f.write(f"Full distribution (sorted): {sorted(counts.tolist())}\n")

print("Observed D1 hits: 40  vs expected mean:", counts.mean())
