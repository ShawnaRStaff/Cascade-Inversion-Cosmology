"""M4 Option-1 analysis: void percolation transition p_v_c.

The standard percolation transition asks: at what occupation density does
the OCCUPIED set first form a spanning cluster? We measured this at
p_c ~ 0.18.

The dual question: at what density does the UNOCCUPIED (void) set
stop forming a spanning cluster? Call this p_v_c. For uncorrelated 3D
site percolation, p_v_c = 1 - p_c ~ 0.688.

If our correlated dynamics gives p_v_c ~ 0.65, that's a real natural
attractor: at p_v_c the topology of the system fundamentally changes
(voids stop spanning, the occupied cluster becomes simply
double-connected). The dynamics is genuinely special at p_v_c. If the
universe lives at the void-percolation transition, the gamma match at
p=0.65 stops being a free parameter and becomes a consequence of
substrate topology.

We compute this on the existing M4 snapshot masks. No new simulation
required.

Run:
    .venv/bin/python scripts/run_milestone4_void_percolation.py [npz]
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from void_cascade.percolation import check_spanning  # noqa: E402


def latest_xi_npz() -> Path:
    outdir = REPO_ROOT / "data" / "outputs"
    candidates = sorted(outdir.glob("manna_3d_xi_data_*.npz"))
    if not candidates:
        raise FileNotFoundError("No manna_3d_xi_data_*.npz in data/outputs/")
    return candidates[-1]


def main() -> None:
    if len(sys.argv) > 1:
        path = Path(sys.argv[1])
    else:
        path = latest_xi_npz()
    print(f"Loading {path}")
    data = np.load(path)
    L = int(data["L"])
    p_targets = sorted([float(p) for p in data["p_targets"]])
    print(f"L = {L}, p_targets = {p_targets}")
    print()

    # For each p, for each seed, compute spanning state of:
    #   (a) the occupied set (ever_toppled itself)
    #   (b) the void set (NOT ever_toppled)
    print(f"{'p':>6}  {'seed':>6}  {'occ.span':>8}  {'occ.large':>10}  "
          f"{'occ.frac':>9}  {'void.span':>9}  {'void.large':>11}  "
          f"{'void.frac':>10}")

    results: dict[float, dict] = {}
    for p in p_targets:
        pp = int(round(p * 100))
        occ_spans = []
        void_spans = []
        occ_largest_frac = []
        void_largest_frac = []
        seed_count = 0
        for k in data.files:
            if k.startswith("mask_") and k.endswith(f"_p{pp:02d}"):
                seed = k.split("_p")[0].split("mask_s")[1]
                mask = data[k]
                void = ~mask

                occ_r = check_spanning(mask)
                void_r = check_spanning(void)

                occ_total = mask.sum()
                void_total = void.sum()
                occ_largest_frac.append(
                    occ_r.largest_cluster_size / max(occ_total, 1)
                )
                void_largest_frac.append(
                    void_r.largest_cluster_size / max(void_total, 1)
                )
                occ_spans.append(int(occ_r.percolates))
                void_spans.append(int(void_r.percolates))

                print(f"{p:>6.2f}  {seed:>6}  {bool(occ_r.percolates)!s:>8}  "
                      f"{occ_r.largest_cluster_size:>10d}  "
                      f"{occ_largest_frac[-1]:>9.3f}  "
                      f"{bool(void_r.percolates)!s:>9}  "
                      f"{void_r.largest_cluster_size:>11d}  "
                      f"{void_largest_frac[-1]:>10.3f}")
                seed_count += 1

        results[p] = {
            "occ_spans_fraction": np.mean(occ_spans) if occ_spans else np.nan,
            "void_spans_fraction": np.mean(void_spans) if void_spans else np.nan,
            "occ_largest_frac": np.mean(occ_largest_frac) if occ_largest_frac else np.nan,
            "void_largest_frac": np.mean(void_largest_frac) if void_largest_frac else np.nan,
            "n_seeds": seed_count,
        }

    # Find transition: smallest p at which the void set NO LONGER spans for
    # the majority of seeds. (Looking for the p where void.spans flips
    # from True to False across the snapshot grid.)
    print()
    print("Aggregated:")
    print(f"{'p':>6}  {'occ.span':>9}  {'void.span':>10}  {'occ.dom':>9}  {'void.dom':>10}")
    p_values = sorted(results.keys())
    for p in p_values:
        r = results[p]
        print(f"{p:>6.2f}  {r['occ_spans_fraction']:>9.2f}  "
              f"{r['void_spans_fraction']:>10.2f}  "
              f"{r['occ_largest_frac']:>9.3f}  "
              f"{r['void_largest_frac']:>10.3f}")

    # Find the transition: largest p at which all seeds still have void
    # spanning, vs. smallest p at which no seed has void spanning.
    last_void_span = None
    first_no_void_span = None
    for p in p_values:
        r = results[p]
        if r["void_spans_fraction"] >= 0.99:
            last_void_span = p
        if r["void_spans_fraction"] < 0.01 and first_no_void_span is None:
            first_no_void_span = p

    print()
    if last_void_span is not None:
        print(f"Voids span for all seeds up through p = {last_void_span:.2f}")
    if first_no_void_span is not None:
        print(f"Voids fail to span starting at p = {first_no_void_span:.2f}")
    if last_void_span is not None and first_no_void_span is not None:
        pvc = 0.5 * (last_void_span + first_no_void_span)
        print(f"Void percolation transition p_v_c is between {last_void_span:.2f}"
              f" and {first_no_void_span:.2f}  (midpoint estimate {pvc:.3f})")
        print()
        print(f"Compare to gamma=1.8 crossing at p ~ 0.65")
        print(f"Compare to uncorrelated 3D site percolation p_v_c = 1 - 0.3116 = 0.688")

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    occ_spans = [results[p]["occ_spans_fraction"] for p in p_values]
    void_spans = [results[p]["void_spans_fraction"] for p in p_values]
    occ_frac = [results[p]["occ_largest_frac"] for p in p_values]
    void_frac = [results[p]["void_largest_frac"] for p in p_values]
    ax1.plot(p_values, occ_spans, "o-", color="C0", label="occupied spans (fraction of seeds)")
    ax1.plot(p_values, void_spans, "s-", color="C3", label="void spans (fraction of seeds)")
    ax1.axvline(0.65, color="0.5", linestyle="--", alpha=0.6, label="γ=1.8 crossing at p=0.65")
    ax1.axvline(0.688, color="0.3", linestyle=":", alpha=0.6, label="uncorrelated 3D dual p_v=0.688")
    ax1.set_xlabel("$p$")
    ax1.set_ylabel("fraction of seeds spanning")
    ax1.set_title("Spanning transitions: occupied (fast) and void (late)")
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8)

    ax2.plot(p_values, occ_frac, "o-", color="C0",
             label="largest occupied / total occupied")
    ax2.plot(p_values, void_frac, "s-", color="C3",
             label="largest void / total void")
    ax2.axvline(0.65, color="0.5", linestyle="--", alpha=0.6)
    ax2.set_xlabel("$p$")
    ax2.set_ylabel("largest cluster fraction within its component class")
    ax2.set_title("Cluster dominance fraction")
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=8)

    fig.tight_layout()
    outdir = REPO_ROOT / "data" / "outputs"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_void_percolation_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print()
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
