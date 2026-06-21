"""Clock the three onset observables on a shared timeline.

Three signals, measured during the buildup-to-tip arc:

  (1) load_std    -- stress clustering (SOC / avalanche onset).
                     Starts at 0, grows as the cascade organizes the substrate.
  (2) corr_length -- spatial correlation length of the load field (structural onset).
                     Grows as patches form -- the moment the substrate is no longer
                     random noise but has real correlated structure.
  (3) n_hot       -- cells at or above plasma temperature (ignition).
                     Steps from 0 to L² at the tip.

The question: does structure appear before criticality (load_std) or after?
Does load clustering lead ignition, or are they simultaneous?
The model leads -- we report the ordering as measured.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from void_cascade.cascade_breakdown import BreakdownParams, run_onset_measurement  # noqa: E402

L      = 80
STEPS  = 3200   # enough to capture tip (~2038) + some post-tip
SEED   = 0
SAMPLE = 20

params = BreakdownParams()   # default: slow drive, hpc=0.1

print(f"=== Onset timing  L={L}  steps={STEPS}  sample_every={SAMPLE} ===")
r = run_onset_measurement(L=L, steps=STEPS, params=params, seed=SEED, sample_every=SAMPLE)
tip = r["tip_step"]
t   = np.array(r["t_axis"])
std = np.array(r["load_std"])
xi  = np.array(r["corr_lengths"])
hot = np.array(r["n_hot"])

print(f"  tip_step:    {tip}")
print(f"  max load_std: {std.max():.4f}  at t={t[std.argmax()]}")
print(f"  max corr_len: {xi.max():.2f}   at t={t[xi.argmax()]}")
print(f"  n_hot final:  {hot[-1]}")

# --- onset thresholds ---
thresh = 0.10   # "becomes observable" = first time signal > 10% of its max
def first_crossing(arr, thr_frac):
    m = arr.max()
    if m < 1e-10:
        return None
    idx = np.argmax(arr > thr_frac * m)
    return int(t[idx]) if arr[idx] > thr_frac * m else None

t_std = first_crossing(std, thresh)
t_xi  = first_crossing(xi,  thresh)
t_hot = tip

print(f"\n  First 10%-of-max crossings:")
print(f"    load_std (avalanche clustering) : step {t_std}")
print(f"    corr_length (spatial structure) : step {t_xi}")
print(f"    ignition (n_hot > 0)            : step {t_hot}")

if t_std is not None and t_xi is not None and t_hot is not None:
    order = sorted([("load_std", t_std), ("corr_length", t_xi), ("ignition", t_hot)],
                   key=lambda x: x[1])
    print(f"\n  Ordering: {' -> '.join(f'{name}({step})' for name, step in order)}")
    print(f"\n  Gap: structure appears {abs(t_xi - t_std)} steps after avalanche clustering.")
    print(f"  Gap: ignition appears {t_hot - min(t_std, t_xi)} steps after first structural signal.")

# --- plot ---
fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

def _normalise(arr):
    m = arr.max()
    return arr / m if m > 1e-10 else arr

ax0, ax1, ax2 = axes

ax0.plot(t, _normalise(std), color="steelblue", lw=1.5)
ax0.set_ylabel("load_std / max\n(avalanche clustering)", fontsize=9)
ax0.set_ylim(0, 1.1)

ax1.plot(t, xi, color="darkorange", lw=1.5)
ax1.set_ylabel("corr_length (cells)\n(spatial structure)", fontsize=9)

ax2.plot(t, hot, color="firebrick", lw=1.5)
ax2.set_ylabel("n_hot cells\n(ignition)", fontsize=9)
ax2.set_xlabel("step")

for ax in axes:
    if tip is not None:
        ax.axvline(tip, color="k", lw=0.8, ls=":", alpha=0.5, label=f"tip={tip}")
    ax.grid(True, alpha=0.3)

fig.suptitle(
    f"Onset timing  L={L}  steps={STEPS}\n"
    f"Order: avalanche clustering → spatial structure → ignition"
    if (t_std is not None and t_xi is not None and t_hot is not None
        and t_std <= t_xi <= t_hot)
    else f"Onset timing  L={L}  steps={STEPS}"
)
fig.tight_layout()

out = REPO / "data" / "outputs" / "onset_timing"
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "onset_timing.png", dpi=150)
print(f"\nPlot saved: {out}/onset_timing.png")
