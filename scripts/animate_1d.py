"""Render a GIF of 1D Oslo dynamics.

Run with the repo venv:
    .venv/bin/python scripts/animate_1d.py

Defaults are tuned to produce a short, watchable animation showing the pile
filling from empty and a handful of avalanches once steady state is reached.
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from void_cascade.visualize_1d import record_animation  # noqa: E402


def main() -> None:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out = REPO_ROOT / "data" / "outputs" / f"oslo_1d_dynamics_{stamp}.gif"
    print(f"Rendering Oslo dynamics animation to {out}")
    path = record_animation(L=24, max_frames=600, seed=42, out_path=out, fps=20)
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
