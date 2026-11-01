#!/usr/bin/env python3
"""Plot the available follow-up probe diffusion and charge-sign estimates."""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
STATES = ["B3", "U4", "H5", "U6", "H7"]
CAMPAIGNS = ["long_N4000", "N8000"]


def main():
    summary = json.loads((ROOT / "results" / "summary.json").read_text())
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    for offset, campaign in enumerate(CAMPAIGNS):
        systems = summary["systems"][campaign]
        for j, code in enumerate(STATES):
            probe = systems[code]["D_probe"]
            if probe["mean"] is None:
                continue
            x = j + (offset - 0.5) * 0.24
            lo, hi = (probe["CI"] or [probe["mean"]] * 2)
            axes[0].errorbar(x, probe["mean"],
                             yerr=[[probe["mean"] - lo], [hi - probe["mean"]]],
                             fmt="o", capsize=3, color=f"C{offset}")
            for sign, marker, shift in [("D_plus2", "^", -0.04), ("D_minus2", "v", 0.04)]:
                item = systems[code][sign]
                if item["mean"] is not None:
                    low, high = (item["CI"] or [item["mean"]] * 2)
                    axes[1].errorbar(x + shift, item["mean"],
                                     yerr=[[item["mean"] - low], [high - item["mean"]]],
                                     fmt=marker, capsize=2, color=f"C{offset}")
    for ax in axes:
        ax.set_xticks(np.arange(len(STATES)), STATES)
        ax.set_ylabel("Diffusion coefficient (LJ units)")
        ax.grid(alpha=0.2)
    axes[0].set_title("D_probe: circles, 95% replica CI")
    axes[1].set_title("D(+2): triangles up; D(−2): triangles down")
    fig.text(0.5, 0.01, "Blue: long N=4000; orange: N=8000", ha="center")
    target = ROOT / "results" / "followup_transport.pdf"
    fig.savefig(target)
    plt.close(fig)
    print(target)


if __name__ == "__main__":
    main()
