#!/usr/bin/env python3
"""Real-data benchmark figures for the ClinVar 2-star+ evaluation.

Reads outputs from `03_evaluate_concordance.py` (under
`data/benchmark/output_v5/` by default) and emits 5 PDF/PNG panels
under `<output_dir>/figures/`:

  fig_concordance_matrix       row-normalised heatmap of truth × predicted
  fig_recall_by_class          per-class same-direction recall (v4)
  fig_v1_vs_v5_recall          paired bars showing the lift from loading
                               PhyloP+SpliceAI+ClinGen GDV (the v1
                               baseline is hard-coded from the prior run
                               whose results are recorded in METHODS.md;
                               it predates the SA-wiring fixes)
  fig_criterion_fires          top criteria by total fire count (split
                               by truth direction)
  fig_bp7_pvs1_delta           the two single-criterion deltas that
                               drove the recall lift

Usage:
  generate_figures.py                     # uses ../../data/benchmark/output_v5
  generate_figures.py <out_dir>
"""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

plt.rcParams.update(
    {
        "font.size": 13,
        "axes.titlesize": 17,
        "axes.labelsize": 15,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 12,
        "legend.title_fontsize": 13,
    }
)

C = {
    "P": "#dc2626",
    "LP": "#f97316",
    "VUS": "#6b7280",
    "LB": "#3b82f6",
    "B": "#10b981",
    "v1": "#94a3b8",
    "v4": "#6c7aee",
    "delta_up": "#10b981",
    "delta_down": "#ef4444",
}
CLASSES = ["Pathogenic", "Likely_pathogenic", "VUS", "Likely_benign", "Benign"]
CLASS_SHORT = {
    "Pathogenic": "P",
    "Likely_pathogenic": "LP",
    "VUS": "VUS",
    "Likely_benign": "LB",
    "Benign": "B",
}
CLASS_COLOR = {
    "Pathogenic": C["P"],
    "Likely_pathogenic": C["LP"],
    "VUS": C["VUS"],
    "Likely_benign": C["LB"],
    "Benign": C["B"],
}

# ──────────────────────────────────────────────────────────────────────
# v1 baseline — measured on the same 673,660-variant ClinVar 2-star+ set
# *before* PhyloP / SpliceAI / ClinGen GDV were loaded and before the
# SpliceAI-camelCase / PhyloP-routing wiring fixes. Captured here so the
# v1↔v4 comparison panel doesn't need a re-run of the prior pipeline.
# Source: METHODS.md "Real-Data Concordance" section as written for v1.
# ──────────────────────────────────────────────────────────────────────
V1_RECALL = {
    "Pathogenic": 15.7,
    "Likely_pathogenic": 20.9,
    "VUS": 96.6,
    "Likely_benign": 3.2,
    "Benign": 33.2,
}
V1_HEADLINE = {
    "exact_match": 52.7,
    "same_direction": 54.7,
    "opposite_direction": 0.005,
    "no_call": 0.0,
}
# v1 BP7 fired 0 times; PVS1 fired 5,233 (Pathogenic) + 403 (LP) = 5,636.
# Used in the BP7+PVS1 delta panel.
V1_FIRES = {
    "BP7": 0,
    "PVS1": 5_233 + 403,
    "PVS1_Supporting": 47 + 2,
    "PS1": 9_068 + 3_240,
    "BS2": 333 + 34 + 12_866 + 75_504,
    "BA1": 1 + 874 + 41_183,
}


def read_matrix(out_dir: Path):
    """Return a dict[(truth, predicted)] -> count from concordance_matrix.csv."""
    rows = list(csv.reader((out_dir / "concordance_matrix.csv").open()))
    header = rows[0]  # ["truth"] + CLASSES + ["NoCall"]
    matrix = {}
    for row in rows[1:]:
        truth = row[0]
        for col, val in zip(header[1:], row[1:]):
            matrix[(truth, col)] = int(val)
    return matrix


def parse_summary(out_dir: Path) -> dict[str, float]:
    """Pull headline metrics out of concordance_summary.txt."""
    text = (out_dir / "concordance_summary.txt").read_text()
    out = {}
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("Exact-match rate:"):
            out["exact_match"] = float(s.split(":")[1].strip().rstrip("%"))
        elif s.startswith("Same-direction rate:"):
            out["same_direction"] = float(s.split(":")[1].strip().rstrip("%"))
        elif s.startswith("Opposite-direction rate:"):
            out["opposite_direction"] = float(s.split(":")[1].strip().rstrip("%"))
        elif s.startswith("NoCall rate:"):
            out["no_call"] = float(s.split(":")[1].strip().rstrip("%"))
    return out


def read_criterion_fires(out_dir: Path) -> dict[str, dict[str, int]]:
    """Return {criterion: {Pathogenic: n, Likely_pathogenic: n, ...}}."""
    fires = {}
    with (out_dir / "criterion_firing_rates.csv").open() as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            code = row["criterion"]
            fires[code] = {tcl: int(row.get(f"{tcl}_fired", 0)) for tcl in CLASSES}
    return fires


def fig_concordance_matrix(matrix, out_dir: Path, fig_dir: Path):
    fig, ax = plt.subplots(figsize=(11, 8))
    cols = CLASSES + ["NoCall"]
    mat = np.zeros((len(CLASSES), len(cols)))
    for i, t in enumerate(CLASSES):
        row_total = sum(matrix.get((t, c), 0) for c in cols)
        for j, c in enumerate(cols):
            mat[i, j] = 100 * matrix.get((t, c), 0) / row_total if row_total else 0

    im = ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(
        [CLASS_SHORT.get(c, c) for c in cols], fontsize=14, fontweight="bold"
    )
    ax.set_yticks(range(len(CLASSES)))
    ax.set_yticklabels([CLASS_SHORT[c] for c in CLASSES], fontsize=14, fontweight="bold")
    ax.set_xlabel("fastVEP predicted", fontweight="bold")
    ax.set_ylabel("ClinVar 2-star+ truth", fontweight="bold")
    ax.set_title(
        "ACMG concordance matrix (row-normalised %, full SA stack)",
        fontweight="bold",
    )

    for i, t in enumerate(CLASSES):
        row_total = sum(matrix.get((t, c), 0) for c in cols)
        for j, c in enumerate(cols):
            cnt = matrix.get((t, c), 0)
            pct = mat[i, j]
            color = "white" if pct > 50 else "black"
            label = f"{cnt:,}\n({pct:.0f}%)"
            ax.text(
                j,
                i,
                label,
                ha="center",
                va="center",
                fontsize=11,
                color=color,
                fontweight="bold" if cols[j] == t else "normal",
            )
        # diagonal frame
        if i < len(cols):
            ax.add_patch(
                plt.Rectangle(
                    (i - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=C["delta_up"], lw=3
                )
            )

    cbar = plt.colorbar(im, ax=ax, shrink=0.85)
    cbar.set_label("Row %")
    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(fig_dir / f"fig_concordance_matrix.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  fig_concordance_matrix")


def fig_recall_by_class(matrix, out_dir: Path, fig_dir: Path):
    """Stacked bar: per-class outcome share (Same / Opp / NoCall / VUS-call)."""
    fig, ax = plt.subplots(figsize=(11, 6))
    n_per = {t: sum(matrix.get((t, c), 0) for c in CLASSES + ["NoCall"]) for t in CLASSES}

    def share(truth, mask):
        n = n_per[truth] or 1
        return 100 * sum(matrix.get((truth, c), 0) for c in mask) / n

    same_dir_mask = {
        "Pathogenic": ["Pathogenic", "Likely_pathogenic"],
        "Likely_pathogenic": ["Pathogenic", "Likely_pathogenic"],
        "VUS": ["VUS"],
        "Likely_benign": ["Likely_benign", "Benign"],
        "Benign": ["Likely_benign", "Benign"],
    }
    opp_mask = {
        "Pathogenic": ["Benign", "Likely_benign"],
        "Likely_pathogenic": ["Benign", "Likely_benign"],
        "VUS": [],
        "Likely_benign": ["Pathogenic", "Likely_pathogenic"],
        "Benign": ["Pathogenic", "Likely_pathogenic"],
    }

    same = [share(t, same_dir_mask[t]) for t in CLASSES]
    nocall = [share(t, ["NoCall"]) for t in CLASSES]
    opp = [share(t, opp_mask[t]) for t in CLASSES]
    other = [100 - s - n - o for s, n, o in zip(same, nocall, opp)]

    x = np.arange(len(CLASSES))
    ax.bar(x, same, color=C["delta_up"], label="Same direction", alpha=0.9)
    ax.bar(x, other, bottom=same, color=C["VUS"], alpha=0.7, label="VUS / off-direction non-opposite")
    ax.bar(
        x,
        nocall,
        bottom=[s + o for s, o in zip(same, other)],
        color="#fbbf24",
        alpha=0.85,
        label="NoCall (no ACMG returned)",
    )
    ax.bar(
        x,
        opp,
        bottom=[s + o + n for s, o, n in zip(same, other, nocall)],
        color=C["delta_down"],
        alpha=0.9,
        label="Opposite direction",
    )

    ax.set_xticks(x)
    ax.set_xticklabels([CLASS_SHORT[c] for c in CLASSES], fontweight="bold", fontsize=14)
    ax.set_xlabel("ClinVar 2-star+ truth class", fontweight="bold")
    ax.set_ylabel("% of class")
    ax.set_ylim(0, 105)
    ax.set_title("Per-class outcome breakdown", fontweight="bold")
    ax.legend(loc="upper right")
    ax.grid(axis="y", alpha=0.15)

    for i, (s, n) in enumerate(zip(same, [n_per[t] for t in CLASSES])):
        ax.text(
            i, s + 1.5, f"{s:.0f}%\nn={n:,}", ha="center", fontsize=11, fontweight="bold"
        )

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(fig_dir / f"fig_recall_by_class.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  fig_recall_by_class")


def fig_v1_vs_v5_recall(matrix, fig_dir: Path):
    """Paired bars of per-class same-direction recall, v1 baseline vs v4."""
    fig, ax = plt.subplots(figsize=(11, 6))

    def v4_recall(truth):
        n = sum(matrix.get((truth, c), 0) for c in CLASSES + ["NoCall"]) or 1
        if truth in ("Pathogenic", "Likely_pathogenic"):
            same = matrix.get((truth, "Pathogenic"), 0) + matrix.get(
                (truth, "Likely_pathogenic"), 0
            )
        elif truth == "VUS":
            same = matrix.get((truth, "VUS"), 0)
        else:
            same = matrix.get((truth, "Benign"), 0) + matrix.get(
                (truth, "Likely_benign"), 0
            )
        return 100 * same / n

    v1 = [V1_RECALL[c] for c in CLASSES]
    v4 = [v4_recall(c) for c in CLASSES]

    x = np.arange(len(CLASSES))
    width = 0.38
    bars_v1 = ax.bar(
        x - width / 2,
        v1,
        width,
        color=C["v1"],
        alpha=0.9,
        label="v1: REVEL + gnomAD + ClinVar only",
    )
    bars_v4 = ax.bar(
        x + width / 2,
        v4,
        width,
        color=C["v4"],
        alpha=0.95,
        label="v5: + PhyloP + SpliceAI + ClinGen GDV + indel allele fix (current)",
    )
    ax.set_xticks(x)
    ax.set_xticklabels([CLASS_SHORT[c] for c in CLASSES], fontweight="bold", fontsize=14)
    ax.set_ylabel("Same-direction recall (%)")
    ax.set_ylim(0, 105)
    ax.set_title(
        "Recall lift from loading PhyloP + SpliceAI + ClinGen Gene-Disease Validity",
        fontweight="bold",
    )
    ax.grid(axis="y", alpha=0.15)
    ax.legend(loc="upper center")

    for x0, b1, b4 in zip(x, v1, v4):
        ax.text(x0 - width / 2, b1 + 1, f"{b1:.0f}%", ha="center", fontsize=10)
        ax.text(x0 + width / 2, b4 + 1, f"{b4:.0f}%", ha="center", fontsize=10, fontweight="bold")
        delta = b4 - b1
        color = C["delta_up"] if delta > 0 else C["delta_down"]
        ax.annotate(
            f"{delta:+.0f} pp",
            xy=(x0, max(b1, b4) + 7),
            ha="center",
            fontsize=11,
            color=color,
            fontweight="bold",
        )

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(fig_dir / f"fig_v1_vs_v5_recall.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  fig_v1_vs_v5_recall")


def fig_criterion_fires(fires: dict[str, dict[str, int]], fig_dir: Path):
    """Stacked horizontal bars: top criteria, fires colored by truth class."""
    totals = sorted(
        fires.items(), key=lambda kv: -sum(kv[1].values())
    )[:18]
    codes = [c for c, _ in totals]

    fig, ax = plt.subplots(figsize=(13, 9))
    y = np.arange(len(codes))
    left = np.zeros(len(codes))
    for tcl in CLASSES:
        vals = np.array([fires[c].get(tcl, 0) for c in codes])
        ax.barh(y, vals, left=left, color=CLASS_COLOR[tcl], alpha=0.9, label=CLASS_SHORT[tcl])
        left += vals

    ax.set_yticks(y)
    ax.set_yticklabels(codes, fontfamily="monospace", fontsize=12)
    ax.invert_yaxis()
    ax.set_xlabel("Times the criterion was met (across all 627k classified variants)")
    ax.set_title(
        "Criterion fire counts by truth class\n(stacked: each colour = ClinVar truth label)",
        fontweight="bold",
    )
    ax.legend(loc="lower right", title="Truth class")
    ax.grid(axis="x", alpha=0.15)

    for yi, code in enumerate(codes):
        total = sum(fires[code].values())
        ax.text(total + total * 0.005, yi, f"{total:,}", va="center", fontsize=10, color="#374151")

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(fig_dir / f"fig_criterion_fires.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  fig_criterion_fires")


def fig_bp7_pvs1_delta(fires: dict[str, dict[str, int]], fig_dir: Path):
    """Lollipop comparing v1 vs v4 fire counts for the criteria most
    affected by the SA-source additions (BP7 from PhyloP+SpliceAI;
    PVS1 from ClinGen GDV)."""
    keys = ["BP7", "PVS1", "PVS1_Supporting", "PS1", "BS2", "BA1"]
    v4_total = {k: sum(fires.get(k, {}).get(t, 0) for t in CLASSES) for k in keys}
    v1_total = {k: V1_FIRES.get(k, 0) for k in keys}

    fig, ax = plt.subplots(figsize=(11, 6))
    y = np.arange(len(keys))
    for yi, k in enumerate(keys):
        v1n = v1_total[k]
        v4n = v4_total[k]
        ax.plot([v1n, v4n], [yi, yi], color="#cbd5e1", lw=3, zorder=1)
        ax.scatter([v1n], [yi], color=C["v1"], s=130, zorder=2, label="v1" if yi == 0 else None)
        ax.scatter([v4n], [yi], color=C["v4"], s=160, zorder=3, label="v4" if yi == 0 else None)
        delta = v4n - v1n
        rate = (v4n / max(v1n, 1)) if v1n else float("inf")
        annot = f"  Δ {delta:+,}" + (f" ({rate:.1f}×)" if v1n else "  (new)")
        ax.text(max(v1n, v4n), yi, annot, va="center", fontsize=11, color=C["delta_up"], fontweight="bold")

    ax.set_yticks(y)
    ax.set_yticklabels(keys, fontfamily="monospace", fontsize=13)
    ax.invert_yaxis()
    ax.set_xscale("symlog", linthresh=10)
    ax.set_xlabel("Times the criterion fired (log scale)")
    ax.set_title(
        "Single-criterion lift: v1 (REVEL + gnomAD only) → v4 (full SA stack)",
        fontweight="bold",
    )
    ax.grid(axis="x", alpha=0.15)
    ax.legend(loc="lower right")

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(fig_dir / f"fig_bp7_pvs1_delta.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  fig_bp7_pvs1_delta")


def fig_headline_v1_vs_v5(out_dir: Path, fig_dir: Path):
    headline = parse_summary(out_dir)
    fig, ax = plt.subplots(figsize=(10, 5.5))
    metrics = [
        ("Same\ndirection", "same_direction"),
        ("Exact\nmatch", "exact_match"),
        ("Opposite\ndirection", "opposite_direction"),
        ("NoCall", "no_call"),
    ]
    x = np.arange(len(metrics))
    width = 0.38
    v1 = [V1_HEADLINE[k] for _, k in metrics]
    v4 = [headline.get(k, 0.0) for _, k in metrics]
    ax.bar(x - width / 2, v1, width, color=C["v1"], alpha=0.9, label="v1")
    ax.bar(x + width / 2, v4, width, color=C["v4"], alpha=0.95, label="v5")

    ax.set_xticks(x)
    ax.set_xticklabels([lab for lab, _ in metrics], fontsize=13, fontweight="bold")
    ax.set_ylabel("%")
    ax.set_title("Headline concordance metrics: v1 vs v5", fontweight="bold")
    ax.legend()
    ax.grid(axis="y", alpha=0.15)

    for x0, a, b in zip(x, v1, v4):
        ax.text(x0 - width / 2, a + 1, f"{a:.1f}", ha="center", fontsize=11)
        ax.text(x0 + width / 2, b + 1, f"{b:.1f}", ha="center", fontsize=11, fontweight="bold")

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(fig_dir / f"fig_headline_v1_vs_v5.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  fig_headline_v1_vs_v5")


def main():
    if len(sys.argv) > 1:
        out_dir = Path(sys.argv[1])
    else:
        out_dir = Path(__file__).resolve().parents[3] / "data/benchmark/output_v5"
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading {out_dir}")
    matrix = read_matrix(out_dir)
    fires = read_criterion_fires(out_dir)

    print("Generating figures...")
    fig_concordance_matrix(matrix, out_dir, fig_dir)
    fig_recall_by_class(matrix, out_dir, fig_dir)
    fig_v1_vs_v5_recall(matrix, fig_dir)
    fig_headline_v1_vs_v5(out_dir, fig_dir)
    fig_criterion_fires(fires, fig_dir)
    fig_bp7_pvs1_delta(fires, fig_dir)

    print(f"\nDone. {len(list(fig_dir.glob('*.png')))} PNG / {len(list(fig_dir.glob('*.pdf')))} PDF in {fig_dir}")


if __name__ == "__main__":
    main()
