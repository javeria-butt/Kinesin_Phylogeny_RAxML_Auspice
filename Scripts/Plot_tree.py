"""
plot_tree.py
────────────
Visualise a Newick phylogenetic tree with family / directionality colouring.

Usage:
    python scripts/plot_tree.py \
        --tree    results/kinesin_tree.nwk \
        --mode    kinesin \
        --outdir  figures

    python scripts/plot_tree.py \
        --tree    results/kinesin14_tree.nwk \
        --mode    kinesin14 \
        --outdir  figures
"""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from Bio import Phylo

# ── Colour palettes ────────────────────────────────────────────────
FAMILY_COLORS = {
    "Kinesin-1":  "#2196F3",
    "Kinesin-2":  "#4CAF50",
    "Kinesin-3":  "#FF9800",
    "Kinesin-4":  "#9C27B0",
    "Kinesin-5":  "#00BCD4",
    "Kinesin-6":  "#FF5722",
    "Kinesin-7":  "#607D8B",
    "Kinesin-8":  "#8BC34A",
    "Kinesin-13": "#FFC107",
    "Kinesin-14": "#F44336",
}

DIRECTION_COLORS = {
    "minus":   "#D32F2F",
    "unknown": "#78909C",
}

SPECIES_COLORS = {
    "Hs": "#1565C0",
    "Dm": "#2E7D32",
    "Sc": "#F57F17",
    "At": "#6A1B9A",
    "Ce": "#00695C",
    "Cg": "#BF360C",
}


# ── Label parsers ──────────────────────────────────────────────────
def label_to_family(label: str) -> str:
    for part in (label or "").split("_"):
        m = re.match(r"Kinesin(\d+)", part)
        if m:
            return f"Kinesin-{m.group(1)}"
    return "Unknown"


def label_to_direction(label: str) -> str:
    parts = (label or "").split("_")
    return parts[2] if len(parts) > 2 else "unknown"


def label_to_species(label: str) -> str:
    return (label or "").split("_")[0]


# ── Rectangular tree ───────────────────────────────────────────────
def draw_rectangular(tree, color_fn, legend_handles, title, output):
    n_tips = tree.count_terminals()
    fig, ax = plt.subplots(figsize=(14, max(8, n_tips * 0.38)))

    Phylo.draw(
        tree, axes=ax, do_show=False,
        label_func=lambda c: c.name if c.is_terminal() else "",
        branch_labels=lambda c: (
            f"{c.confidence:.0f}" if c.confidence and not c.is_terminal() else ""
        ),
        label_colors=color_fn,
    )

    ax.set_title(title, fontsize=12, fontweight="bold", pad=10)
    ax.set_xlabel("Branch Length (substitutions per site)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(handles=legend_handles, fontsize=8, loc="lower right",
              framealpha=0.9, ncol=2)

    plt.tight_layout()
    plt.savefig(output, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {output}")


# ── Circular tree ──────────────────────────────────────────────────
def draw_circular(tree, color_map_fn, color_palette, legend_handles, title, output):
    terminals = tree.get_terminals()
    n = len(terminals)
    angles = {t.name: 2 * np.pi * i / n for i, t in enumerate(terminals)}

    fig = plt.figure(figsize=(12, 12))
    ax  = fig.add_subplot(111, projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.axis("off")
    ax.set_title(title, fontsize=12, fontweight="bold", pad=20)

    for t in terminals:
        ang = angles[t.name]
        col = color_map_fn(t.name)
        ax.text(ang, 1.06, t.name, ha="center", va="center",
                fontsize=6, color=col,
                rotation=np.degrees(ang), rotation_mode="anchor")
        ax.plot([ang, ang], [0.96, 1.01], color=col, lw=2)

    # Arcs per group
    groups = {}
    for t in terminals:
        key = color_map_fn(t.name)
        groups.setdefault(key, []).append(angles[t.name])

    for col, ang_list in groups.items():
        if len(ang_list) < 2:
            continue
        arc = np.linspace(min(ang_list), max(ang_list), 60)
        ax.plot(arc, [0.93] * 60, color=col, lw=3, alpha=0.75,
                solid_capstyle="round")

    ax.legend(handles=legend_handles, fontsize=7, loc="lower center",
              bbox_to_anchor=(0.5, -0.04), ncol=5, frameon=False)

    plt.savefig(output, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {output}")


# ── Main ───────────────────────────────────────────────────────────
def plot_kinesin(tree_file: str, outdir: str):
    tree = Phylo.read(tree_file, "newick")
    out  = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    color_fn  = lambda c: FAMILY_COLORS.get(label_to_family(c.name or ""), "#888")
    handles   = [mpatches.Patch(color=c, label=f) for f, c in FAMILY_COLORS.items()]

    draw_rectangular(tree, color_fn, handles,
                     "Kinesin Motor Protein Superfamily — Phylogenetic Tree",
                     out / "kinesin_phylotree.png")

    draw_circular(tree,
                  lambda lbl: FAMILY_COLORS.get(label_to_family(lbl), "#888"),
                  FAMILY_COLORS, handles,
                  "Kinesin Superfamily — Circular Phylogeny",
                  out / "kinesin_circular_tree.png")


def plot_kinesin14(tree_file: str, outdir: str):
    tree = Phylo.read(tree_file, "newick")
    out  = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    color_fn = lambda c: DIRECTION_COLORS.get(label_to_direction(c.name or ""), "#888")
    handles  = [mpatches.Patch(color=c, label=f"Direction: {d}")
                for d, c in DIRECTION_COLORS.items()]

    draw_rectangular(tree, color_fn, handles,
                     "Kinesin-14 Family — Phylogenetic Tree\n"
                     "(Red = confirmed minus-end; numbers = bootstrap)",
                     out / "kinesin14_phylotree.png")

    sp_color_fn = lambda lbl: SPECIES_COLORS.get(label_to_species(lbl), "#888")
    sp_handles  = [mpatches.Patch(color=c, label=sp)
                   for sp, c in SPECIES_COLORS.items()]

    draw_circular(tree, sp_color_fn, SPECIES_COLORS, sp_handles,
                  "Kinesin-14 Family — Circular Phylogeny (by Species)",
                  out / "kinesin14_circular_tree.png")


def main():
    parser = argparse.ArgumentParser(description="Plot phylogenetic tree figures")
    parser.add_argument("--tree",   required=True, help="Input Newick tree file")
    parser.add_argument("--mode",   choices=["kinesin", "kinesin14"], required=True)
    parser.add_argument("--outdir", default="figures")
    args = parser.parse_args()

    print(f"Plotting {args.mode} tree → {args.outdir}/")
    if args.mode == "kinesin":
        plot_kinesin(args.tree, args.outdir)
    else:
        plot_kinesin14(args.tree, args.outdir)
    print("Done.")


if __name__ == "__main__":
    main()
