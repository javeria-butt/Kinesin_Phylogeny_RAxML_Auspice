"""
auspice_json_builder.py
───────────────────────
Convert a Newick tree into an Auspice v2 JSON for interactive viewing at
https://auspice.us  (drag-and-drop the output JSON file).

Usage:
    python scripts/auspice_json_builder.py \
        --tree   results/kinesin_tree.nwk \
        --mode   kinesin \
        --output auspice/kinesin_auspice.json

    python scripts/auspice_json_builder.py \
        --tree   results/kinesin14_tree.nwk \
        --mode   kinesin14 \
        --output auspice/kinesin14_auspice.json
"""

import argparse
import json
import re
from pathlib import Path
from Bio import Phylo


# ── Label parsers ──────────────────────────────────────────────────
def label_to_family(label: str) -> str:
    for part in (label or "").split("_"):
        m = re.match(r"Kinesin(\d+)", part)
        if m:
            return f"Kinesin-{m.group(1)}"
    return "Unknown"


def label_to_direction(label: str) -> str:
    parts = (label or "").split("_")
    raw   = parts[2] if len(parts) > 2 else "unknown"
    return "Minus-end" if raw == "minus" else "Unknown"


SPECIES_FULL = {
    "Hs": "Homo sapiens",
    "Dm": "Drosophila melanogaster",
    "Sc": "Saccharomyces cerevisiae",
    "At": "Arabidopsis thaliana",
    "Ce": "Caenorhabditis elegans",
    "Cg": "Cricetulus griseus",
}


def label_to_species(label: str) -> str:
    sp = (label or "").split("_")[0]
    return SPECIES_FULL.get(sp, sp)


def label_to_gene(label: str) -> str:
    parts = (label or "").split("_")
    return parts[1] if len(parts) > 1 else "?"


# ── Clade → dict ───────────────────────────────────────────────────
def clade_to_dict_kinesin(clade, depth=0):
    label = clade.name or f"node_{depth}_{id(clade)}"
    node  = {
        "name": label,
        "node_attrs": {
            "div":           clade.branch_length or 0,
            "family":        {"value": label_to_family(label)},
            "species":       {"value": label_to_species(label)},
            "gene":          {"value": label_to_gene(label)},
            "directionality": {
                "value": "Minus-end" if label_to_family(label) == "Kinesin-14" else "Plus-end"
            },
        },
    }
    if clade.confidence is not None:
        node["node_attrs"]["bootstrap"] = {"value": float(clade.confidence)}
    if clade.clades:
        node["children"] = [clade_to_dict_kinesin(c, depth + 1) for c in clade.clades]
    return node


def clade_to_dict_kinesin14(clade, depth=0):
    label = clade.name or f"node_{depth}_{id(clade)}"
    node  = {
        "name": label,
        "node_attrs": {
            "div":            clade.branch_length or 0,
            "directionality": {"value": label_to_direction(label)},
            "species":        {"value": label_to_species(label)},
            "gene":           {"value": label_to_gene(label)},
        },
    }
    if clade.confidence is not None:
        node["node_attrs"]["bootstrap"] = {"value": float(clade.confidence)}
    if clade.clades:
        node["children"] = [clade_to_dict_kinesin14(c, depth + 1) for c in clade.clades]
    return node


# ── Auspice JSON builders ──────────────────────────────────────────
def build_kinesin_json(tree_file: str, output: str):
    tree = Phylo.read(tree_file, "newick")
    auspice = {
        "version": "v2",
        "meta": {
            "title": "Kinesin Motor Protein Superfamily Phylogeny",
            "description": (
                "Maximum-likelihood phylogeny of kinesin motor proteins "
                "(Kinesin-1 to Kinesin-14). Built with Biopython and RAxML. "
                "Reference: Duke Kinesin Site — https://sites.duke.edu/kinesin/kinesin-tree/"
            ),
            "colorings": [
                {"key": "family",         "title": "Kinesin Family",  "type": "categorical"},
                {"key": "species",        "title": "Species",         "type": "categorical"},
                {"key": "directionality", "title": "Motor Direction", "type": "categorical"},
                {"key": "gene",           "title": "Gene Name",       "type": "categorical"},
            ],
            "display_defaults": {"color_by": "family"},
            "panels": ["tree"],
        },
        "tree": clade_to_dict_kinesin(tree.root),
    }
    _write(auspice, output)


def build_kinesin14_json(tree_file: str, output: str):
    tree = Phylo.read(tree_file, "newick")
    auspice = {
        "version": "v2",
        "meta": {
            "title": "Kinesin-14 Family Phylogeny",
            "description": (
                "Phylogeny of the Kinesin-14 (C-terminal motor) family. "
                "Minus-end directed motors: DmNcd, ScKAR3, CgCHO2, AtKCBP. "
                "Reference: Duke Kinesin Site — https://sites.duke.edu/kinesin/the-kinesin-14-family/"
            ),
            "colorings": [
                {"key": "directionality", "title": "Motor Direction", "type": "categorical"},
                {"key": "species",        "title": "Species",         "type": "categorical"},
                {"key": "gene",           "title": "Gene Name",       "type": "categorical"},
            ],
            "display_defaults": {"color_by": "directionality"},
            "panels": ["tree"],
        },
        "tree": clade_to_dict_kinesin14(tree.root),
    }
    _write(auspice, output)


def _write(data: dict, output: str):
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w") as fh:
        json.dump(data, fh, indent=2)
    print(f"Auspice JSON → {output}")
    print("View at: https://auspice.us  (drag and drop the file)")


# ── Main ───────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Build Auspice v2 JSON from Newick tree")
    parser.add_argument("--tree",   required=True)
    parser.add_argument("--mode",   choices=["kinesin", "kinesin14"], required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if args.mode == "kinesin":
        build_kinesin_json(args.tree, args.output)
    else:
        build_kinesin14_json(args.tree, args.output)


if __name__ == "__main__":
    main()
