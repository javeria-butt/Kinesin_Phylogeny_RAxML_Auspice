"""
run_raxml.py
────────────
Run phylogenetic tree inference.
Tries RAxML → FastTree → Biopython NJ (in that order).

Usage:
    python scripts/run_raxml.py \
        --alignment data/alignments/kinesin_aligned.fasta \
        --prefix    kinesin \
        --outdir    results/raxml_kinesin \
        --bootstrap 100
"""

import argparse
import shutil
import subprocess
from pathlib import Path
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


def run_raxml(aln: str, prefix: str, outdir: str, n_boot: int, model: str) -> tuple:
    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)

    # Remove old RAxML files with the same prefix
    for f in outdir_p.glob(f"RAxML_*{prefix}*"):
        f.unlink()

    raxml_bin = "raxmlHPC" if shutil.which("raxmlHPC") else "raxml"
    cmd = (
        f"{raxml_bin} -f a "
        f"-m {model} "
        f"-p 42 -x 42 "
        f"-# {n_boot} "
        f"-s {aln} "
        f"-n {prefix} "
        f"-w {outdir_p.resolve()}"
    )
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    tree_file = outdir_p / f"RAxML_bipartitions.{prefix}"
    if result.returncode != 0 or not tree_file.exists():
        raise RuntimeError(result.stderr)
    return str(tree_file), "RAxML"


def run_fasttree(aln: str, tree_file: str) -> tuple:
    ft = "FastTree" if shutil.which("FastTree") else "fasttree"
    cmd = f"{ft} -lg -gamma -quiet {aln} > {tree_file}"
    subprocess.run(cmd, shell=True, check=True)
    return tree_file, "FastTree (LG+Gamma)"


def run_nj(aln: str, tree_file: str) -> tuple:
    alignment = AlignIO.read(aln, "fasta")
    calc  = DistanceCalculator("blosum62")
    ctor  = DistanceTreeConstructor(calc, "nj")
    tree  = ctor.build_tree(alignment)
    Phylo.write(tree, tree_file, "newick")
    return tree_file, "Neighbour-Joining / BLOSUM62 (Biopython)"


def infer_tree(aln: str, prefix: str, outdir: str, n_boot: int, model: str) -> tuple:
    final_nwk = str(Path(outdir) / f"{prefix}_tree.nwk")

    if shutil.which("raxmlHPC") or shutil.which("raxml"):
        try:
            tree_file, method = run_raxml(aln, prefix, outdir, n_boot, model)
            shutil.copy(tree_file, final_nwk)
            return final_nwk, method
        except RuntimeError as e:
            print(f"RAxML failed: {e} — falling back to FastTree")

    if shutil.which("FastTree") or shutil.which("fasttree"):
        return run_fasttree(aln, final_nwk)

    return run_nj(aln, final_nwk)


def main():
    parser = argparse.ArgumentParser(description="Infer phylogenetic tree")
    parser.add_argument("--alignment",  required=True)
    parser.add_argument("--prefix",     default="kinesin")
    parser.add_argument("--outdir",     default="results")
    parser.add_argument("--bootstrap",  type=int, default=100)
    parser.add_argument("--model",      default="PROTGAMMALG",
                        help="RAxML substitution model")
    args = parser.parse_args()

    tree_file, method = infer_tree(
        args.alignment, args.prefix, args.outdir, args.bootstrap, args.model
    )

    tree = Phylo.read(tree_file, "newick")
    print(f"Method    : {method}")
    print(f"Leaves    : {tree.count_terminals()}")
    print(f"Int nodes : {len(tree.get_nonterminals())}")
    print(f"Tree file : {tree_file}")


if __name__ == "__main__":
    main()
