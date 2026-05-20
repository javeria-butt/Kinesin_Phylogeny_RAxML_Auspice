"""
run_alignment.py
────────────────
Run multiple sequence alignment on a FASTA file.
Tries MUSCLE → MAFFT → Biopython (in that order).

Usage:
    python scripts/run_alignment.py \
        --input  data/sequences/kinesin_all_families.fasta \
        --output data/alignments/kinesin_aligned.fasta
"""

import argparse
import shutil
import subprocess
from pathlib import Path
from Bio import SeqIO, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator


def align_muscle(in_fasta: str, out_fasta: str) -> str:
    """MUSCLE v5 alignment."""
    cmd = f"muscle -align {in_fasta} -output {out_fasta}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(result.stderr)
    return "MUSCLE v5"


def align_mafft(in_fasta: str, out_fasta: str) -> str:
    """MAFFT --auto alignment."""
    cmd = f"mafft --auto {in_fasta} > {out_fasta}"
    subprocess.run(cmd, shell=True, check=True)
    return "MAFFT (--auto)"


def align_biopython(in_fasta: str, out_fasta: str) -> str:
    """
    No-op 'alignment' using Biopython — sequences written as-is.
    Only suitable for short demo runs; not for publication.
    """
    recs = list(SeqIO.parse(in_fasta, "fasta"))
    SeqIO.write(recs, out_fasta, "fasta")
    return "Biopython identity (no aligner available — install MUSCLE)"


def run_alignment(in_fasta: str, out_fasta: str) -> str:
    Path(out_fasta).parent.mkdir(parents=True, exist_ok=True)

    if shutil.which("muscle"):
        method = align_muscle(in_fasta, out_fasta)
    elif shutil.which("mafft"):
        method = align_mafft(in_fasta, out_fasta)
    else:
        method = align_biopython(in_fasta, out_fasta)

    aln = AlignIO.read(out_fasta, "fasta")
    print(f"Method          : {method}")
    print(f"Sequences       : {len(aln)}")
    print(f"Alignment length: {aln.get_alignment_length()} columns")
    print(f"Output          : {out_fasta}")
    return method


def main():
    parser = argparse.ArgumentParser(description="Run multiple sequence alignment")
    parser.add_argument("--input",  required=True, help="Input FASTA")
    parser.add_argument("--output", required=True, help="Output aligned FASTA")
    args = parser.parse_args()
    run_alignment(args.input, args.output)


if __name__ == "__main__":
    main()
