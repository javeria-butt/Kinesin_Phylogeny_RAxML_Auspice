"""
fetch_sequences.py
──────────────────
Fetch kinesin protein sequences from NCBI and save to FASTA.

Usage:
    python scripts/fetch_sequences.py --family all   --out data/sequences/kinesin_all_families.fasta
    python scripts/fetch_sequences.py --family kin14 --out data/sequences/kinesin14_family.fasta
"""

import argparse
import time
from pathlib import Path
from Bio import Entrez, SeqIO

# ── Set your email — required by NCBI ─────────────────────────────
Entrez.email = "user@example.com"

# ── Accession tables ───────────────────────────────────────────────
ALL_ACCESSIONS = [
    ("NP_004983", "Kinesin-1",  "Hs", "KIF5A"),
    ("NP_004521", "Kinesin-1",  "Hs", "KIF5B"),
    ("NP_004954", "Kinesin-1",  "Hs", "KIF5C"),
    ("NP_524153", "Kinesin-1",  "Dm", "KHC"),
    ("NP_499043", "Kinesin-1",  "Ce", "UNC-116"),
    ("NP_009376", "Kinesin-2",  "Hs", "KIF3A"),
    ("NP_000219", "Kinesin-2",  "Hs", "KIF3B"),
    ("NP_524582", "Kinesin-2",  "Dm", "KLP64D"),
    ("NP_004270", "Kinesin-3",  "Hs", "KIF1A"),
    ("NP_055046", "Kinesin-3",  "Hs", "KIF1B"),
    ("NP_493636", "Kinesin-3",  "Ce", "UNC-104"),
    ("NP_002253", "Kinesin-4",  "Hs", "KIF4A"),
    ("NP_524914", "Kinesin-4",  "Dm", "KLP3A"),
    ("NP_004514", "Kinesin-5",  "Hs", "KIF11"),
    ("NP_524907", "Kinesin-5",  "Dm", "KLP61F"),
    ("NP_013710", "Kinesin-5",  "Sc", "CIN8"),
    ("NP_006255", "Kinesin-6",  "Hs", "KIF23"),
    ("NP_492484", "Kinesin-6",  "Ce", "ZEN-4"),
    ("NP_001800", "Kinesin-7",  "Hs", "KIF10"),
    ("NP_055732", "Kinesin-8",  "Hs", "KIF18A"),
    ("NP_013490", "Kinesin-8",  "Sc", "KIP3"),
    ("NP_006058", "Kinesin-13", "Hs", "KIF2A"),
    ("NP_006836", "Kinesin-13", "Hs", "KIF2C"),
    ("NP_477177", "Kinesin-14", "Dm", "NCD"),
    ("NP_013491", "Kinesin-14", "Sc", "KAR3"),
    ("NP_003147", "Kinesin-14", "Hs", "KIFC1"),
    ("NP_055854", "Kinesin-14", "Hs", "KIFC2"),
    ("NP_009118", "Kinesin-14", "Hs", "KIFC3"),
]

KIN14_ACCESSIONS = [
    ("NP_477177", "Kinesin-14", "Dm", "NCD"),
    ("NP_013491", "Kinesin-14", "Sc", "KAR3"),
    ("NP_198067", "Kinesin-14", "At", "KCBP"),
    ("AAA35501",  "Kinesin-14", "Cg", "CHO2"),
    ("NP_003147", "Kinesin-14", "Hs", "KIFC1"),
    ("NP_055854", "Kinesin-14", "Hs", "KIFC2"),
    ("NP_009118", "Kinesin-14", "Hs", "KIFC3"),
    ("NP_492310", "Kinesin-14", "Ce", "KLP-15"),
    ("NP_498907", "Kinesin-14", "Ce", "KLP-16"),
    ("NP_510926", "Kinesin-14", "Ce", "KLP-17"),
    ("NP_172918", "Kinesin-14", "At", "KATA"),
    ("NP_001031", "Kinesin-14", "At", "KATB"),
    ("NP_187943", "Kinesin-14", "At", "KATC"),
]


def fetch(accessions, output_fasta, delay=0.4):
    """Fetch sequences from NCBI and write to FASTA."""
    records = []
    failed  = []
    for acc, family, species, gene in accessions:
        label = f"{species}_{gene}_{family.replace('-', '')}"
        try:
            handle = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
            rec    = SeqIO.read(handle, "fasta")
            handle.close()
            rec.id          = label
            rec.description = ""
            records.append(rec)
            print(f"  OK  {acc:15s}  {label}")
        except Exception as exc:
            failed.append(acc)
            print(f"  FAIL {acc:15s}  {exc}")
        time.sleep(delay)

    Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, output_fasta, "fasta")
    print(f"\n{len(records)} sequences → {output_fasta}")
    if failed:
        print(f"Failed: {failed}")
    return records


def main():
    parser = argparse.ArgumentParser(description="Fetch kinesin sequences from NCBI")
    parser.add_argument("--family", choices=["all", "kin14"], default="all",
                        help="Which family to fetch")
    parser.add_argument("--out", required=True, help="Output FASTA path")
    parser.add_argument("--delay", type=float, default=0.4,
                        help="Delay between NCBI requests (seconds)")
    args = parser.parse_args()

    table = ALL_ACCESSIONS if args.family == "all" else KIN14_ACCESSIONS
    fetch(table, args.out, delay=args.delay)


if __name__ == "__main__":
    main()
