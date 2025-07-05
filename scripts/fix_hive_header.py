#!/usr/bin/env python3
"""
Supplementary Script S1: fix_hive_header.py

Problem:
    The CUTG HIVE dump (e.g. refseq_cds.tsv) has a malformed header
    in which adjacent codon names have accidentally merged. All
    downstream parsing shifts columns left, so Taxid/Organelle get
    overwritten by species strings, etc.

Solution:
    Read the original file, skip its first (broken) header line, and
    write out a new file whose first line is the correct, tab-delimited
    list of 64 codon columns plus metadata fields.

Usage:
    python fix_hive_header.py refseq_cds.tsv 
Produces:
    refseq_cds_fixed.tsv
"""

import sys
from pathlib import Path
import pandas as pd
 
# The exact, correct header line for the CUTG refseq_cds dump:

FIXED_HEADER = [
    "Division", "Assembly", "Taxid", "Species", "Organelle",
    "Translation Table", "# CDS", "# Codons",
    "GC%", "GC1%", "GC2%", "GC3%",
    # the 64 DNA codons in the proper order:
    "TTT","TTC","TTA","TTG",
    "CTT","CTC","CTA","CTG",
    "ATT","ATC","ATA","ATG",
    "GTT","GTC","GTA","GTG",
    "TAT","TAC","TAA","TAG",
    "CAT","CAC","CAA","CAG",
    "AAT","AAC","AAA","AAG",
    "GAT","GAC","GAA","GAG",
    "TCT","TCC","TCA","TCG",
    "CCT","CCC","CCA","CCG",
    "ACT","ACC","ACA","ACG",
    "GCT","GCC","GCA","GCG",
    "TGT","TGC","TGA","TGG",
    "CGT","CGC","CGA","CGG",
    "AGT","AGC","AGA","AGG",
    "GGT","GGC","GGA","GGG",
]
 
def main(raw_path: Path):
    if not raw_path.exists():
        print(f"PROBLEM: File not found: {raw_path}")
        sys.exit(f"File not found: {raw_path}")

    print(f"entered Main")

    fixed_path = raw_path.with_name(raw_path.stem + "_fixed.tsv")
    print(f"Reading broken file: {raw_path}")
    print(f"Writing fixed file:  {fixed_path}")

    body = pd.read_csv(raw_path, sep="\t", skiprows=1, header=None,
                       dtype=str, low_memory=False)

    # Drop the extraneous final column (column index 76):
    if len(body.columns) == len(FIXED_HEADER) + 1:
        print("Extra column detected, dropping the last column...")
        body = body.drop(columns=[76])

    # Ensure columns now match the FIXED_HEADER length:
    if len(body.columns) != len(FIXED_HEADER):
        print(f"Final column mismatch detected after adjustment. Data columns: {len(body.columns)}, Header columns: {len(FIXED_HEADER)}")
        sys.exit("Column mismatch persists after adjustment. Check input data manually.")

    body.columns = FIXED_HEADER
    body = body[FIXED_HEADER]

    body.to_csv(fixed_path, sep="\t", index=False, header=FIXED_HEADER)

    print("Header fixed. Downstream workflows should now use:")
    print(f"    {fixed_path}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"Usage: python {sys.argv[0]} <broken_tsv_file>")
    main(Path(sys.argv[1]))

