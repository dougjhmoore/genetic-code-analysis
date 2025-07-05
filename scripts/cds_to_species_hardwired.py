#!/usr/bin/env python3
r"""cds_to_species_hardwired.py
---------------------------------
Hard‑wired one‑shot conversion of your CUTG *CDS* tables into per‑species
summary tables.  **No command‑line arguments needed.**  Just run:

    python cds_to_species_hardwired.py

Prerequisites
-------------
* pandas
* tqdm  (pip install tqdm)

What it does
------------
1. Loads the two big CDS tables located in
   C:\Users\djhmo\OneDrive\Projects\ACF\CUTG\CDS\
   * genbank_species_cds.tsv
   * refseq_species_cds.tsv

2. Streams them in batches of 50 000 rows, coercing blanks to 0.
3. Sums every numeric column by the "Species" field.
4. Writes:
   * genbank_species.tsv
   * refseq_species.tsv
   in **the same folder** as the inputs.

Each run prints a tqdm progress bar + a final row count.
"""

import pandas as pd
from pathlib import Path
from tqdm import tqdm

# -------------------------------------------------------------
BASE = Path(r"C:\Users\djhmo\OneDrive\Projects\ACF\CUTG\CDS")
IN_FILES = [BASE / "genbank_species_cds.tsv", BASE / "refseq_species_cds.tsv"]
CHUNK = 50_000  # rows per chunk

CODONS_DNA = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
CODONS_RNA = [c.replace("T", "U") for c in CODONS_DNA]

# mapping DNA->RNA column names once
DNA2RNA = dict(zip(CODONS_DNA, CODONS_RNA))

for in_path in IN_FILES:
    if not in_path.exists():
        print(f"ERROR: {in_path} not found – skipping.")
        continue

    print(f"Processing {in_path.name} …")
    agg = None  # will hold the running per‑species totals

    for chunk in tqdm(pd.read_csv(in_path, sep="\t", engine="python", chunksize=CHUNK),
                      unit="rows"):
        # normalise column names if needed
        if "Species" not in chunk.columns:
            if "SpeciesName" in chunk.columns:
                chunk.rename(columns={"SpeciesName": "Species"}, inplace=True)
            else:
                raise RuntimeError("Input file missing 'Species' column")
        if "# Codons" not in chunk.columns and "TotalCodons" in chunk.columns:
            chunk.rename(columns={"TotalCodons": "# Codons"}, inplace=True)

        # convert DNA columns to RNA and coerce blanks to 0
        chunk.rename(columns=DNA2RNA, inplace=True)
        chunk[CODONS_RNA] = chunk[CODONS_RNA].apply(pd.to_numeric, errors="coerce").fillna(0)

        # aggregate this batch
        batch_sum = chunk.groupby("Species", as_index=False).sum(numeric_only=True)
        agg = batch_sum if agg is None else pd.concat([agg, batch_sum]).groupby("Species", as_index=False).sum()

    # finished streaming – write output
    out_path = in_path.with_name(in_path.name.replace("_cds", "_species"))
    agg.to_csv(out_path, sep="\t", index=False)
    print(f"Written {len(agg):,} rows → {out_path.name}\n")
