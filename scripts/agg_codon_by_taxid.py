#!/usr/bin/env python3
r"""agg_codon_by_taxid.py – aggregate *fixed* CDS tables
=======================================================
Collapses a header‑fixed RefSeq or GenBank CDS TSV into **one row per
(Taxid, Organelle)** containing

* the 64 DNA‑codon columns
* `#CDS` – number of CDS rows aggregated
* `#Codons` – total codons (row‑wise sum of the 64 counts)

Quick usage
-----------
    python agg_codon_by_taxid.py refseq  ACF\CUTG               
    python agg_codon_by_taxid.py genBank ACF\CUTG --progress 50k

The script prints a ✓ when writing the final file and, when `--progress`
is given, an updating counter every *N* rows so huge inputs do not look
stall‑bound.

CLI reference
-------------
<dataset>      "refseq" | "genBank" (case‑insensitive)
[root]         project root (default: current directory)
--cds-dir DIR  folder with *fixed* CDS tables   (default: "CDS")
--agg-dir DIR  folder for aggregated output     (default: "AGG")
--progress N   emit a progress line every **N** input rows
               value may use k/M suffixes: 10k, 2M etc.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

# ------------------------------------------------------------
# helpers
# ------------------------------------------------------------

def _resolve_dir(raw: str | Path, root: Path) -> Path:
    p = Path(raw)
    return p if p.is_absolute() else (root / p)


def _find_input(dataset: str, cds_dir: Path) -> Path:
    tag = dataset.lower()
    cand = [f"{tag}_species_cds_fixed.tsv", f"{tag}_cds_fixed.tsv",
            f"{tag}_species_cds.tsv", f"{tag}_cds.tsv"]
    for name in cand:
        p = cds_dir / name
        if p.exists():
            return p
    raise FileNotFoundError(f"None of {', '.join(cand)} found in {cds_dir}")


def _codon_columns() -> List[str]:
    bases = "TCAG"  # lexicographic order used by CUTG
    return [a + b + c for a in bases for b in bases for c in bases]


def _parse_int_with_suffix(text: str) -> int:
    """Accept  "500", "10k", "2M", etc., return integer."""
    m = re.fullmatch(r"(\d+)([kKmM]?)", text)
    if not m:
        raise ValueError(text)
    num, suf = m.groups()
    n = int(num)
    if suf in {"k", "K"}:
        n *= 1_000
    elif suf in {"m", "M"}:
        n *= 1_000_000
    return n

# ------------------------------------------------------------
# main
# ------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Aggregate (Taxid, Organelle) codon counts with optional progress feedback.")
    ap.add_argument("dataset", choices=["refseq", "genBank", "genbank", "RefSeq"],
                    help="Which dataset to process (refseq | genBank)")
    ap.add_argument("root", nargs="?", default=".", help="Project root (default: current dir)")
    ap.add_argument("--cds-dir", default="CDS", help="Folder containing fixed CDS TSV files")
    ap.add_argument("--agg-dir", default="AGG", help="Folder to write aggregated output")
    ap.add_argument("--progress", metavar="N", default=None,
                    help="Print progress every N rows (accepts k/M suffix)")

    args = ap.parse_args()

    root    = Path(args.root).resolve()
    cds_dir = _resolve_dir(args.cds_dir, root)
    agg_dir = _resolve_dir(args.agg_dir, root)
    agg_dir.mkdir(parents=True, exist_ok=True)

    dataset  = args.dataset.lower()
    raw_path = _find_input(dataset, cds_dir)
    out_path = agg_dir / f"{dataset}_codon_species.tsv"

    step = None
    if args.progress is not None:
        try:
            step = _parse_int_with_suffix(args.progress)
            if step <= 0:
                step = None
        except ValueError:
            ap.error("--progress expects integer or k/M‑suffix value, e.g. 10000 or 10k")

    codon_cols = _codon_columns()
    counts: Dict[Tuple[str, str], Dict[str, int]] = {}

    # ---------- aggregation loop ----------
    row_no = 0
    with raw_path.open(newline="") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for row in rdr:
            key = (row["Taxid"], row["Organelle"])
            acc = counts.setdefault(key, {c: 0 for c in codon_cols})
            acc.setdefault("#CDS", 0)
            for c in codon_cols:
                acc[c] += int(row[c])
            acc["#CDS"] += 1

            row_no += 1
            if step and row_no % step == 0:
                print(f"…{row_no:,} rows", file=sys.stderr)

    # ---------- write result ----------
    with out_path.open("w", newline="") as f:
        fieldnames = ["Taxid", "Organelle", *codon_cols, "#CDS", "#Codons"]
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        w.writeheader()
        for (taxid, org), acc in counts.items():
            row = {c: acc[c] for c in codon_cols}
            row["#CDS"]    = acc["#CDS"]
            row["#Codons"] = sum(acc[c] for c in codon_cols)
            w.writerow({"Taxid": taxid, "Organelle": org, **row})

    print(f"\u2713 wrote {out_path.relative_to(root)}")


if __name__ == "__main__":
    main()
