#!/usr/bin/env python3
r"""
fc_checker.py – First‑Classness (FC) audit with progress bar + null model
-----------------------------------------------------------------------
$ py fc_checker.py genbank_species.tsv refseq_species.tsv \
      --map orbit_map.csv --plots --null 10000 --progress 2000

Flags
-----
--map FILE        two‑column CSV, no header (codon,orbit).  Default: orbit_map.csv
--progress INT    print a heartbeat every INT rows             [default 1000]
--plots           write violin + null‑hist PNGs to ./fc_out
--null [N]        add shuffled‑orbit null model; optional N    [default 10000]
--quiet           suppress per‑row warnings/NaN notes

Requires numpy pandas scipy matplotlib
"""

import argparse, math, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
import matplotlib.pyplot as plt

CODON_ORDER = [a + b + c for a in "UCAG" for b in "UCAG" for c in "UCAG"]

# ---------------------------------------------------------------------------

def load_map(path):
    df = pd.read_csv(path, header=None, names=["codon", "orbit"])
    if len(df) != 64 or len(set(df.codon)) != 64:
        raise ValueError("orbit_map.csv must list each RNA codon exactly once")
    return dict(zip(df.codon, df.orbit.astype(int)))


def read_cutg(path):
    df = pd.read_csv(path, sep="\t", header=0, engine="python", on_bad_lines="skip")
    df = df.iloc[:, [0] + list(range(df.shape[1] - 64, df.shape[1]))]
    df.rename(columns=lambda c: c.replace("T", "U") if len(c) == 3 else c, inplace=True)
    df.columns = ["species"] + CODON_ORDER
    df[CODON_ORDER] = df[CODON_ORDER].apply(pd.to_numeric, errors="coerce").fillna(0)
    return df


def to_rscu(counts, degen):
    total = counts.sum()
    if total == 0:
        return np.full(64, np.nan)
    freq = counts / total
    return np.array([freq[i] * degen[CODON_ORDER[i]] for i in range(64)])


def sigma_ratio(rscu, orbit):
    groups = {}
    for v, c in zip(rscu, CODON_ORDER):
        groups.setdefault(orbit[c], []).append(v)
    intra = [x - np.mean(v) for v in groups.values() for x in v]
    inter = [np.mean(v) for v in groups.values() for _ in v]
    den = np.std(inter, ddof=1)
    return np.nan if den == 0 else np.std(intra, ddof=1) / den

# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tables", nargs="+", help="CUTG species TSV files")
    ap.add_argument("--map", default="orbit_map.csv")
    ap.add_argument("--plots", action="store_true")
    ap.add_argument("--null", nargs="?", const=10000, type=int, default=0,
                    help="add null model with optional shuffle count [10000]")
    ap.add_argument("--progress", type=int, default=1000)
    ap.add_argument("--out", default="fc_out")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    Path(args.out).mkdir(exist_ok=True)
    ORBIT = load_map(args.map)
    DEGEN = {c: ORBIT[c] for c in CODON_ORDER}

    ratios = []
    rows = 0
    for tbl in args.tables:
        print(f"Reading {tbl} …")
        for _, row in read_cutg(tbl).iterrows():
            rows += 1
            if rows % args.progress == 0:
                print(f"  processed {rows:,} rows …", flush=True)
            rscu = to_rscu(row[1:], DEGEN)
            if np.isnan(rscu).any():
                if not args.quiet:
                    print(f"    skipped NaN row {row[0]}")
                continue
            r = sigma_ratio(rscu, ORBIT)
            if math.isfinite(r):
                ratios.append(r)

    if not ratios:
        sys.exit("No valid rows – check orbit map and input tables.")

    med = float(np.median(ratios))
    print(f"Median σ_intra/σ_inter = {med:.3f}  (n = {len(ratios):,})")

    if args.plots:
        plt.figure(figsize=(2.4, 4))
        plt.violinplot(ratios, showmedians=True)
        plt.ylabel("σ_intra / σ_inter"); plt.tight_layout()
        plt.savefig(Path(args.out, "violin_sigma_ratio.png"), dpi=300)
        print("Violin saved →", Path(args.out, "violin_sigma_ratio.png"))

    if args.null:
        rng = np.random.default_rng(2025)
        null = []
        base_orbits = list(ORBIT.values())
        for _ in range(args.null):
            rng.shuffle(base_orbits)
            fake = dict(zip(CODON_ORDER, base_orbits))
            groups = {}
            for c in CODON_ORDER:
                groups.setdefault(fake[c], []).append(1.0)
            intra = [x - np.mean(v) for v in groups.values() for x in v]
            inter = [np.mean(v) for v in groups.values() for _ in v]
            d = np.std(inter, ddof=1)
            null.append(np.std(intra, ddof=1) / d if d else np.nan)
        null = [n for n in null if math.isfinite(n)]
        z = (med - np.mean(null)) / np.std(null)
        p = (sum(n <= med for n in null) + 1) / (len(null) + 1)
        print(f"Null‑model Z = {z:.2f},  empirical P = {p:.4g}")
        if args.plots:
            plt.figure(figsize=(3.5, 3))
            plt.hist(null, bins=40, alpha=0.7)
            plt.axvline(med, color="red"); plt.xlabel("median σ_intra/σ_inter (null)")
            plt.tight_layout(); name = Path(args.out, "null_hist.png")
            plt.savefig(name, dpi=300); print("Null hist →", name)

if __name__ == "__main__":
    main()
