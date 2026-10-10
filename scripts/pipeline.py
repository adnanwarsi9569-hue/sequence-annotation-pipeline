#!/usr/bin/env python3
"""Sequence alignment & annotation pipeline.

FASTA -> BLASTP vs Swiss-Prot -> filter (e-value, identity, coverage)
      -> annotate best hit -> CSV + summary plots
"""
import argparse
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO

COLS = ["qseqid", "sseqid", "pident", "length", "qlen", "evalue", "bitscore", "stitle"]


def run_blast(query, db, evalue, threads, out):
    cmd = [
        "blastp", "-query", str(query), "-db", str(db),
        "-evalue", str(evalue), "-max_target_seqs", "10",
        "-num_threads", str(threads),
        "-outfmt", "6 " + " ".join(COLS), "-out", str(out),
    ]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_title(title):
    """Pull protein name, gene, organism, UniProt ID out of a Swiss-Prot header."""
    first = title.split()[0]
    uniprot = first.split("|")[1] if "|" in first else first
    rest = title[len(first):].strip()
    name = re.split(r"\sOS=", rest)[0]
    org = re.search(r"OS=(.*?)\sOX=", rest)
    gene = re.search(r"GN=(\S+)", rest)
    return pd.Series({
        "uniprot_id": uniprot,
        "annotation": name,
        "gene": gene.group(1) if gene else "NA",
        "organism": org.group(1) if org else "NA",
    })


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--query", required=True, help="Input FASTA")
    p.add_argument("--db", required=True, help="BLAST protein DB prefix")
    p.add_argument("--evalue", type=float, default=1e-5)
    p.add_argument("--identity", type=float, default=40.0, help="Min %% identity")
    p.add_argument("--coverage", type=float, default=50.0, help="Min %% query coverage")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--outdir", default="results")
    a = p.parse_args()

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    raw = out / "raw_blast.tsv"

    # 1. Alignment
    run_blast(a.query, a.db, a.evalue, a.threads, raw)
    hits = pd.read_csv(raw, sep="\t", names=COLS)
    print(f"Raw hits: {len(hits)}")

    # 2. Filtering
    hits["query_coverage"] = (hits["length"] / hits["qlen"] * 100).clip(upper=100).round(1)
    passed = hits[
        (hits["evalue"] <= a.evalue)
        & (hits["pident"] >= a.identity)
        & (hits["query_coverage"] >= a.coverage)
    ]
    print(f"Hits passing filters: {len(passed)}")

    # 3. Best hit per query + annotation
    best = (passed.sort_values(["evalue", "bitscore"], ascending=[True, False])
                  .groupby("qseqid", as_index=False).first())
    best = pd.concat([best, best["stitle"].apply(parse_title)], axis=1)

    # keep every query, even those without a good hit
    queries = pd.DataFrame(
        [(r.id, r.description) for r in SeqIO.parse(a.query, "fasta")],
        columns=["qseqid", "query_description"],
    )
    final = queries.merge(best, on="qseqid", how="left")
    final["status"] = np.where(final["sseqid"].isna(), "No significant hit", "Annotated")
    final = final[["qseqid", "query_description", "status", "uniprot_id", "gene",
                   "annotation", "organism", "pident", "query_coverage",
                   "evalue", "bitscore"]]
    final.to_csv(out / "annotated_hits.csv", index=False)
    print(f"Annotated {(final.status == 'Annotated').sum()} / {len(final)} sequences")

    # 4. Visualisation
    ann = final[final.status == "Annotated"]
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))

    ax[0].hist(ann["pident"], bins=15, color="#2a9d8f", edgecolor="white")
    ax[0].set(title="Identity of best hits", xlabel="% identity", ylabel="Sequences")

    logE = -np.log10(ann["evalue"].astype(float).clip(lower=1e-300))
    ax[1].hist(logE, bins=15, color="#e76f51", edgecolor="white")
    ax[1].set(title="Significance of best hits", xlabel="-log10(e-value)", ylabel="Sequences")

    def category(name):
        n = name.lower()
        if "tyrosine" in n:
            return "Tyrosine kinase"
        if "serine" in n or "threonine" in n:
            return "Ser/Thr kinase"
        if any(k in n for k in ["lipid", "glycerol", "inositol", "phosphatidyl"]):
            return "Lipid kinase"
        return "Other / metabolic"

    cats = ann["annotation"].apply(category).value_counts().sort_values()
    ax[2].barh(cats.index, cats.values, color="#264653")
    ax[2].set(title="Annotated function classes", xlabel="Sequences")

    plt.tight_layout()
    plt.savefig(out / "summary_plots.png", dpi=200)
    print("Saved: raw_blast.tsv, annotated_hits.csv, summary_plots.png ->", out)


if __name__ == "__main__":
    main()
