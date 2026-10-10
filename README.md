# Sequence Alignment & Annotation Pipeline

Python pipeline that aligns protein sequences (FASTA) against UniProt Swiss-Prot
using BLASTP, filters hits by e-value, identity and query coverage, annotates each
sequence with its gene function, and writes a CSV plus summary plots.

Built for the OpenBluff Bioinformatics Internship (Task 1).

## Workflow

FASTA -> BLASTP vs Swiss-Prot -> filters -> best hit per query -> annotation -> CSV + plots

## Requirements

- Python 3.8+
- NCBI BLAST+ (sudo apt install ncbi-blast+)
- pip install -r requirements.txt

## Usage

    bash scripts/setup_db.sh
    python scripts/fetch_sequences.py
    python scripts/pipeline.py --query data/query_sequences.fasta --db data/sprot_db \
      --evalue 1e-5 --identity 40 --coverage 50 --outdir results

| Argument | Default | Description |
|---|---|---|
| --query | required | Input FASTA |
| --db | required | BLAST database prefix |
| --evalue | 1e-5 | Max e-value |
| --identity | 40 | Min percent identity |
| --coverage | 50 | Min percent query coverage |
| --threads | 4 | BLAST threads |
| --outdir | results | Output folder |

## Sample dataset

40 unique human kinase proteins (one isoform each) from NCBI RefSeq,
downloaded with scripts/fetch_sequences.py.

## Outputs

- results/raw_blast.tsv : raw BLAST hits
- results/annotated_hits.csv : best hit per query (gene, function, organism, identity, coverage, e-value, bitscore)
- results/summary_plots.png : identity, e-value and function-class plots

## Sample run results

- 433 raw hits, 298 passed filters (e-value <= 1e-5, identity >= 40%, coverage >= 50%)
- 40 / 40 sequences annotated, 40 unique genes

## Author

Adnan Warsi
