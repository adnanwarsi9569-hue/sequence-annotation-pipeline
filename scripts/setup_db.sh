#!/usr/bin/env bash
# Downloads Swiss-Prot and builds a local BLAST protein database
set -e
mkdir -p data && cd data
wget -nc https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip -kf uniprot_sprot.fasta.gz
makeblastdb -in uniprot_sprot.fasta -dbtype prot -out sprot_db
