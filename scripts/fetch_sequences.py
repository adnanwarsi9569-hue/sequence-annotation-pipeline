from Bio import Entrez, SeqIO
import re, sys

Entrez.email = "adnanwarsi9569@gmail.com"
term = "Homo sapiens[Organism] AND RefSeq[Filter] AND kinase[Title] AND 300:700[SLEN]"
h = Entrez.esearch(db="protein", term=term, retmax=400)
ids = Entrez.read(h)["IdList"]
h.close()
print(f"Found {len(ids)} IDs", file=sys.stderr)

h = Entrez.efetch(db="protein", id=",".join(ids), rettype="fasta", retmode="text")
records = list(SeqIO.parse(h, "fasta"))
h.close()

# keep one isoform per protein name
seen, unique = set(), []
for r in records:
    name = r.description.split(" ", 1)[1]
    name = re.sub(r",?\s*(isoform \w+|precursor)", "", name).lower()
    if name not in seen:
        seen.add(name)
        unique.append(r)
unique = unique[:40]

SeqIO.write(unique, "data/query_sequences.fasta", "fasta")
print(f"Saved {len(unique)} unique proteins", file=sys.stderr)
