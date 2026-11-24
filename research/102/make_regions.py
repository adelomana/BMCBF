"""Step 1. Build gene-body and promoter BED files from the GENCODE vM25 GTF (mm10 / GRCm38).

Input:  <workdir>/data/gencode.vM25.annotation.gtf.gz
Output: <workdir>/results/gene_bodies.bed and <workdir>/results/promoters.bed
        (BED6: chrom, start, end, gene_id|gene_name, 0, strand)

Only protein-coding genes on autosomes are kept (chrX, chrY, chrM dropped, see README).

Usage:
    python make_regions.py
    python make_regions.py --promoter-flank 2500 --workdir /other/path
"""
import argparse
import gzip
import os
import re

WORKDIR = "/Users/adrian/research/bmcbf/102"
GTF = "gencode.vM25.annotation.gtf.gz"
KEEP = {f"chr{i}" for i in range(1, 20)}


def attr(field, key):
    m = re.search(rf'{key} "([^"]+)"', field)
    return m.group(1) if m else None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--workdir", default=WORKDIR)
    p.add_argument("--promoter-flank", type=int, default=2500)
    p.add_argument("--gene-type", default="protein_coding")
    a = p.parse_args()

    gtf = os.path.join(a.workdir, "data", GTF)
    if not os.path.exists(gtf):
        gtf = gtf.removesuffix(".gz")
    results = os.path.join(a.workdir, "results")
    os.makedirs(results, exist_ok=True)

    opener = gzip.open if gtf.endswith(".gz") else open
    bodies, proms = [], []
    with opener(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "gene" or f[0] not in KEEP:
                continue
            if attr(f[8], "gene_type") != a.gene_type:
                continue
            chrom, start, end, strand = f[0], int(f[3]) - 1, int(f[4]), f[6]
            name = f"{attr(f[8], 'gene_id')}|{attr(f[8], 'gene_name')}"
            bodies.append((chrom, start, end, name, 0, strand))
            tss = start if strand == "+" else end
            proms.append((chrom, max(0, tss - a.promoter_flank), tss + a.promoter_flank, name, 0, strand))

    for fname, rows in (("gene_bodies.bed", bodies), ("promoters.bed", proms)):
        rows.sort(key=lambda r: (r[0], r[1]))
        path = os.path.join(results, fname)
        with open(path, "w") as out:
            for r in rows:
                out.write("\t".join(map(str, r)) + "\n")
        print(f"{path}: {len(rows)} regions")


if __name__ == "__main__":
    main()
