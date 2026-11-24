"""Step 0. Download the input data into <workdir>/data.

Downloads the 8 H3K9me2 bedGraphs listed in <workdir>/metadata/samples.tsv from GEO, and the GENCODE vM25 GTF.
Each GEO supplementary file lives at:
https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1725nnn/GSM1725725/suppl/<file>

Usage:
    python download_data.py
    python download_data.py --workdir /other/path

Files already present (compressed or uncompressed) are skipped, so the script can be rerun
after an interruption.
"""
import argparse
import os
import shutil
import urllib.request

import pandas as pd

WORKDIR = "/Users/adrian/research/bmcbf/102"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/samples/{prefix}nnn/{gsm}/suppl/{file}"
GTF = ("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/"
       "gencode.vM25.annotation.gtf.gz")


def fetch(url, dest):
    if os.path.exists(dest) or os.path.exists(dest.removesuffix(".gz")):
        print(f"skip (exists): {dest}")
        return
    print(f"downloading {url}")
    tmp = dest + ".part"
    with urllib.request.urlopen(url) as r, open(tmp, "wb") as out:
        shutil.copyfileobj(r, out)
    os.replace(tmp, dest)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--workdir", default=WORKDIR)
    a = p.parse_args()

    data = os.path.join(a.workdir, "data")
    os.makedirs(data, exist_ok=True)
    os.makedirs(os.path.join(a.workdir, "results"), exist_ok=True)

    samples = pd.read_csv(os.path.join(a.workdir, "metadata", "samples.tsv"), sep="\t")
    for name in samples["file"]:
        if not name.endswith(".gz"):
            name += ".gz"  # GEO only has the compressed files
        gsm = name.split("_")[0]
        fetch(GEO.format(prefix=gsm[:-3], gsm=gsm, file=name), os.path.join(data, name))

    fetch(GTF, os.path.join(data, os.path.basename(GTF)))


if __name__ == "__main__":
    main()
