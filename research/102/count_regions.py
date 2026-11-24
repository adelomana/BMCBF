"""Step 2. Approximate read counts per region from the GEO coverage bedGraphs.

The GEO bedGraphs are integer pileups of reads extended to the fragment length
(300 nt for H3K9me2 native ChIP, per the paper's Methods). So for any region:

    reads ~= sum(coverage * overlap_bp) / fragment_length

Counts are rounded so DESeq2 can use them. The script also writes the approximate
total reads per sample (whole bedGraph), used for library-size normalization in Step 3.

Input:  <workdir>/metadata/samples.tsv, the bedGraphs in <workdir>/data (compressed or not),
        <workdir>/results/<regions>.bed from Step 1
Output: <workdir>/results/counts_<regions>.tsv and <workdir>/results/counts_<regions>_totals.tsv

Usage:
    python count_regions.py gene_bodies
    python count_regions.py promoters
    python count_regions.py gene_bodies --workdir /other/path
"""
import argparse
import os

import numpy as np
import pandas as pd

WORKDIR = "/Users/adrian/research/bmcbf/102"


def find_bedgraph(data, name):
    for candidate in (name, name.removesuffix(".gz"), name + ".gz"):
        path = os.path.join(data, candidate)
        if os.path.exists(path):
            return path
    raise FileNotFoundError(f"{name} not found in {data}, run download_data.py first")


def load_bedgraph(path):
    return pd.read_csv(
        path, sep=r"\s+", header=None, names=["chrom", "start", "end", "cov"],
        dtype={"chrom": "category", "start": np.int32, "end": np.int32, "cov": np.float32},
    )


def region_area(bg_chrom, starts, ends):
    """Sum of coverage * bp over [starts, ends) for one chromosome.

    bg_chrom must be non-overlapping intervals (true for bedGraph pileups).
    """
    s = bg_chrom["start"].to_numpy()
    e = bg_chrom["end"].to_numpy()
    v = bg_chrom["cov"].to_numpy()
    order = np.argsort(s)
    s, e, v = s[order], e[order], v[order]
    cum = np.concatenate([[0.0], np.cumsum(v.astype(np.float64) * (e - s))])

    def area_upto(x):
        k = np.searchsorted(s, x, side="left")  # intervals starting before x
        total = cum[k]
        last = np.clip(k - 1, 0, None)
        overshoot = np.where((k > 0) & (e[last] > x), v[last] * (e[last] - x), 0.0)
        return total - overshoot

    return area_upto(ends) - area_upto(starts)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("regions", choices=["gene_bodies", "promoters"])
    p.add_argument("--workdir", default=WORKDIR)
    p.add_argument("--fraglen", type=int, default=300)
    a = p.parse_args()

    data = os.path.join(a.workdir, "data")
    results = os.path.join(a.workdir, "results")
    samples = pd.read_csv(os.path.join(a.workdir, "metadata", "samples.tsv"), sep="\t")
    regions = pd.read_csv(os.path.join(results, f"{a.regions}.bed"), sep="\t", header=None,
                          names=["chrom", "start", "end", "name", "score", "strand"])

    counts = pd.DataFrame(index=regions["name"])
    totals = {}
    for row in samples.itertuples():
        bg = load_bedgraph(find_bedgraph(data, row.file))
        totals[row.sample] = round(float(((bg["end"] - bg["start"]).astype(np.float64) * bg["cov"]).sum()) / a.fraglen)
        col = np.zeros(len(regions))
        for chrom, reg in regions.groupby("chrom", sort=False):
            sub = bg[bg["chrom"] == chrom]
            if len(sub) == 0:
                continue
            col[reg.index.to_numpy()] = region_area(sub, reg["start"].to_numpy(), reg["end"].to_numpy())
        counts[row.sample] = np.rint(col / a.fraglen).astype(np.int64)
        print(f"{row.sample}: ~{totals[row.sample]:,} reads total, "
              f"{counts[row.sample].sum():,} in regions")
        del bg

    out = os.path.join(results, f"counts_{a.regions}.tsv")
    counts.to_csv(out, sep="\t")
    pd.Series(totals, name="total_reads").to_csv(out.replace(".tsv", "_totals.tsv"), sep="\t")
    print(f"written: {out}")


if __name__ == "__main__":
    main()
