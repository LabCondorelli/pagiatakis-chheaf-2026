#!/usr/bin/env python3
"""
Fig. 1E, coverage tracks - strand resolved RNA-seq coverage at the syntenic locus, averaged
within each group (non-failing and hypertrophic cardiomyopathy).

The merged bigWig tracks are already normalised across groups, so the two are comparable in
height. The Uniq versions are used, matching the featureCounts quantification, which counts
uniquely mapping reads.

Usage: python3 fig1E_coverage_tracks.py
Output: results/coverage_{hcm,nf}_{plus,minus}.tsv
"""
import math, os
import pyBigWig

# Directory with the merged bigWigs: {group}_merged_{strand}Uniq.bw
WIGS = os.environ.get("CHH_WIGS", "<PATH_TO>/wigs")
OUT  = os.environ.get("CHH_OUT", "results")
os.makedirs(OUT, exist_ok=True)

CHROM, START, END = "chr2", 55105000, 55145000   # displayed window, with margin

for group in ("hcm", "nf"):
    for strand in ("plus", "minus"):
        src = f"{WIGS}/{group}_merged_{strand}Uniq.bw"
        if not os.path.exists(src):
            print(f"  missing: {src}")
            continue
        bw = pyBigWig.open(src)
        vals = bw.values(CHROM, START, END, numpy=False)
        bw.close()
        cov = [0.0 if (v is None or math.isnan(v)) else float(v) for v in vals]
        dst = f"{OUT}/coverage_{group}_{strand}.tsv"
        with open(dst, "w") as fh:
            for i, v in enumerate(cov):
                fh.write(f"{START + i}\t{v:.6f}\t{strand}\n")
        print(f"  {os.path.basename(dst):<28} {len(cov):>6} bases | non-zero "
              f"{sum(1 for v in cov if v > 0):>6} | max {max(cov, default=0.0):.4f}")
