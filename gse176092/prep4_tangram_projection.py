#!/usr/bin/env python3
"""
Tangram projection (mode='clusters') of the 19 genes of Fig. 1D and SX: published map x summed
counts per cell type. --map none fits Tangram again (seed 42; close to, not the published one).

Usage:
  python3 prep4_tangram_projection.py --dir work/tangram --map tangram_map_manuscript.tsv \
      --counts tangram_celltype_counts_manuscript.tsv
Output: <dir>/tangram_projection_19genes.csv, <dir>/tangram_celltype_abundance.csv
"""
import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc

GENES = ["gm12092", "rtn4", "gm13481", "zeb2", "acvr2a", "gm5086", "enc1",
         "gm30085", "pdzrn4", "gxylt1", "a730036i17rik", "slc20a1", "ckap2l",
         "gm17501", "ccn1", "ddah1", "gm42047", "klhdc4", "slc7a5"]
SEED = 42

ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True, help="folder with spatial.h5ad (and scrna.h5ad with --map none)")
ap.add_argument("--map", default="none", help="published cell-type map, spots x cell types; 'none' fits Tangram")
ap.add_argument("--counts", help="published summed counts per cell type, genes x cell types")
a = ap.parse_args()
d = a.dir

if a.map != "none":
    tmap = pd.read_csv(a.map, sep="\t", index_col=0).astype(np.float32)
    counts = pd.read_csv(a.counts, sep="\t", index_col=0)[tmap.columns]
    spots = sc.read_h5ad(os.path.join(d, "spatial.h5ad"), backed="r").obs_names
    if list(spots) != list(tmap.index):
        raise SystemExit("the spots of the map are not those of spatial.h5ad")
    missing = [g for g in GENES if g not in counts.index]
    if missing:
        raise SystemExit(f"genes missing from the counts: {missing}")
    out = pd.DataFrame(tmap.values @ counts.loc[GENES].values.T, index=tmap.index, columns=GENES)
    abundance = tmap
else:
    import tangram as tg
    ad_sc = sc.read_h5ad(os.path.join(d, "scrna.h5ad"))
    ad_sp = sc.read_h5ad(os.path.join(d, "spatial.h5ad"))
    tg.pp_adatas(ad_sc, ad_sp, genes=None)
    ad_map = tg.map_cells_to_space(ad_sc, ad_sp, mode="clusters", cluster_label="cell.type",
                                   random_state=SEED)
    proj = tg.project_genes(ad_map, ad_sc, cluster_label="cell.type").to_df()
    missing = [g for g in GENES if g not in proj.columns]
    if missing:
        raise SystemExit(f"genes missing from the projection: {missing}")
    out = proj[GENES]
    abundance = pd.DataFrame(ad_map.X.T, index=ad_map.var_names, columns=ad_map.obs["cell.type"])

out.index.name = "spot_id"
out.to_csv(os.path.join(d, "tangram_projection_19genes.csv"))
abundance.to_csv(os.path.join(d, "tangram_celltype_abundance.csv"))
print(f"projection ({'published map' if a.map != 'none' else 'Tangram fitted, seed 42'}): "
      f"{out.shape[0]} spots x {out.shape[1]} genes -> {d}/tangram_projection_19genes.csv")
print(out.describe().loc[["mean", "max"]].T.round(4).to_string())
