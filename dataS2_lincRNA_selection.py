#!/usr/bin/env python3
"""
Data S2 - the 40 lincRNA modulated by TAC, their flanking coding genes and the
cis-candidate selection rule.

Rule: a modulated lincRNA is a cis-regulatory candidate if it is up-regulated and at least
one of its two flanking coding genes is itself modulated in at least one TAC vs sham
comparison. Direction comes from the sign of the differential calls, not from clustering.
Trajectory concordance is recorded but does not enter the selection.

Flanking gene status has three values, and they are not interchangeable: 'modulated',
'tested, not modulated', and 'not expressed' (below the edgeR expression filter, therefore
never tested).

Usage: python3 dataS2_lincRNA_selection.py   (run dataS2_00_neighbours.sh first)
Output: results/dataS2_lincRNA_selection.tsv, results/dataS2_methods.tsv
"""
import csv, os
import pandas as pd

# Differential expression table of the mouse TAC time course (limma-voom, GENCODE vM25), with
# every gene that passed the expression filter, not only the significant ones: a flanking gene
# absent from it was never tested, which is not the same as tested and not modulated.
# 00_de_and_cpm.R writes it as <prefix>_DE_calls.tsv.
DE  = os.environ.get("CHH_TAC_DE_ALL", "<PATH_TO>/mouse_tac_DE_calls.tsv")
OUT = os.environ.get("CHH_OUT", "results")
RES=["Results.tac_1d_vs_sham_1d","Results.tac_2d_vs_sham_2d",
     "Results.tac_4d_vs_sham_4d","Results.tac_7d_vs_sham_7d"]
TP=["1d","2d","4d","7d"]

def read_de(fn):
    """Every required column has to be there. A missing one used to read as a call of zero,
    which turns a renamed contrast into a table of genes that look not modulated."""
    out={}
    with open(fn) as fh:
        r=csv.reader(fh,delimiter="\t"); hdr=next(r); off=None
        need=set(RES)|{"Genes.ensembl_gene_id","Genes.gene_biotype","Genes.external_gene_name"}
        absent=need-set(hdr)
        if absent: raise SystemExit(f"{fn}: columns absent: {', '.join(sorted(absent))}")
        for row in r:
            if off is None: off=len(row)-len(hdr)
            d=dict(zip(hdr,row[off:] if off>0 else row))
            gid=d.get("Genes.ensembl_gene_id","")
            if not gid: continue
            if gid in out: raise SystemExit(f"{fn}: gene {gid} appears twice")
            out[gid]=d
    return out
allg=read_de(DE)                                   # every gene that passed the filter

def mods(d):
    v=[]
    for k in RES:
        c=int(float(d[k]))
        if c not in (-1,0,1): raise SystemExit(f"{k}: call {c}, expected -1, 0 or 1")
        v.append(c)
    return v
sign={g: d for g, d in allg.items() if any(mods(d))}   # the significant ones

def is_mod(gid, src):
    d=src.get(gid); return bool(d) and any(x!=0 for x in mods(d))
def mod_status(gid, src):
    """Three states. A neighbour missing from the DE table is not 'not modulated':
    it failed filterByExpr and was never tested."""
    if not gid: return "no neighbour annotated"
    d=src.get(gid)
    if d is None: return "not expressed"
    return "modulated" if any(x!=0 for x in mods(d)) else "tested, not modulated"

# --- flanking genes ---
def read_near(f, col_prefix):
    out={}
    for line in open(f):
        p=line.rstrip("\n").split("\t")
        nc=p[3].split(".")[0]
        if p[11]==".": out[nc]=(None,None,None,None,None); continue
        # Column order: p[13] is the strand, p[14] the gene name (from GENCODE, so present
        # also for genes absent from the DE table).
        # Distance between the two TSS, taken according to strand, 1-based, negative on the left.
        tss_a = int(p[1]) + 1 if p[5] == "+" else int(p[2])
        tss_b = int(p[9]) + 1 if p[13] == "+" else int(p[10])
        out[nc]=(p[11], p[11].split(".")[0], p[14], int(p[16]), tss_b - tss_a)   # id_ver, id, name, gap, TSS
    return out
# -io excluded overlapping genes, so read-through models do not hide the real neighbours.
up=read_near(f"{OUT}/neighbours_left.txt","up"); dn=read_near(f"{OUT}/neighbours_right.txt","dn")

# --- the 40 modulated lincRNA ---
linc=[g for g,d in sign.items() if d.get("Genes.gene_biotype")=="lincRNA" and any(x!=0 for x in mods(d))]
print(f"  modulated lincRNA: {len(linc)}  (expected 40)")

# distances in kb, one decimal, half-up rounding
from decimal import Decimal, ROUND_HALF_UP
def kb(bp):
    if bp is None: return ""
    return float((Decimal(bp) / 1000).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))

rows=[]
for g in linc:
    d=sign[g]; m=mods(d)
    u=up.get(g,(None,)*5); w=dn.get(g,(None,)*5)
    u_mod = is_mod(u[1], allg) if u[1] else False
    w_mod = is_mod(w[1], allg) if w[1] else False
    passes = u_mod or w_mod
    # the neighbour name is always reported, from GENCODE when the gene is not in the DE table
    def nm(x, fallback=""):
        v = allg.get(x,{}).get("Genes.external_gene_name","") if x else ""
        return v or (fallback or "")
    def mv(x): return "/".join(str(v) for v in mods(allg[x])) if x and x in allg else ""
    rows.append({
        "lincRNA_ensembl":g, "lincRNA_name":d.get("Genes.external_gene_name",""),
        **{f"lincRNA_{t}":m[i] for i,t in enumerate(TP)},
        "left_PC_ensembl":u[0] or "", "left_PC_name":nm(u[1], u[2]),
        "left_TSS_distance_kb":kb(u[4]),
        "left_PC_status":mod_status(u[1], allg), "left_PC_results_1d/2d/4d/7d":mv(u[1]),
        "right_PC_ensembl":w[0] or "", "right_PC_name":nm(w[1], w[2]),
        "right_TSS_distance_kb":kb(w[4]),
        "right_PC_status":mod_status(w[1], allg), "right_PC_results_1d/2d/4d/7d":mv(w[1]),
        "selected_as_cis_candidate":"YES" if passes else "no",
        "reason": ("at least one flanking coding gene modulated" if passes
                   else "neither flanking coding gene modulated"),
    })
df=pd.DataFrame(rows)

# --- direction, from the sign of the differential calls -----------------------
# Up-regulated when more time points call the lincRNA up than down: the direction comes from
# the calls themselves, so the assignment is deterministic.
import numpy as np
call = df[[f"lincRNA_{t}" for t in TP]].to_numpy()
n_up, n_dn = (call == 1).sum(1), (call == -1).sum(1)
df["direction"] = np.where(n_up > n_dn, "up-regulated", "down-regulated")
print(f"  direction: {(df.direction=='up-regulated').sum()} up-regulated, "
      f"{(df.direction=='down-regulated').sum()} down-regulated")

# --- selection rule -----------------------------------------------------------
mod = lambda c: df[c].eq("modulated")
df["flanking_PC_modulated"] = np.where(mod("left_PC_status") | mod("right_PC_status"), "yes", "no")
df["cis_candidate"] = np.where((df.direction=="up-regulated") & (df.flanking_PC_modulated=="yes"), "YES", "no")
df["outcome"] = np.where(df.direction!="up-regulated", "not up-regulated",
                  np.where(df.flanking_PC_modulated=="no",
                           "up-regulated, but neither flanking coding gene modulated",
                           "up-regulated and at least one flanking coding gene modulated"))
n=(df.cis_candidate=="YES").sum()
print(f"  cis-regulatory candidates: {n}")
for k,v in df.outcome.value_counts().items(): print(f"    {k:<52} {v}")

cols=["lincRNA_name","lincRNA_ensembl"]+[f"lincRNA_{t}" for t in TP]+[
 "direction","left_PC_name","left_TSS_distance_kb","left_PC_status",
 "left_PC_results_1d/2d/4d/7d","right_PC_name","right_TSS_distance_kb",
 "right_PC_status","right_PC_results_1d/2d/4d/7d",
 "flanking_PC_modulated","cis_candidate","outcome"]
# Order: candidates first, then the other up-regulated, then the down-regulated, by name
# within each block. Sorting on explicit booleans, because sorting the strings would not work.
out=df[cols].assign(_c=df.cis_candidate.eq("YES"),
                    _d=df.direction.eq("up-regulated")).sort_values(
        ["_c","_d","lincRNA_name"], ascending=[False,False,True]).drop(columns=["_c","_d"])
out.to_csv(f"{OUT}/dataS2_lincRNA_selection.tsv",sep="\t",index=False)
# The method summary that ships with the table.
meta=pd.DataFrame({"field":["Dataset","DE variant","Significance","lincRNA set",
                            "Neighbours",
                            "Selection rule","Result","Reproducibility note"],
 "value":["NCBI SRA PRJNA1372455: RNA-seq of isolated cardiomyocytes at 1, 2, 4 and 7 days "
          "after TAC or sham surgery",
          "limma-voom, TAC versus sham at each time point, on the GENCODE vM25 annotation "
          "(all biotypes)",
          "FDR < 0.05 and |log2FC| > 0.5 (Results columns, coded -1/0/+1)",
          "40 lincRNAs modulated in at least one post-TAC time point, of which 15 up-regulated "
          "and 25 down-regulated",
          "bedtools closest v2.31.1 on GENCODE vM25: -D ref, -id / -iu, -io (overlaps excluded), -t first",
          "retained if up-regulated AND at least one flanking coding gene modulated; "
          "trajectory concordance recorded but not used for selection",
          "7 candidate cis-regulatory lincRNAs, including Gm12092 (Chheaf)",
          "flanking genes are reported in one of three states. 'modulated': significant at one or "
          "more time points. 'tested, not modulated': expressed and tested, but not significant at "
          "any time point. 'not expressed': below the edgeR filterByExpr expression filter and "
          "therefore never tested, which is not the same as tested and found unchanged"]})
meta.to_csv(f"{OUT}/dataS2_methods.tsv",sep="\t",index=False)
print(f"\n  written: {OUT}/dataS2_lincRNA_selection.tsv and dataS2_methods.tsv")
print("\n  candidates:")
print(out[out.cis_candidate=="YES"][["lincRNA_name","left_PC_status","right_PC_status"]].to_string(index=False))
