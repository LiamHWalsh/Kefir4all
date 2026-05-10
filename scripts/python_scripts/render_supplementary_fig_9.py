"""Regenerate Supplementary Figure 9: dRep secondary-cluster assignments
across the Kefir4All study, by species and timepoint.

Inputs:
  data/Cdb.csv
  data/mag_metadata_v3.csv.gz

Outputs:
  figures/Supplementary_Fig_9.png
  figures/Supplementary_Fig_9.pdf
  figures/Supplementary_Fig_9_data.tsv

Usage:
  python scripts/python_scripts/render_supplementary_fig_9.py
"""
from __future__ import annotations
import gzip, math, re
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
CDB  = ROOT / "data" /  "Cdb.csv"
MD   = ROOT / "data" /  "mag_metadata_v3.csv.gz"
OUT  = ROOT / "figures"; OUT.mkdir(exist_ok=True)

cdb = pd.read_csv(CDB)
with gzip.open(MD) as f: md = pd.read_csv(f, low_memory=False)


def strip_fa(s): return re.sub(r"\.fa$", "", str(s))
cdb["g"] = cdb["genome"].apply(strip_fa)
md["g"]  = md["user_genome"].apply(strip_fa)

m = cdb.merge(md[["g", "classification", "Sample", "data_source", "kefir type"]], on="g", how="inner")

def stage_from(g):
    mt = re.search(r"_T(\d+)_", str(g))
    return f"T{mt.group(1)}" if mt else None

m["Stage"] = m["genome"].apply(stage_from)
m = m.dropna(subset=["Stage", "classification"])
m = m[m["data_source"] == "This study"]

def short(c):
    mt = re.search(r"s__([^;]+)$", str(c))
    return mt.group(1).strip() if mt else str(c)

m["species"] = m["classification"].apply(short)

def kg(k):
    if k in ("ML", "MG"): return "Milk kefir"
    if k in ("WL", "WG"): return "Water kefir"
    return "Other"

m["kg"] = m["kefir type"].apply(kg)
m = m[m["kg"].isin(["Milk kefir", "Water kefir"])]

def prev(d, n=10):
    c = d["species"].value_counts()
    return c[c >= n].index.tolist()

milk_sp = prev(m[m["kg"] == "Milk kefir"])
water_sp = prev(m[m["kg"] == "Water kefir"])
stage_order = ["T0", "T1", "T2", "T3", "T4", "T5", "T6"]
stage_labels = ["T0", "wk01", "wk05", "wk09", "wk13", "wk17", "wk21"]

def panel(ax, sub, title):
    if sub.empty: ax.set_visible(False); return
    pivot = (sub.groupby(["Stage", "secondary_cluster"]).size()
                .unstack(fill_value=0).reindex(stage_order, fill_value=0))
    bottoms = np.zeros(len(stage_order))
    cmap = plt.colormaps.get_cmap("tab20")
    for i, c in enumerate(pivot.columns):
        vals = pivot[c].values
        ax.bar(stage_labels, vals, bottom=bottoms, color=cmap(i % 20), edgecolor="white", linewidth=0.3)
        bottoms += vals
    ax.set_title(title, fontsize=8)
    ax.tick_params(axis="x", rotation=45, labelsize=6)
    ax.tick_params(axis="y", labelsize=6)
    ax.set_ylabel("MAGs", fontsize=6)

all_sp = [(s, "milk") for s in milk_sp] + [(s, "water") for s in water_sp]
cols = 4; rows = math.ceil(len(all_sp) / cols)
fig, axes = plt.subplots(rows, cols, figsize=(cols * 3, rows * 2.0))
flat = axes.ravel() if rows * cols > 1 else [axes]
milk = m[m["kg"] == "Milk kefir"]
water = m[m["kg"] == "Water kefir"]
for ax, (sp, kgrp) in zip(flat, all_sp):
    sub = milk[milk["species"] == sp] if kgrp == "milk" else water[water["species"] == sp]
    panel(ax, sub, f"{sp} ({'milk' if kgrp == 'milk' else 'water'} kefir)")
for ax in flat[len(all_sp):]:
    ax.set_visible(False)

plt.suptitle("Supplementary Fig. 9. dRep secondary-cluster assignments across the Kefir4All study, by species and timepoint.", fontsize=10)
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(OUT / "Supplementary_Fig_9.png", dpi=200, bbox_inches="tight")
plt.savefig(OUT / "Supplementary_Fig_9.pdf", bbox_inches="tight")

(m[["species", "kg", "Stage", "secondary_cluster", "primary_cluster", "Sample"]]
 .rename(columns={"kg": "kefir_type", "Stage": "stage"})
 .to_csv(OUT / "Supplementary_Fig_9_data.tsv", sep="\t", index=False))
print(f"Wrote {OUT / 'Supplementary_Fig_9.png'}")
