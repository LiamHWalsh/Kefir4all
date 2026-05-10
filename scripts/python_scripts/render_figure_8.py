"""Regenerate main-text Figure 8 (revised manuscript): inStrain-defined
within-species cluster detections by species and timepoint, milk and water
kefir prevalent species. No cropping; built directly from the inStrain
genome-level summary used in the strain analysis.

Inputs:
  data/instrain_genome_species_primary_data_v4.csv

Outputs:
  figures/Figure_8.png
  figures/Figure_8.pdf
  figures/Figure_8_data.tsv

Usage:
  python scripts/python_scripts/render_figure_8.py
"""
from __future__ import annotations
import math, re
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC  = ROOT / "data" /  "instrain_genome_species_primary_data_v4.csv"
OUT  = ROOT / "figures"
OUT.mkdir(exist_ok=True)

df = pd.read_csv(SRC, low_memory=False)
df = df.dropna(subset=["breadth", "breadth_expected", "Stage"])
df = df[(df["breadth"] >= 0.35) & ((df["breadth"] / df["breadth_expected"]) >= 0.75)]
df = df[df["Stage"].astype(str).str.match(r"^T\d$")]
if "data_source.x" in df.columns:
    df = df[df["data_source.x"] == "This study"]


def kefir_from_genome(g: str) -> str | None:
    m = re.search(r"_(W[GL]|M[GL])_", str(g))
    if not m:
        return None
    return "Milk kefir" if m.group(1).startswith("M") else "Water kefir"


df["kg"] = df["genome"].apply(kefir_from_genome)
df = df.dropna(subset=["kg"])

stage_order = ["T0", "T1", "T2", "T3", "T4", "T5", "T6"]
stage_labels = ["T0", "wk01", "wk05", "wk09", "wk13", "wk17", "wk21"]


def prev(d, n=10):
    c = d["classification"].value_counts()
    return c[c >= n].index.tolist()


milk_sp = prev(df[df["kg"] == "Milk kefir"])
water_sp = prev(df[df["kg"] == "Water kefir"])


def panel(ax, sub, title):
    if sub.empty:
        ax.set_visible(False)
        return
    pivot = (
        sub.groupby(["Stage", "cluster"]).size().unstack(fill_value=0).reindex(stage_order, fill_value=0)
    )
    bottoms = np.zeros(len(stage_order))
    cmap = plt.colormaps.get_cmap("tab20")
    for i, c in enumerate(pivot.columns):
        vals = pivot[c].values
        ax.bar(stage_labels, vals, bottom=bottoms, color=cmap(i % 20), edgecolor="white", linewidth=0.3)
        bottoms += vals
    ax.set_title(title, fontsize=8)
    ax.tick_params(axis="x", rotation=45, labelsize=6)
    ax.tick_params(axis="y", labelsize=6)
    ax.set_ylabel("metagenomes", fontsize=6)


all_sp = [(s, "milk") for s in milk_sp] + [(s, "water") for s in water_sp]
n = len(all_sp)
cols = 4
rows = math.ceil(n / cols)
fig, axes = plt.subplots(rows, cols, figsize=(cols * 3, rows * 2.0))
flat = axes.ravel() if rows * cols > 1 else [axes]
milk = df[df["kg"] == "Milk kefir"]
water = df[df["kg"] == "Water kefir"]
for ax, (sp, kg) in zip(flat, all_sp):
    sub = milk[milk["classification"] == sp] if kg == "milk" else water[water["classification"] == sp]
    panel(ax, sub, f"{sp} ({'milk' if kg == 'milk' else 'water'} kefir)")
for ax in flat[len(all_sp) :]:
    ax.set_visible(False)

plt.suptitle(
    "Figure 8. inStrain-defined within-species cluster detections across the Kefir4All study, by species and timepoint.",
    fontsize=10,
)
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(OUT / "Figure_8.png", dpi=220, bbox_inches="tight")
plt.savefig(OUT / "Figure_8.pdf", bbox_inches="tight")

(df[["classification", "kg", "Stage", "cluster", "sample_id"]]
   .rename(columns={"kg": "kefir_type", "Stage": "stage"})
   .to_csv(OUT / "Figure_8_data.tsv", sep="\t", index=False))
print(f"Wrote {OUT / 'Figure_8.png'}")
