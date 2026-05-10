"""Regenerate Supplementary Figure 8: within-sample polymorphism (per-genome
SNV count and nucleotide diversity) across prevalent kefir species.

Inputs:
  data/instrain_genome_species_primary_data_v4.csv

Outputs:
  figures/Supplementary_Fig_8.png
  figures/Supplementary_Fig_8.pdf
  figures/Supplementary_Fig_8_data.tsv

Usage:
  python scripts/python_scripts/render_supplementary_fig_8.py
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC  = ROOT / "data" /  "instrain_genome_species_primary_data_v4.csv"
OUT  = ROOT / "figures"; OUT.mkdir(exist_ok=True)

df = pd.read_csv(SRC, low_memory=False)
df = df.dropna(subset=["breadth", "breadth_expected", "SNV_count", "nucl_diversity"])
df = df[(df["breadth"] >= 0.35) & ((df["breadth"] / df["breadth_expected"]) >= 0.75)]

milk = df[df["kefir type.x"].isin(["ML", "MG"])]
water = df[df["kefir type.x"].isin(["WL", "WG"])]

def filter_prev(d, n=10):
    c = d["classification"].value_counts()
    return d[d["classification"].isin(c[c >= n].index)]

mp = filter_prev(milk); wp = filter_prev(water)

def order(d):
    return d.groupby("classification")["SNV_count"].median().sort_values(ascending=False).index.tolist()

mo, wo = order(mp), order(wp)

fig, axes = plt.subplots(2, 2, figsize=(11, 10))
def boxp(ax, d, oi, ycol, ylabel, title, log=False):
    groups = [d[d["classification"] == s][ycol].values for s in oi]
    bp = ax.boxplot(groups, tick_labels=oi, vert=True, showfliers=False, patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("#4C9F70"); patch.set_alpha(0.7)
    if log: ax.set_yscale("log")
    ax.set_ylabel(ylabel); ax.set_title(title)
    ax.tick_params(axis="x", rotation=70, labelsize=7)
    for label in ax.get_xticklabels():
        label.set_horizontalalignment("right")

boxp(axes[0,0], mp, mo, "SNV_count", "SNV sites per genome (log)", "(A) Milk kefir – per-genome SNV counts", log=True)
boxp(axes[0,1], wp, wo, "SNV_count", "SNV sites per genome (log)", "(B) Water kefir – per-genome SNV counts", log=True)
boxp(axes[1,0], mp, mo, "nucl_diversity", "Nucleotide diversity",       "(C) Milk kefir – nucleotide diversity")
boxp(axes[1,1], wp, wo, "nucl_diversity", "Nucleotide diversity",       "(D) Water kefir – nucleotide diversity")

plt.suptitle(
    "Supplementary Fig. 8. Within-sample polymorphism (per-genome SNV sites and nucleotide diversity) across prevalent kefir species, indexed by inStrain.",
    fontsize=10,
)
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(OUT / "Supplementary_Fig_8.png", dpi=200, bbox_inches="tight")
plt.savefig(OUT / "Supplementary_Fig_8.pdf", bbox_inches="tight")

# Per-species summary
import pandas as pd
summ = (df[df["classification"].isin(set(mo) | set(wo))]
        .groupby(["kefir type.x", "classification"]).agg(
    n_detections=("SNV_count", "size"),
    median_snv=("SNV_count", "median"),
    median_nucl_div=("nucl_diversity", "median"),
).reset_index())
summ.to_csv(OUT / "Supplementary_Fig_8_data.tsv", sep="\t", index=False, float_format="%.4f")
print(f"Wrote {OUT / 'Supplementary_Fig_8.png'}")
