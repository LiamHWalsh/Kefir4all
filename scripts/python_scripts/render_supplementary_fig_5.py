"""Regenerate Supplementary Figure 5 (revised manuscript): per-species
proportion of within-species MAG pairs reaching ≥99% ANI.

Inputs:
  data/ndb.csv.gz       (dRep ANI table)
  data/Cdb.csv           (dRep cluster assignments)
  data/mag_metadata_v3.csv.gz

Outputs:
  figures/Supplementary_Fig_5.png
  figures/Supplementary_Fig_5.pdf
  figures/Supplementary_Fig_5_data.tsv

Usage:
  python scripts/python_scripts/render_supplementary_fig_5.py
"""
from __future__ import annotations
import gzip, re
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
NDB  = ROOT / "data" /  "ndb.csv.gz"
CDB  = ROOT / "data" /  "Cdb.csv"
MD   = ROOT / "data" /  "mag_metadata_v3.csv.gz"
OUT  = ROOT / "figures"; OUT.mkdir(exist_ok=True)

with gzip.open(NDB) as f: ndb = pd.read_csv(f)
cdb = pd.read_csv(CDB)
with gzip.open(MD) as f: md = pd.read_csv(f, low_memory=False)

def strip_fa(s): return re.sub(r"\.fa$", "", str(s))
ndb["q"] = ndb["querry"].apply(strip_fa)
ndb["r"] = ndb["reference"].apply(strip_fa)
cdb["g"] = cdb["genome"].apply(strip_fa)
md["g"]  = md["user_genome"].apply(strip_fa)

g_pc = dict(zip(cdb["g"], cdb["primary_cluster"]))
g_sp = dict(zip(md["g"], md["classification"]))
this_study = set(md[md["data_source"] == "This study"]["g"]) if "data_source" in md.columns else set(md["g"])

ndb["q_pc"] = ndb["q"].map(g_pc)
ndb["r_pc"] = ndb["r"].map(g_pc)
ndb["q_sp"] = ndb["q"].map(g_sp)
within = ndb[(ndb["q_pc"] == ndb["r_pc"]) & ndb["q_pc"].notna() & (ndb["q"] != ndb["r"])]
within = within[within["q"].isin(this_study) & within["r"].isin(this_study)]


def short(c):
    m = re.search(r"s__([^;]+)$", str(c))
    return m.group(1).strip() if m else str(c)


within = within.assign(species=within["q_sp"].apply(short))
agg = within.groupby("species").agg(n_pairs=("ani", "size"), n_ge99=("ani", lambda s: (s >= 0.99).sum())).reset_index()
agg["pct_ge99"] = 100 * agg["n_ge99"] / agg["n_pairs"]
agg = agg[agg["n_pairs"] >= 3].sort_values("pct_ge99", ascending=False)
agg.to_csv(OUT / "Supplementary_Fig_5_data.tsv", sep="\t", index=False, float_format="%.2f")

top = agg.head(30).sort_values("pct_ge99", ascending=True)
fig, ax = plt.subplots(figsize=(8, 9))
bars = ax.barh(top["species"], top["pct_ge99"], color="#4C9F70", edgecolor="white")
ax.set_xlabel("Pairwise comparisons within species ≥99% ANI (%)")
ax.set_xlim(0, 100)
for bar, n in zip(bars, top["n_pairs"]):
    ax.text(bar.get_width() + 1, bar.get_y() + bar.get_height() / 2, f"n={n}", va="center", fontsize=7)
ax.set_title("Supplementary Fig. 5. Per-species proportion of within-species MAG pairs reaching ≥99% ANI.", fontsize=10)
plt.tight_layout()
plt.savefig(OUT / "Supplementary_Fig_5.png", dpi=200, bbox_inches="tight")
plt.savefig(OUT / "Supplementary_Fig_5.pdf", bbox_inches="tight")
print(f"Wrote {OUT / 'Supplementary_Fig_5.png'}")
