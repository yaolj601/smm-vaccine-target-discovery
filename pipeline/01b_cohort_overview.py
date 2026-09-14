"""
Step 1b: Cohort overview figure.

A Figure-1-style summary of what the cohort actually contains, per sample:
disease stage, SMM risk stratum, how many plasma cells were profiled, and how
those cells split into neoplastic vs normal.

Note on longitudinal data: GSE193531 has NONE. All 35 barcode suffixes map 1:1
to the 35 sample IDs, there are no repeat runs, and every library is a single
CD138+ fraction (".138P"). So a swimmer/timeline layout is not possible - this
is the cohort-structure equivalent, and the absence of serial sampling is
itself worth showing.

Input:  $SMM_DATA_DIR/smm_plasma_cells.h5ad
Output: $SMM_OUT_DIR/figures/cohort_overview.png|.svg
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

from config import DATA, OUT

FIG = f"{OUT}/figures"
os.makedirs(FIG, exist_ok=True)
matplotlib.rcParams["font.family"] = ["Liberation Sans", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"

STAGE_ORDER = ["NBM", "MGUS", "SMM", "MM"]
C_NEO = "#9B2C62"      # neoplastic
C_NORM = "#3E7C77"     # normal
C_STAGE = {"NBM": "#8C7480", "MGUS": "#C79A3E", "SMM": "#9B2C62", "MM": "#5B2E7A"}
C_RISK = {"high": "#A02C2C", "low": "#C7A89B", "n/a": "#E3D8D5"}

print("Loading AnnData...")
adata = sc.read_h5ad(f"{DATA}/smm_plasma_cells.h5ad")
m = adata.obs

tab = pd.crosstab(m["sample_ID"], m["normal_or_neoplastic"])
for c in ("neoplastic", "normal"):
    if c not in tab:
        tab[c] = 0
tab["total"] = tab["neoplastic"] + tab["normal"]
tab["stage"] = m.groupby("sample_ID", observed=True)["disease_stage"].first()
risk = (m.groupby("sample_ID", observed=True)["smm_risk"].first()
        if "smm_risk" in m.columns else pd.Series("n/a", index=tab.index))
tab["risk"] = risk
tab["paired"] = (tab["neoplastic"] > 0) & (tab["normal"] > 0)

# order: stage, then high-risk first within SMM, then by tumour burden
tab["stage_i"] = tab["stage"].astype(str).map(
    {s: i for i, s in enumerate(STAGE_ORDER)}).astype(float)
tab["risk_i"] = tab["risk"].astype(str).map(
    {"high": 0, "low": 1, "n/a": 2}).fillna(2).astype(float)
tab = tab.sort_values(["stage_i", "risk_i", "neoplastic"], ascending=[True, True, False])

n = len(tab)
y = np.arange(n)[::-1]

fig = plt.figure(figsize=(11.5, 0.30 * n + 2.0))
gs = fig.add_gridspec(1, 3, width_ratios=[0.5, 3.1, 1.5], wspace=0.06)
ax_t = fig.add_subplot(gs[0])   # stage + risk tracks
ax_b = fig.add_subplot(gs[1])   # stacked composition
ax_p = fig.add_subplot(gs[2])   # paired / tumour purity

# ---------------------------------------------------------- annotation tracks
for i, (sid, r) in zip(y, tab.iterrows()):
    ax_t.add_patch(Rectangle((0, i - .38), .46, .76,
                             facecolor=C_STAGE[r["stage"]], edgecolor="none"))
    ax_t.add_patch(Rectangle((.52, i - .38), .46, .76,
                             facecolor=C_RISK.get(r["risk"], "#E3D8D5"), edgecolor="none"))
ax_t.set_xlim(-.06, 1.04); ax_t.set_ylim(-.8, n - .2)
ax_t.set_yticks(y); ax_t.set_yticklabels(tab.index, fontsize=8)
ax_t.set_xticks([.23, .75]); ax_t.set_xticklabels(["stage", "risk"], fontsize=8, rotation=90)
for sp in ax_t.spines.values():
    sp.set_visible(False)
ax_t.tick_params(length=0)

# ------------------------------------------------------- cells profiled (log)
ax_b.barh(y, tab["neoplastic"], color=C_NEO, height=.72, label="neoplastic PC")
ax_b.barh(y, tab["normal"], left=tab["neoplastic"], color=C_NORM, height=.72,
          label="normal PC")
for i, (sid, r) in zip(y, tab.iterrows()):
    ax_b.text(r["total"] + 60, i, f"{int(r['total']):,}", va="center",
              fontsize=7.5, color="#6B5560")
ax_b.set_xlim(0, tab["total"].max() * 1.16)
ax_b.set_ylim(-.8, n - .2)
ax_b.set_yticks([])
ax_b.set_xlabel("CD138+ plasma cells profiled", fontsize=9)
ax_b.legend(loc="upper right", fontsize=8, frameon=False,
            bbox_to_anchor=(1.0, 1.0))
for sp in ("top", "right", "left"):
    ax_b.spines[sp].set_visible(False)
ax_b.grid(axis="x", color="#E3D8D5", lw=.6, zorder=0)
ax_b.set_axisbelow(True)

# ------------------------------------------------- tumour fraction + pairing
frac = np.where(tab["total"] > 0, tab["neoplastic"] / tab["total"].clip(lower=1), 0)
ax_p.barh(y, frac, color="#C9B3BE", height=.5, zorder=1)
ax_p.scatter(frac, y, s=26, color=C_NEO, zorder=3)
for i, p in zip(y[tab["paired"].values], frac[tab["paired"].values]):
    ax_p.scatter(1.14, i, marker="D", s=22, color="#2F6B67", zorder=3)
ax_p.set_xlim(-.03, 1.24); ax_p.set_ylim(-.8, n - .2)
ax_p.set_yticks([])
ax_p.set_xticks([0, .5, 1]); ax_p.set_xticklabels(["0", "50%", "100%"], fontsize=8)
ax_p.set_xlabel("tumour fraction of PCs", fontsize=9)
for sp in ("top", "right", "left"):
    ax_p.spines[sp].set_visible(False)
ax_p.grid(axis="x", color="#E3D8D5", lw=.6, zorder=0)
ax_p.set_axisbelow(True)

handles = ([Patch(facecolor=C_STAGE[s], label=s) for s in STAGE_ORDER]
           + [Patch(facecolor=C_RISK["high"], label="SMM high-risk"),
              Patch(facecolor=C_RISK["low"], label="SMM low-risk"),
              Line2D([], [], marker="D", ls="none", color="#2F6B67",
                     label="has BOTH tumour + normal PCs (paired)")])
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .045),
           ncol=4, fontsize=8, frameon=False)

n_pair = int(tab["paired"].sum())
fig.suptitle(f"GSE193531 cohort: {n} patients, one CD138+ aspirate each "
             "\u2014 no longitudinal sampling", fontsize=12.5, y=.995)
fig.text(.5, .972, f"{n_pair} samples carry both neoplastic and normal plasma cells "
                   f"(the only within-patient comparisons available)",
         ha="center", fontsize=9, color="#6B5560")

fig.subplots_adjust(left=.085, right=.985, top=.955, bottom=.085)
plt.savefig(f"{FIG}/cohort_overview.png", dpi=200)
plt.savefig(f"{FIG}/cohort_overview.svg")
plt.close()

print(f"\nsamples: {n}   paired (tumour+normal): {n_pair}")
print(tab.groupby("stage", observed=True)
      .agg(samples=("total", "size"), cells=("total", "sum"),
           paired=("paired", "sum")).to_string())
print(f"\nSaved {FIG}/cohort_overview.png|.svg")
