"""
Step 1: Preprocess GSE193531.

Loads the UMI count matrix and cell metadata, aligns them, applies light QC,
and saves a single AnnData file that all later steps read.

Input:  $SMM_DATA_DIR/counts.csv.gz  (GSE193531_umi-count-matrix.csv.gz)
        $SMM_DATA_DIR/meta.csv.gz    (GSE193531_cell-level-metadata.csv.gz)
Output: $SMM_DATA_DIR/smm_plasma_cells.h5ad

Run 00_fetch_data.py first if the inputs are not present.
"""

import numpy as np
import pandas as pd
import scanpy as sc

from config import COUNTS_CSV, DATA, META_CSV, is_excluded

# pandas 3.x stores strings as Arrow arrays, which AnnData cannot write to
# h5ad. Turn that off so all string columns stay plain object dtype.
pd.options.future.infer_string = False

# ---------------------------------------------------------------- load data
for _p in (COUNTS_CSV, META_CSV):
    if not _p.exists():
        raise SystemExit(f"missing input: {_p}\nRun 00_fetch_data.py first.")

print("Loading count matrix (this is the slow step, ~1-2 min)...")
counts = pd.read_csv(COUNTS_CSV, index_col=0)
meta = pd.read_csv(META_CSV, index_col=0)

print(f"counts: {counts.shape[0]} genes x {counts.shape[1]} cells")
print(f"meta:   {meta.shape[0]} cells")

# cells are rows in meta, columns in counts -> align them
shared = meta.index.intersection(counts.columns)
counts = counts[shared]
meta = meta.loc[shared]
print(f"aligned cells: {len(shared)}")

# ------------------------------------------------------------- clean labels
# 20 cells have a blank normal_or_neoplastic label -> drop them
meta = meta[meta["normal_or_neoplastic"].isin(["normal", "neoplastic"])].copy()
counts = counts[meta.index]
print(f"after dropping blank labels: {counts.shape[1]} cells")

# ------------------------------------------------- recover SMM risk stratum
# The GEO `sample_ID` column flattens SMM-1..12, but the cell barcodes preserve
# the label the authors actually used: SMMh (high-risk) vs SMMl (low-risk).
# That distinction matters here - high-risk SMM (~50% progression at 2 years)
# is the population a preventive vaccine would enrol, low-risk is far more
# indolent - so carry it through rather than let sample_ID discard it.
meta["orig_label"] = [b.split("-", 2)[2].split(".")[0] if b.count("-") >= 2 else ""
                      for b in meta.index]
meta["smm_risk"] = np.where(meta["orig_label"].str.startswith("SMMh"), "high",
                   np.where(meta["orig_label"].str.startswith("SMMl"), "low", "n/a"))
_r = meta.groupby("smm_risk", observed=True)["sample_ID"].nunique()
print(f"SMM risk strata recovered from barcodes: "
      f"{_r.get('high', 0)} high-risk, {_r.get('low', 0)} low-risk samples")

# ------------------------------------------------------------------ QC
# The dataset is pre-curated, so QC is permissive: we only remove cells with
# very few detected genes or very high mitochondrial fraction.
qc_pass = (meta["n_genes"] > 200) & (meta["frac_mito"] < 0.25)
print(f"QC pass: {qc_pass.sum()} / {len(meta)} cells")
meta = meta[qc_pass]
counts = counts[meta.index]

# ---------------------------------------------------------------- build AnnData
# AnnData wants cells x genes, so transpose the counts.
# Convert to float32 FIRST and free the float64 copy to keep memory low.
counts_f32 = counts.T.astype("float32")
del counts

adata = sc.AnnData(counts_f32)
adata.obs = meta.copy()

# Genes we never want as vaccine targets: immunoglobulins, HLA, mito genes.
# We keep them in the object (useful for QC) but tag them. The filter itself
# lives in config.py so steps 1 and 3 cannot drift apart.
adata.var["excluded_from_targets"] = [is_excluded(g) for g in adata.var_names]
print(f"excluded gene tags: {adata.var['excluded_from_targets'].sum()}")

adata.write_h5ad(f"{DATA}/smm_plasma_cells.h5ad")

# ------------------------------------------------------------- summary print
print("\nCells per disease stage:")
print(adata.obs["disease_stage"].value_counts())
print("\nCells per normal/neoplastic label:")
print(adata.obs["normal_or_neoplastic"].value_counts())
print("\nSamples per stage:")
print(adata.obs.groupby("disease_stage")["sample_ID"].nunique())
print("\nSaved:", f"{DATA}/smm_plasma_cells.h5ad")
