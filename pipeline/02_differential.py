"""
Step 2: Differential expression - find genes enriched in neoplastic plasma cells.

Approach: "pseudobulk" - add up all UMI counts per (sample, cell status).
This gives one count vector per sample per status, which is the statistically
correct unit of replication (cells are not independent).

Contrasts:
  A. neoplastic vs normal plasma cells (all samples)
  B. SMM neoplastic vs NBM normal          <- the vaccine indication
  C. paired: neoplastic vs normal within the same patients

A fourth contrast (all-neoplastic vs NBM-normal, "pan-tumour signal") was
removed: it correlated r=0.957 with contrast B, flagged all 15 shortlisted
genes exactly as B did, and was never consumed downstream. The remaining three
are relabelled A/B/C so the sequence has no gap.

Also computes per-gene vaccine-relevant stats:
  - pct_cells_expressing: % of neoplastic cells with >0 counts
  - patient_prevalence: % of patients whose neoplastic cells DETECT the gene
    (the OFF-THE-SHELF criterion: the target must be present in most patients)
  - {mgus,smm,mm}_coverage_frac: % of tumour samples at that stage expressing
    the gene above a real threshold (5 CPM), not merely detecting it.
    SMM coverage is the off-the-shelf gate that actually matters - an antigen
    absent from most SMM tumours cannot work as a shared vaccine. MM coverage
    is the DURABILITY check: the vaccine is given at SMM to block progression,
    so an antigen lost on the way to MM fails exactly when it is needed.
    MGUS coverage is reported but never gated on - only 3 MGUS samples carry
    tumour cells (307 cells; two under 70), so at n=3 it measures which
    translocation subtypes happened to be sampled rather than stage biology
    (CCND1 is 683 / 0.7 / 313 CPM across the three). See WRITEUP section 1.
  - recurrence_frac: % of paired patients in which the gene is actually HIGHER
    in neoplastic than in that same patient's normal plasma cells
    (the REPRODUCIBILITY criterion - detection alone badly overstates
    enrichment; see WRITEUP section 1)

Input:  $SMM_DATA_DIR/smm_plasma_cells.h5ad
Output: $SMM_DATA_DIR/de_results.csv, $SMM_DATA_DIR/wilcoxon_sensitivity.csv
"""

import numpy as np
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from config import DATA

print("Loading AnnData...")
adata = sc.read_h5ad(f"{DATA}/smm_plasma_cells.h5ad")
counts_sc = pd.DataFrame(adata.X.T, index=adata.var_names, columns=adata.obs_names)
meta = adata.obs.copy()

# ------------------------------------------------------- per-gene cell stats
print("Computing per-gene expression stats from single cells...")
neoplastic_cells = counts_sc.loc[:, meta["normal_or_neoplastic"] == "neoplastic"]
normal_cells = counts_sc.loc[:, meta["normal_or_neoplastic"] == "normal"]
patient_samples = [s for s in meta["sample_ID"].unique() if not s.startswith("NBM")]

pct_neoplastic = (neoplastic_cells > 0).mean(axis=1)
pct_normal = (normal_cells > 0).mean(axis=1)

# The pooled fractions above are 57% MM cells and 42% SMM, and the normal pool
# is 90% healthy-donor cells - so a gene's apparent "tumour specificity" can be
# an MM phenomenon. IFITM1 is expressed in 46% of pooled neoplastic cells but
# only 21% of SMM ones. The specificity term in the score therefore uses the
# indication-matched pair instead: SMM tumour cells against healthy-donor
# plasma cells, the same two groups contrast B compares. The pooled fractions
# stay as reported columns.
_smm_neo = ((meta["normal_or_neoplastic"] == "neoplastic")
            & (meta["disease_stage"] == "SMM")).values
_nbm_nor = ((meta["normal_or_neoplastic"] == "normal")
            & (meta["disease_stage"] == "NBM")).values
pct_smm_neoplastic = (counts_sc.loc[:, _smm_neo] > 0).mean(axis=1)
pct_nbm_normal = (counts_sc.loc[:, _nbm_nor] > 0).mean(axis=1)
print(f"cell fractions: pooled {neoplastic_cells.shape[1]} neoplastic / "
      f"{normal_cells.shape[1]} normal; indication-matched "
      f"{int(_smm_neo.sum())} SMM-tumour / {int(_nbm_nor.sum())} healthy-donor")

# patient prevalence: gene detected (>1 UMI total) in neoplastic pseudobulk
# of at least X% of patients. Some low-burden patients (esp. MGUS) have no
# neoplastic cells at all, so we only count patients that do.
neo_pseudobulk = neoplastic_cells.T.groupby(meta.loc[neoplastic_cells.columns, "sample_ID"]).sum()
patients_with_tumor = [s for s in patient_samples if s in neo_pseudobulk.index]
print(f"patients with neoplastic cells: {len(patients_with_tumor)} / {len(patient_samples)}")
detected_per_patient = (neo_pseudobulk.loc[patients_with_tumor] > 1)
patient_prevalence = detected_per_patient.mean(axis=0)

# --------------------------------------------- within-patient recurrence
# "Is this gene genuinely enriched in the tumour, in THIS patient?" - asked
# once per sample that has both neoplastic and normal plasma cells, then
# counted across samples. This is the direction half of the recurrence
# criterion in Yao et al. (Cancer Res 2023); their Bonferroni-significance
# half is unreachable at this depth (WRITEUP section 1 explains why).
#
# Two gates, both necessary. Direction alone is NOT enough: with only 6-166
# normal plasma cells per sample, a near-tie on near-zero values counts as
# "enriched" (NEB in MGUS-3 is 1.3 vs 1.2 CPM; ATP10B in SMM-8 is 0.2 vs 0.0).
# That artefact put two length/low-expression candidates into the top 15.
#   - RECUR_LOG2FC mirrors Seurat's logfc.threshold, which Yao et al. used.
#   - RECUR_MIN_CPM additionally requires the gene to be genuinely expressed in
#     the tumour compartment, which the +1 pseudocount alone does not enforce.
# Together these cut genome-wide perfect-recurrence genes from 112 to 7.
RECUR_LOG2FC = 0.25
RECUR_MIN_CPM = 5.0

print("Computing within-patient recurrence...")
per_sample = pd.crosstab(meta["sample_ID"], meta["normal_or_neoplastic"])
paired_samples = per_sample[(per_sample.get("normal", 0) > 0)
                            & (per_sample.get("neoplastic", 0) > 0)].index.tolist()

# The 11 paired samples are MGUS 3 / SMM 5 / MM 3, so a gene's recurrence
# could in principle be carried by the non-SMM pairs. Recurrence COUNTS rather
# than pools, so one MM patient is one vote and cannot dominate the way a
# pooled regression can - empirically no candidate has zero SMM-pair support,
# minimum is 2 of 5. Reported rather than gated for exactly that reason: a gate
# here would filter nothing and imply a check that is not doing any work.
_stage = meta.groupby("sample_ID", observed=True)["disease_stage"].first()
smm_paired = [s for s in paired_samples if _stage[s] == "SMM"]

higher = pd.Series(0, index=counts_sc.index, dtype=int)
higher_smm = pd.Series(0, index=counts_sc.index, dtype=int)
for s in paired_samples:
    in_s = (meta["sample_ID"] == s)
    neo = counts_sc.loc[:, (in_s & (meta["normal_or_neoplastic"] == "neoplastic")).values]
    nor = counts_sc.loc[:, (in_s & (meta["normal_or_neoplastic"] == "normal")).values]
    # CPM within the sample, so library size does not drive the comparison
    cpm_neo = neo.sum(axis=1) / max(neo.values.sum(), 1) * 1e6
    cpm_nor = nor.sum(axis=1) / max(nor.values.sum(), 1) * 1e6
    lfc = np.log2((cpm_neo + 1) / (cpm_nor + 1))
    hit = ((lfc >= RECUR_LOG2FC) & (cpm_neo >= RECUR_MIN_CPM)).astype(int)
    higher += hit
    if _stage[s] == "SMM":
        higher_smm += hit

recurrence_n = higher
recurrence_frac = higher / len(paired_samples)
print(f"recurrence computed over {len(paired_samples)} paired samples "
      f"({len(smm_paired)} of them SMM): "
      f"{(recurrence_frac >= 0.8).sum()} genes enriched in >=80% of them")

# ------------------------------------------------- per-stage tumour coverage
# "In how many tumours of this stage is the antigen actually ON?" - a real
# expression threshold, not the >1 UMI detection bar used by patient_prevalence
# (which is why that metric rated EDNRB as broadly shared when it clears 5 CPM
# in only 3 of 12 SMM tumours).
COVERAGE_CPM = 5.0

print("Computing per-stage tumour coverage...")
neo_mask = (meta["normal_or_neoplastic"] == "neoplastic")
stage_of = meta.groupby("sample_ID", observed=True)["disease_stage"].first()
# SMMh / SMMl are reported separately but NOT gated on. High-risk SMM (~50%
# progression at 2 years) is the population a vaccine would actually enrol, so
# coverage there is the honest denominator - but n=9 vs n=12 is a real loss of
# power, and the choice of which to gate on is a scope decision, not a
# technical one. Both are reported so the tradeoff is visible.
risk_of = meta.groupby("sample_ID", observed=True)["smm_risk"].first()
STRATA = {"MGUS": lambda s: stage_of[s] == "MGUS",
          "SMM":  lambda s: stage_of[s] == "SMM",
          "MM":   lambda s: stage_of[s] == "MM",
          "SMMh": lambda s: stage_of[s] == "SMM" and risk_of.get(s) == "high",
          "SMMl": lambda s: stage_of[s] == "SMM" and risk_of.get(s) == "low"}

coverage = {}
for stg, belongs in STRATA.items():
    samples = [s for s in meta["sample_ID"].unique()
               if belongs(s) and (neo_mask & (meta["sample_ID"] == s)).sum() > 0]
    hits = pd.Series(0, index=counts_sc.index, dtype=int)
    for s in samples:
        X = counts_sc.loc[:, (neo_mask & (meta["sample_ID"] == s)).values]
        hits += (X.sum(axis=1) / max(X.values.sum(), 1) * 1e6 >= COVERAGE_CPM).astype(int)
    coverage[stg] = (hits, len(samples))
    print(f"  {stg}: {len(samples)} tumour samples")

# durability: is the antigen retained as disease progresses SMM -> MM?
smm_f = coverage["SMM"][0] / coverage["SMM"][1]
mm_f = coverage["MM"][0] / coverage["MM"][1]
mm_retained = mm_f >= (smm_f - 0.2)

# ------------------------------------------------------------- pseudobulk
print("Building pseudobulk matrix...")
group_key = meta["sample_ID"].astype(str) + "|" + meta["normal_or_neoplastic"].astype(str)
pseudobulk = counts_sc.T.groupby(group_key).sum().T  # genes x groups
pseudobulk = pseudobulk.loc[:, pseudobulk.sum(axis=0) > 0]

# pseudobulk metadata table
pb_meta = pd.DataFrame(
    [g.split("|") for g in pseudobulk.columns],
    index=pseudobulk.columns,
    columns=["sample_ID", "status"],
)
pb_meta["stage"] = [meta.loc[meta["sample_ID"] == s, "disease_stage"].iloc[0]
                    for s in pb_meta["sample_ID"]]
print(f"pseudobulk: {pseudobulk.shape[0]} genes x {pseudobulk.shape[1]} sample-groups")

# keep genes with at least 10 total counts (DESeq2 needs low-count filtering)
keep = pseudobulk.sum(axis=1) >= 10
pb = pseudobulk[keep].T.astype(int)  # samples x genes, integer counts


def run_deseq2(counts_df, meta_df, design, contrast, name):
    """Run one DESeq2 contrast and return a results DataFrame."""
    dds = DeseqDataSet(counts=counts_df, metadata=meta_df, design=design)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=contrast)
    stats.summary()
    res = stats.results_df
    res["design"] = design
    res["contrast"] = str(contrast)
    res["contrast_name"] = name
    return res


results = []

# Contrast A: neoplastic vs normal, all samples
print("Contrast A: neoplastic vs normal (all samples)...")
meta_a = pb_meta[["status"]].copy()
results.append(run_deseq2(pb, meta_a, "~status", ["status", "neoplastic", "normal"], "A_neoplastic_vs_normal"))

# Contrast B: SMM neoplastic vs NBM normal (the vaccine indication)
print("Contrast B: SMM-neoplastic vs NBM-normal...")
sel_b = ((pb_meta["stage"] == "SMM") & (pb_meta["status"] == "neoplastic")) | \
        ((pb_meta["stage"] == "NBM") & (pb_meta["status"] == "normal"))
sub_b = pb[sel_b]
meta_b = pb_meta[sel_b][["status"]].copy()
results.append(run_deseq2(sub_b, meta_b, "~status", ["status", "neoplastic", "normal"], "B_SMM_vs_NBM"))

# Contrast C: paired within patients (only patient samples that have BOTH
# normal and neoplastic pseudobulks)
print("Contrast C: paired neoplastic vs normal within patients...")
both = pb_meta.groupby("sample_ID")["status"].nunique()
paired_samples = both[both == 2].index
sel_d = pb_meta["sample_ID"].isin(paired_samples)
sub_d = pb[sel_d]
meta_d = pb_meta[sel_d].copy()
meta_d["patient"] = meta_d["sample_ID"]  # one biopsy per sample ID
results.append(run_deseq2(sub_d, meta_d, "~patient + status",
                          ["status", "neoplastic", "normal"], "C_paired"))

# ------------------------------------------------------- merge and save
de = pd.concat(results)
de.index.name = "gene"
de = de.reset_index()

# attach the per-gene vaccine stats
de["pct_cells_neoplastic"] = de["gene"].map(pct_neoplastic)
de["pct_cells_normal"] = de["gene"].map(pct_normal)
de["pct_cells_smm_neoplastic"] = de["gene"].map(pct_smm_neoplastic)
de["pct_cells_nbm_normal"] = de["gene"].map(pct_nbm_normal)
de["patient_prevalence"] = de["gene"].map(patient_prevalence)
de["recurrence_n_samples"] = de["gene"].map(recurrence_n)
de["recurrence_frac"] = de["gene"].map(recurrence_frac)
de["recurrence_denom"] = len(paired_samples)
de["recurrence_smm_n"] = de["gene"].map(higher_smm)
de["recurrence_smm_denom"] = len(smm_paired)
for stg, (hits, n) in coverage.items():
    de[f"{stg.lower()}_coverage_n"] = de["gene"].map(hits)
    de[f"{stg.lower()}_coverage_frac"] = de["gene"].map(hits / n)
    de[f"{stg.lower()}_coverage_denom"] = n
de["mm_retained"] = de["gene"].map(mm_retained)

de.to_csv(f"{DATA}/de_results.csv", index=False)
print(f"saved {len(de)} DE rows to de_results.csv")

# ------------------------------------------- sensitivity check: Wilcoxon
print("Wilcoxon rank-sum sensitivity check (single-cell level)...")
adata.obs["status"] = adata.obs["normal_or_neoplastic"]
sc.tl.rank_genes_groups(adata, groupby="status", groups=["neoplastic"],
                        reference="normal", method="wilcoxon")
wil = sc.get.rank_genes_groups_df(adata, group="neoplastic")
wil.to_csv(f"{DATA}/wilcoxon_sensitivity.csv", index=False)
print("saved wilcoxon_sensitivity.csv")

# quick sanity check: known myeloma/cancer-testis genes should be up
print("\nSanity check (contrast A log2FC for known MM/CT antigens):")
a = de[de["contrast_name"] == "A_neoplastic_vs_normal"].set_index("gene")
for g in ["PRAME", "MAGEA3", "MAGEA6", "XAGE1", "GAGE4", "CTAG1B", "BCMA", "TNFRSF17"]:
    if g in a.index:
        print(f"  {g}: log2FC={a.loc[g, 'log2FoldChange']:.2f}, padj={a.loc[g, 'padj']:.2e}")
