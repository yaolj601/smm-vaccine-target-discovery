"""
Step 6: Composite ranking and final shortlist.

Merges all evidence and computes a weighted composite score (each component
scaled 0-1):

  1. DE evidence (weight 0.30):
     log2FC (SMM vs healthy, capped at 5) x within-patient recurrence x
     specificity (fraction of neoplastic cells expressing MINUS fraction of
     normal cells). Recurrence, not detection-prevalence: a gene detected in
     every patient is not necessarily enriched in any of them.
  2. HLA-I presentation (weight 0.20):
     allele coverage (fraction of the 27-allele panel with a strong binder)
     x epitope density (strong binders per 100 aa - length-normalized so big
     proteins don't win just by being big)
  3. Normal-tissue safety (weight 0.25):
     NOT tumour specificity - that is the third term inside the DE component
     above, and it compares tumour plasma cells against the same patient's
     normal plasma cells. This asks the different question of where ELSE in
     the body the gene is expressed, from the tier assigned in step 2b by
     outlier testing against ~50 HPA consensus tissues: 1a confined to bone
     marrow/lymphoid, 1b confined to testis (cancer-testis antigen), 2
     restricted to some other non-vital tissue, 4 broadly expressed. Neural
     and vital-organ outliers are rejected in step 3 and never reach here.
     Both questions are needed: CADM1 was 38% of neoplastic cells vs 11% of
     normal - excellent tumour specificity - and retina-enriched.
  4. novelty_confidence (weight 0.25):
     two categorical LLM judgements - is this plausibly a myeloma antigen
     (A-D), and is it still unclaimed (development tier) - combined by the
     lookup grid in step 5b. Genes marked confidence D
     (housekeeping / broadly essential) are REJECTED here rather than
     penalised: a known housekeeping gene does not belong on a vaccine
     shortlist at any rank, whatever its expression data says. The D gate
     sits in this step because step 5 runs after epitope prediction.

Input:  $SMM_DATA_DIR/{candidates_annotated,epitope_results,de_results}.csv
        $SMM_DATA_DIR/llm_scores.json (optional)
Output: $SMM_OUT_DIR/ranked_targets.csv
        $SMM_OUT_DIR/ranked_targets.md
        $SMM_OUT_DIR/all_candidates_scored.csv
        $SMM_OUT_DIR/evidence/*            (intermediate tables, audit trail)
        $SMM_OUT_DIR/figures/*.png|.svg
"""

import json
import os
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


from config import DATA, OUT

FIG = f"{OUT}/figures"
EVIDENCE = f"{OUT}/evidence"
os.makedirs(FIG, exist_ok=True)
os.makedirs(EVIDENCE, exist_ok=True)

matplotlib.rcParams["font.family"] = ["Liberation Sans", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"

# Epitope prediction gets the lowest weight: it is the least validated layer
# (affinity models, no elution data, self-tolerance unmodeled). LLM knowledge
# (existing programs, safety signals) is weighted equal to safety.
WEIGHTS = {"de": 0.30, "epitope": 0.20, "safety": 0.25, "llm": 0.25}
WEIGHT_LABELS = {"de": "DE", "epitope": "HLA presentation",
                 "safety": "tissue safety", "llm": "novelty-confidence"}
assert abs(sum(WEIGHTS.values()) - 1.0) < 1e-9, "composite weights must sum to 1"

# Single source of truth for the human-readable formula, so the number printed
# in ranked_targets.md can never drift from the number actually applied.
WEIGHT_FORMULA = " + ".join(f"{WEIGHTS[k]:.2f}*{WEIGHT_LABELS[k]}" for k in WEIGHTS)

# ------------------------------------------------------------- load evidence
annot = pd.read_csv(f"{DATA}/candidates_annotated.csv")
epi = pd.read_csv(f"{DATA}/epitope_results.csv")
de = pd.read_csv(f"{DATA}/de_results.csv")

de_b = de[de["contrast_name"] == "B_SMM_vs_NBM"].set_index("gene")
de_a = de[de["contrast_name"] == "A_neoplastic_vs_normal"].set_index("gene")
de_c = de[de["contrast_name"] == "C_paired"].set_index("gene")

df = annot.merge(epi, on="gene", how="inner")

# per-cell stats (identical across contrasts; take from A)
for col in ["pct_cells_neoplastic", "pct_cells_normal",
            "pct_cells_smm_neoplastic", "pct_cells_nbm_normal", "patient_prevalence",
            "recurrence_frac", "recurrence_n_samples", "recurrence_denom",
            "recurrence_smm_n", "recurrence_smm_denom",
            "mgus_coverage_n", "mgus_coverage_denom", "smm_coverage_n",
            "smm_coverage_denom", "mm_coverage_n", "mm_coverage_denom",
            "smm_coverage_frac", "mm_retained",
            "smmh_coverage_n", "smmh_coverage_denom", "smmh_coverage_frac",
            "smml_coverage_n", "smml_coverage_denom"]:
    df[col] = df["gene"].map(de_a[col])

# DE stats from the SMM contrast (the vaccine indication)
df["log2fc_smm"] = df["gene"].map(de_b["log2FoldChange"])
df["padj_smm"] = df["gene"].map(de_b["padj"])
# The score ranks on contrast B's fold change, but a gene may have qualified on
# A or C instead. Flag where the scored quantity comes from a contrast the gene
# did not itself clear, so that is visible rather than buried.
df["significant_in_B"] = df["padj_smm"] < 0.05
# paired-contrast support (within-patient evidence)
df["log2fc_paired"] = df["gene"].map(de_c["log2FoldChange"])
df["padj_paired"] = df["gene"].map(de_c["padj"])

# drop genes with no protein (lncRNAs etc.) - already dropped by the merge
print(f"{len(df)} candidates with protein + epitope evidence")

# ------------------------------------------------------- component scores
def scale(s):
    """Scale a series to 0-1 (min-max)."""
    s = s.astype(float)
    if s.max() == s.min():
        return pd.Series(0.5, index=s.index)
    return (s - s.min()) / (s.max() - s.min())

# 1. DE evidence
# Indication-matched: SMM tumour cells vs healthy-donor plasma cells, the same
# two groups contrast B compares. The pooled version is 57% MM cells and would
# reward genes whose specificity is an MM phenomenon.
df["specificity"] = (df["pct_cells_smm_neoplastic"]
                     - df["pct_cells_nbm_normal"]).clip(lower=0)
df["de_component"] = (df["log2fc_smm"].clip(lower=0, upper=5) / 5
                      * df["recurrence_frac"]
                      * df["specificity"])
df["de_component"] = scale(df["de_component"])

# 2. HLA presentation (length-normalized)
df["epitope_density"] = df["n_strong_binders"] / df["protein_length_aa"] * 100
df["epitope_component"] = (scale(df["allele_coverage_frac"])
                           * scale(df["epitope_density"]))
df["epitope_component"] = scale(df["epitope_component"])

# 3. Safety: HPA tissue restriction (vital-tissue cap applied in step 3)
df["safety_component"] = scale(df["hpa_safety_score"])

# 4. novelty_confidence, from step 5b
_nc = f"{DATA}/novelty_confidence.csv"
if not os.path.exists(_nc):
    raise SystemExit(f"missing {_nc}\nRun 05b_novelty_confidence.py first.")
nc = pd.read_csv(_nc).set_index("gene")
for col in ["mm_confidence", "novelty_tier", "novelty_confidence", "rejected",
            "existing_programs"]:
    df[col] = df["gene"].map(nc[col])
df["existing_programs"] = df["existing_programs"].fillna("")
df["llm_available"] = df["mm_confidence"].notna() & df["novelty_tier"].notna()

# confidence D is disqualifying, not merely low-scoring: a known housekeeping
# gene does not belong on a vaccine shortlist at any rank, whatever its
# expression data says - the same logic that rejects vital-tissue outliers.
n_before = len(df)
df = df[df["rejected"] != True].reset_index(drop=True)
if n_before != len(df):
    print(f"rejected {n_before - len(df)} housekeeping/essential genes (confidence D)")

# an unassessed gene sits mid-grid rather than being rewarded or buried
df["novelty_confidence"] = df["novelty_confidence"].astype(float).fillna(4.0)
df["llm_component"] = scale(df["novelty_confidence"])

# ------------------------------------------------------- composite + rank
df["composite_score"] = (WEIGHTS["de"] * df["de_component"]
                         + WEIGHTS["epitope"] * df["epitope_component"]
                         + WEIGHTS["safety"] * df["safety_component"]
                         + WEIGHTS["llm"] * df["llm_component"])

df = df.sort_values("composite_score", ascending=False).reset_index(drop=True)
df["rank"] = df.index + 1

# rationale per gene (assembled from the evidence columns)
def rationale(row):
    parts = []
    parts.append(f"SMM neoplastic vs healthy plasma cells log2FC="
                 f"{row['log2fc_smm']:.1f} (padj={row['padj_smm']:.1e})"
                 if pd.notna(row["log2fc_smm"]) else "DE n/a")
    parts.append(f"expressed in {row['pct_cells_smm_neoplastic']:.0%} of SMM tumour cells "
                 f"vs {row['pct_cells_nbm_normal']:.0%} of healthy-donor plasma cells")
    parts.append(f"shared by {row['patient_prevalence']:.0%} of patients")
    parts.append(f"enriched vs own normal PCs in {int(row['recurrence_n_samples'])}"
                 f"/{int(row['recurrence_denom'])} paired patients "
                 f"({int(row['recurrence_smm_n'])}/{int(row['recurrence_smm_denom'])} SMM)")
    parts.append(f"expressed (>=5 CPM) in {int(row['smm_coverage_n'])}"
                 f"/{int(row['smm_coverage_denom'])} SMM tumours "
                 f"({int(row['smmh_coverage_n'])}/{int(row['smmh_coverage_denom'])} "
                 f"high-risk), "
                 f"{int(row['mm_coverage_n'])}/{int(row['mm_coverage_denom'])} MM, "
                 f"{int(row['mgus_coverage_n'])}/{int(row['mgus_coverage_denom'])} MGUS")
    if not row["mm_retained"]:
        parts.append("NOT retained on progression to MM (durability risk)")
    parts.append(f"{row['n_strong_binders']} strong HLA-I binders on "
                 f"{row['n_alleles_with_binder']}/27 alleles")
    if not row["significant_in_B"]:
        parts.append("NOT significant in the SMM contrast (entered via A or C)")
    parts.append(f"novelty-confidence {row['novelty_confidence']:.0f}/10 "
                 f"(confidence {row['mm_confidence']}, {row['novelty_tier']})")
    parts.append(f"tissue tier {row['tissue_tier']}"
                 + (f" (outlier in {row['hpa_outlier_tissues']})"
                    if row["hpa_outlier_tissues"] else ""))
    if pd.notna(row["log2fc_paired"]) and row["log2fc_paired"] > 0.5:
        parts.append(f"also up in within-patient paired contrast (log2FC={row['log2fc_paired']:.1f})")
    return "; ".join(parts)

df["rationale"] = df.apply(rationale, axis=1)
shortlist = df.head(15).copy()

# ------------------------------------------------------------- save tables
cols = ["rank", "gene", "composite_score", "de_component", "epitope_component",
        "safety_component", "llm_component", "log2fc_smm", "padj_smm",
        "pct_cells_smm_neoplastic", "pct_cells_nbm_normal",
        "pct_cells_neoplastic", "pct_cells_normal", "patient_prevalence",
        "recurrence_frac", "recurrence_n_samples", "recurrence_smm_n",
        "smm_coverage_n", "smmh_coverage_n", "smml_coverage_n",
        "mm_coverage_n", "mgus_coverage_n", "mm_retained",
        "n_strong_binders", "n_alleles_with_binder", "allele_coverage_frac",
        "best_affinity_nM", "tissue_tier", "hpa_outlier_tissues",
        "hpa_lymphoid_fdr", "hpa_subcellular", "localisation",
        "significant_in_B", "surface_or_secreted", "novelty_confidence", "mm_confidence",
        "novelty_tier", "existing_programs", "llm_available", "rationale"]
shortlist[cols].to_csv(f"{OUT}/ranked_targets.csv", index=False)
df[cols].to_csv(f"{OUT}/all_candidates_scored.csv", index=False)

with open(f"{OUT}/ranked_targets.md", "w") as f:
    f.write("# Ranked vaccine target shortlist (SMM)\n\n")
    f.write(f"Composite score = {WEIGHT_FORMULA}\n\n")
    f.write("| Rank | Gene | Score | Rationale |\n|---|---|---|---|\n")
    for _, r in shortlist.iterrows():
        f.write(f"| {r['rank']} | {r['gene']} | {r['composite_score']:.3f} | {r['rationale']} |\n")

print(shortlist[["rank", "gene", "composite_score"]].to_string(index=False))

# ------------------------------------------------------------- figures
# Fig 1: composite score stacked by component (top 15)
top = shortlist.iloc[::-1]  # reverse for horizontal bars
fig, ax = plt.subplots(figsize=(8, 6))
ax.barh(top["gene"], top["de_component"] * WEIGHTS["de"], color="#0279EE", label=WEIGHT_LABELS["de"])
ax.barh(top["gene"], top["epitope_component"] * WEIGHTS["epitope"],
        left=top["de_component"] * WEIGHTS["de"], color="#75A025", label=WEIGHT_LABELS["epitope"])
ax.barh(top["gene"], top["safety_component"] * WEIGHTS["safety"],
        left=top["de_component"] * WEIGHTS["de"] + top["epitope_component"] * WEIGHTS["epitope"],
        color="#FF9400", label=WEIGHT_LABELS["safety"])
ax.barh(top["gene"], top["llm_component"] * WEIGHTS["llm"],
        left=top["de_component"] * WEIGHTS["de"] + top["epitope_component"] * WEIGHTS["epitope"]
        + top["safety_component"] * WEIGHTS["safety"],
        color="#FD9BED", label=WEIGHT_LABELS["llm"])
ax.set_xlabel("Weighted component contribution")
ax.set_title("Top 15 candidate vaccine targets")
ax.legend(loc="lower right", fontsize=8)
plt.tight_layout()
plt.savefig(f"{FIG}/composite_scores.png", dpi=200)
plt.savefig(f"{FIG}/composite_scores.svg")
plt.close()

# Fig 2: specificity vs HLA coverage, labeled scatter
fig, ax = plt.subplots(figsize=(7, 6))
sc_ax = ax.scatter(df["specificity"], df["allele_coverage_frac"],
                   c=df["composite_score"], cmap="viridis", s=60, alpha=0.8)
for _, r in df.head(15).iterrows():
    ax.annotate(r["gene"], (r["specificity"], r["allele_coverage_frac"]),
                fontsize=7, xytext=(3, 3), textcoords="offset points")
ax.set_xlabel("Tumour specificity: % neoplastic cells expressing - % normal plasma cells\n(within-marrow; distinct from the normal-tissue safety component)")
ax.set_ylabel("HLA allele coverage (fraction of 27-allele panel)")
plt.colorbar(sc_ax, label="Composite score")
plt.tight_layout()
plt.savefig(f"{FIG}/specificity_vs_hla.png", dpi=200)
plt.savefig(f"{FIG}/specificity_vs_hla.svg")
plt.close()

# Fig 3: heatmap of key evidence for top 15
# the six quantities the composite actually runs on. The earlier version
# showed patient_prevalence (the superseded >1-UMI detection metric, only 6
# distinct values across the candidate set) and allele_coverage_frac (still
# length-confounded) while omitting recurrence, coverage and novelty entirely.
heat_cols = ["log2fc_smm", "specificity", "recurrence_frac", "smm_coverage_frac",
             "epitope_density", "novelty_confidence"]
heat = shortlist.set_index("gene")[heat_cols]
heat_norm = (heat - heat.min()) / (heat.max() - heat.min())
fig, ax = plt.subplots(figsize=(7, 7))
sns.heatmap(heat_norm, annot=heat.round(2), fmt="", cmap="YlGnBu",
            cbar_kws={"label": "min-max scaled"}, ax=ax)
ax.set_title("Evidence matrix, top 15 targets")
plt.tight_layout()
plt.savefig(f"{FIG}/evidence_heatmap.png", dpi=200)
plt.savefig(f"{FIG}/evidence_heatmap.svg")
plt.close()

# --------------------------------------------------- evidence trail (audit)
# Copy the intermediate tables next to the headline result so a reader can
# trace every number from the raw matrix through to the ranked shortlist
# without re-running the pipeline.
copied = []
for name in ["de_results.csv", "wilcoxon_sensitivity.csv", "tissue_tiers.csv",
             "candidates_annotated.csv", "epitope_results.csv",
             "novelty_confidence.csv", "llm_scores.json", "llm_scores_insession.json"]:
    src = f"{DATA}/{name}"
    if os.path.exists(src):
        shutil.copy2(src, f"{EVIDENCE}/{name}")
        copied.append(name)
print(f"\nEvidence trail -> {EVIDENCE}: {', '.join(copied) if copied else '(nothing found)'}")

print(f"\nSaved ranked_targets.csv, ranked_targets.md, all_candidates_scored.csv "
      f"and 3 figures to {OUT}")
print(f"Composite score = {WEIGHT_FORMULA}")
