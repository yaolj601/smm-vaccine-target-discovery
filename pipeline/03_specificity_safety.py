"""
Step 3: Tumor specificity and safety annotation for DE candidates.

For each candidate gene we add three evidence columns:

  1. Normal-tissue restriction, by outlier testing against HPA consensus RNA
     across ~50 tissues (step 2b). Not just "is it restricted?" but
     "restricted to WHICH tissue?" - genes whose normal expression is confined
     to bone marrow / lymphoid tissue, or to testis (cancer-testis antigens),
     are good vaccine targets; genes that are outliers in neural or vital
     tissue are rejected outright.
  HPA's own `RNA tissue specificity` classification is deliberately NOT
  carried: it answers "how many tissues" where step 2b answers "which
  tissues", the two disagree in the cases that matter (CADM1 is "Tissue
  enhanced" by HPA and retina-enriched by the outlier test), and keeping both
  invited the reader to assume the wrong one drove the score.

  2. HPA subcellular location + protein class:
     is the protein cell-surface or secreted (bonus: also enables antibodies)?
DepMap CRISPR essentiality was deliberately dropped (see WRITEUP section 1).
Short version: it answers "does the tumour die without this gene", which is the
question for a degrader or small molecule, not for a vaccine. Its escape
argument also inverts - an essential antigen is HARDER to lose under immune
pressure - and its normal-toxicity signal duplicates HPA tissue specificity.
On this dataset it changed nothing: the top 15 were identical with and without.

Input:  $SMM_DATA_DIR/de_results.csv
        $SMM_DATA_DIR/tissue_tiers.csv            (from step 2b)
        $SMM_REF_DIR/proteinatlas.tsv             (see README for download)
Output: $SMM_DATA_DIR/candidates_annotated.csv
"""

from pathlib import Path

import pandas as pd

from config import DATA, hpa_tsv, is_excluded

# tier -> safety score. Tier 3 is gated out below; the value is kept so a gene
# that somehow reaches scoring cannot benefit from it. "no_hpa" is penalised at
# 0.3 rather than treated as neutral: for a preventive vaccine, "cannot assess"
# is not the same as "probably fine".
TIER_SAFETY = {"1a_lymphoid": 1.00, "1b_testis_CT": 1.00, "2_restricted": 0.60,
               "4_broad": 0.20, "no_hpa": 0.30, "3_vital": 0.00}
REJECT_TIERS = {"3_vital"}

# ------------------------------------------------- build candidate universe
print("Selecting candidate genes from DE results...")
de = pd.read_csv(f"{DATA}/de_results.csv")

# genes passing DE filters in ANY of the main contrasts
cand_a = de[(de["contrast_name"] == "A_neoplastic_vs_normal")
            & (de["padj"] < 0.05) & (de["log2FoldChange"] > 1)]
cand_b = de[(de["contrast_name"] == "B_SMM_vs_NBM")
            & (de["padj"] < 0.05) & (de["log2FoldChange"] > 1)]
cand_c = de[(de["contrast_name"] == "C_paired")
            & (de["padj"] < 0.05) & (de["log2FoldChange"] > 0.5)]

candidates = sorted(set(cand_a["gene"]) | set(cand_b["gene"]) | set(cand_c["gene"]))
print(f"union of DE candidates: {len(candidates)} genes")

# INDICATION CHECK. Contrast A pools every stage, and MM contributes both more
# cells (10,790 vs 8,431) and more extreme tumours, so a gene can enter on
# MM-driven signal alone - A-only genes have median SMM coverage 0.17 against
# MM coverage 0.38. This is an SMM vaccine, not an MM vaccine, so a gene
# admitted by A alone must at least point the right way in the indication
# contrast. Deliberately a directional bar, not a significance one: at 12 vs 9
# samples B cannot confirm genes that are genuinely up in SMM (136 of the
# A-only set have B log2FC > 1 yet miss padj), and excluding those would trade
# real recall for tidiness.
MIN_B_LFC = 0.5
b_lfc = de[de["contrast_name"] == "B_SMM_vs_NBM"].set_index("gene")["log2FoldChange"]
a_only = set(cand_a["gene"]) - set(cand_b["gene"]) - set(cand_c["gene"])
fails = {g for g in a_only if not (b_lfc.get(g, float("-inf")) > MIN_B_LFC)}
candidates = [g for g in candidates if g not in fails]
print(f"after requiring A-only genes to show log2FC > {MIN_B_LFC} in SMM: "
      f"{len(candidates)} genes ({len(fails)} dropped as MM-driven)")

# drop immunoglobulin / HLA / mitochondrial genes (tagged in step 1;
# is_excluded is shared with step 1 via config.py)
candidates = [g for g in candidates if not is_excluded(g)]
print(f"after excluding IG/HLA/mito genes: {len(candidates)} genes")

# OFF-THE-SHELF GATE: the antigen must actually be ON in at least half of the
# SMM tumours. This replaces the old `patient_prevalence >= 0.3` gate, which
# used a >1-UMI detection bar across all stages and therefore passed genes that
# are barely expressed anywhere (EDNRB reached rank 5 on 3/12 SMM coverage).
# patient_prevalence is still reported, as detection is a different question.
MIN_SMM_COVERAGE = 0.5
stats = de[de["contrast_name"] == "A_neoplastic_vs_normal"].set_index("gene")
keep = [g for g in candidates
        if g in stats.index
        and stats.loc[g, "smm_coverage_frac"] >= MIN_SMM_COVERAGE
        and stats.loc[g, "pct_cells_smm_neoplastic"] >= 0.05]
print(f"after SMM coverage >={MIN_SMM_COVERAGE:.0%} of tumours and >=5% cells: {len(keep)} genes")
candidates = keep

# ------------------------------------------------ normal-tissue restriction
# Applied BEFORE the top-60 cut, so epitope prediction is spent on genes with a
# defensible normal-tissue profile rather than on whatever scored highest on DE.
_tiers = Path(f"{DATA}/tissue_tiers.csv")
if not _tiers.exists():
    raise SystemExit(f"missing {_tiers}\nRun 02b_tissue_restriction.py first.")
tissue = pd.read_csv(_tiers).set_index("gene")
print(f"Applying normal-tissue restriction tiers ({len(tissue)} genes tested)...")
tiers = tissue["tier"].reindex(candidates).fillna("no_hpa")
print("  tier breakdown of the candidate pool:")
for tname, cnt in tiers.value_counts().sort_index().items():
    print(f"    {tname:<14} {cnt:>4}")
candidates = [g for g in candidates if tiers[g] not in REJECT_TIERS]
print(f"after rejecting neural/vital-tissue outliers: {len(candidates)} genes")

# preliminary DE rank to cap the list for epitope prediction.
# The score rewards: big effect size, REPRODUCIBLY enriched within patients,
# expressed in many neoplastic cells, and NOT expressed in normal plasma cells.
# Note this uses recurrence, not detection-prevalence - the prevalence >=30%
# filter above is the separate "shared across patients" gate.
stats = stats.loc[[g for g in candidates if g in stats.index]]
# indication-matched: SMM tumour cells vs healthy-donor plasma cells
stats["specificity"] = (stats["pct_cells_smm_neoplastic"]
                        - stats["pct_cells_nbm_normal"]).clip(lower=0)
stats["de_score"] = (stats["log2FoldChange"].clip(lower=0)
                     * stats["recurrence_frac"]
                     * stats["specificity"])
candidates = stats.sort_values("de_score", ascending=False).head(60).index.tolist()
print(f"top 60 by preliminary DE score -> carrying forward")

# ---------------------------------------------------------- HPA annotation
print("Annotating from Human Protein Atlas...")
_hpa_path = hpa_tsv()
if not _hpa_path.exists():
    raise SystemExit(f"missing HPA table: {_hpa_path}\nRun 00_fetch_data.py first.")
hpa = pd.read_csv(_hpa_path, sep="\t")
hpa.columns = [c.strip('"') for c in hpa.columns]
hpa_idx = hpa.set_index("Gene")

rows = []
for gene in candidates:
    if gene in hpa_idx.index:
        r = hpa_idx.loc[gene]
        if isinstance(r, pd.DataFrame):  # duplicate rows -> take first
            r = r.iloc[0]
        rows.append({
            "gene": gene,
            "hpa_subcellular": r.get("Subcellular main location", ""),
            "hpa_protein_class": r.get("Protein class", ""),
            "hpa_secretome": r.get("Secretome location", ""),
            "hpa_uniprot_id": r.get("Uniprot", ""),
        })
    else:
        rows.append({"gene": gene, "hpa_uniprot_id": ""})

annot = pd.DataFrame(rows)

# safety score from the tissue tier (step 2b). Tier 1 - restricted to
# bone marrow/lymphoid, or to testis - scores 1.0; tier 2 (restricted, but to
# some other non-vital tissue) scores 0.6; genes absent from the consensus
# table are penalised at 0.3 rather than treated as neutral, because for a
# preventive vaccine "cannot assess" is not the same as "probably fine".
annot["tissue_tier"] = annot["gene"].map(tissue["tier"]).fillna("no_hpa")
annot["hpa_outlier_tissues"] = annot["gene"].map(tissue["sig_tissues"]).fillna("")
annot["hpa_lymphoid_fdr"] = annot["gene"].map(tissue["lymphoid_fdr"])
annot["hpa_safety_score"] = annot["tissue_tier"].map(TIER_SAFETY)

# ---------------------------------------------------------- localisation
# ANNOTATION ONLY - never filtered or scored on. HLA-I presents peptides from
# proteasomal degradation of cytosolic and nuclear proteins just as readily as
# membrane ones, so subcellular location is a hard requirement for antibodies,
# CAR-T and ADCs but NOT for a peptide vaccine. Five of the fifteen shortlisted
# genes are intracellular (CCND1, FCRLA, POU2F2, MLLT3, HEY2); a surface-only
# filter would delete the top TWO ranked genes.
#
# It is still worth reporting, because it says which OTHER modalities a target
# is open to - and because the reverse case is a genuine argument for this one.
# FCRLA (rank 2) is HPA-classified intracellular, which is exactly why the
# surface-targeting programmes in myeloma go after its family member FCRL5
# (cevostamab) instead. A peptide vaccine can use a target an antibody cannot.
#
# HPA's "Protein class" carries three mutually exclusive predictions, and they
# are the authoritative signal. Note that the "Secretome location" value
# "Intracellular and membrane" is a NEGATIVE - it marks a protein that is not
# secreted. Substring-matching the concatenated fields for "membrane" reads
# that as a positive and mislabels intracellular proteins as surface ones.
def classify_localisation(row):
    cls = str(row.get("hpa_protein_class", "") or "").lower()
    sec = str(row.get("hpa_secretome", "") or "").lower()
    sub = str(row.get("hpa_subcellular", "") or "").lower()

    if "predicted secreted proteins" in cls or sec.startswith("secreted"):
        return "secreted"
    if "predicted membrane proteins" in cls:
        return "membrane"
    if "predicted intracellular proteins" in cls:
        return "intracellular"
    # fall back to the IHC-based main location when no protein-class call exists
    if "plasma membrane" in sub:
        return "membrane"
    if sub and sub != "nan":
        return "intracellular"
    return "unknown"

annot["localisation"] = annot.apply(classify_localisation, axis=1)
# kept for backwards compatibility, now derived from the corrected call
annot["surface_or_secreted"] = annot["localisation"].isin(["membrane", "secreted"])

annot.to_csv(f"{DATA}/candidates_annotated.csv", index=False)
print(f"\nsaved {len(annot)} annotated candidates to candidates_annotated.csv")
print("\nPreview:")
print(annot.head(10).to_string(index=False))
