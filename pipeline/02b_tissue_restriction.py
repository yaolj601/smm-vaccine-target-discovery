"""
Step 2b: Normal-tissue restriction by outlier testing (HPA consensus RNA).

Replaces an earlier five-bucket mapping of HPA's `RNA tissue specificity`
string, which asked only "is this gene restricted to few tissues?" and never
"restricted to WHICH tissues?". That let CADM1 score a perfect 1.0 while being
retina-enriched, and missed retina, choroid plexus, hypothalamus, cerebral
cortex, cerebellum, liver, kidney and adrenal entirely, because the vital-tissue
cap only substring-matched brain / heart / skeletal muscle / spinal cord.

Method (after Yao et al., github.com/Lijunyao8/tumor_specific_marker):
for each gene, march through its ~50 consensus tissue values; hold one out,
compute the mean and variance of the rest, and run a one-sample t-test asking
whether the held-out tissue is a HIGH outlier against that background. Correct
across every gene x tissue test with the same BH-style procedure the original
Perl uses (p * N / rank, capped at 1, forced monotonic) and cut at FDR <= 0.05.
A gene's "outlier tissues" are the ones that survive.

Tiers, adapted for a peptide/HLA-I vaccine rather than a surface CAR/ADC target:

  1a_lymphoid   outlier in bone marrow / thymus / spleen / appendix / tonsil /
                lymph node. On-target damage lands in the haematopoietic
                compartment, which is monitorable and survivable.

                All six are pooled, as in the source method - do NOT be tempted
                to narrow this to bone marrow on the reasoning that myeloma
                lives there. HPA's bone-marrow sample is whole marrow, where
                plasma cells are ~1% of cells and granulopoietic/erythroid
                precursors dominate, so "bone marrow outlier" is a MYELOID
                signature. Tested against known markers:
                  myeloid/erythroid  MPO, ELANE, PRTN3, CTSG, GYPA, AHSP,
                                     ITGAM, DEFA4    -> bone marrow, 8/8
                  B cell             CD19, MS4A1, CD79A/B, PAX5, BANK1,
                                     AICDA           -> lymph node/tonsil, 0/7
                  plasma cell        BCMA, SLAMF7, FCRL5, MZB1
                                                     -> tonsil/LN/spleen, 0/4
                The three most clinically validated myeloma surface targets -
                BCMA, SLAMF7, FCRL5 - are all tonsil or lymph-node outliers and
                none is a bone-marrow outlier. Filtering to marrow would promote
                MPO-like genes, where off-target damage means marrow failure,
                and demote the B/plasma compartment myeloma actually belongs to.
  1b_testis_CT  outlier in testis only - the cancer-testis antigen class. Germ
                cells are immune-privileged and do not present HLA-I, which is
                why NY-ESO-1 and MAGE-A are the gold standard for cancer
                vaccines. A surface-marker pipeline rejects testis; a vaccine
                pipeline should not.
  2_restricted  restricted, but to some other non-vital tissue.
  3_vital       outlier in neural or vital tissue. REJECTED in step 3.
  4_broad       no significant outlier anywhere - broadly expressed.

This runs once and caches to disk; step 3 reads the table rather than
recomputing a ~1M-test correction on every invocation.

Input:  $SMM_REF_DIR/rna_tissue_consensus.tsv
Output: $SMM_DATA_DIR/tissue_tiers.csv
"""

import numpy as np
import pandas as pd
from scipy import stats

from config import DATA, HPA_CONSENSUS

LYMPHOID = {"bone marrow", "thymus", "spleen", "appendix", "tonsil", "lymph node"}
IMMUNE_PRIVILEGED = {"testis"}
NEURAL = {"retina", "choroid plexus", "hypothalamus", "cerebral cortex", "cerebellum",
          "amygdala", "basal ganglia", "midbrain", "hippocampal formation",
          "spinal cord", "pituitary gland"}
VITAL = {"heart muscle", "skeletal muscle", "liver", "kidney", "lung", "pancreas",
         "adrenal gland"}

FDR_CUTOFF = 0.05

# Safety score per tier. Tier 3 is gated out upstream; the value is kept so a
# gene that somehow reaches scoring cannot benefit from it.
TIER_SAFETY = {
    "1a_lymphoid": 1.00,
    "1b_testis_CT": 1.00,
    "2_restricted": 0.60,
    "4_broad": 0.20,
    "no_hpa": 0.30,   # unknown is penalised, not treated as neutral
    "3_vital": 0.00,
}
REJECT_TIERS = {"3_vital"}


def outlier_tissues(consensus_tsv):
    """Return a DataFrame indexed by gene with the significant outlier tissues.

    Columns: sig_tissues (comma-joined), n_sig_tissues, lymphoid_fdr, tier.
    """
    df = pd.read_csv(consensus_tsv, sep="\t")
    mat = df.pivot_table(index="Gene name", columns="Tissue", values="nTPM",
                         aggfunc="max").dropna(axis=0, how="any")
    X = mat.values.astype(float)
    tissues = np.array(mat.columns)
    n = X.shape[1]
    m = n - 1                                   # background size after hold-out

    # leave-one-out mean/variance without building n copies of the matrix
    tot = X.sum(1, keepdims=True)
    tot2 = (X ** 2).sum(1, keepdims=True)
    mean_b = (tot - X) / m
    var_b = np.clip(((tot2 - X ** 2) - m * mean_b ** 2) / (m - 1), 1e-12, None)
    t = (X - mean_b) / np.sqrt(var_b * (m + 1) / m)
    p = stats.t.sf(t, df=m - 1)                 # one-sided: HIGH outliers only

    # BH-style correction over every gene x tissue test, as in the original
    flat = p.ravel()
    order = np.argsort(flat, kind="mergesort")
    N = flat.size
    fdr_sorted = np.maximum.accumulate(
        np.minimum(flat[order] * N / np.arange(1, N + 1), 1.0))
    fdr = np.empty(N)
    fdr[order] = fdr_sorted
    fdr = fdr.reshape(p.shape)

    sig = fdr <= FDR_CUTOFF
    lym = np.array([x in LYMPHOID for x in tissues])
    lym_fdr = np.where(sig[:, lym].any(1),
                       np.where(sig[:, lym], fdr[:, lym], np.inf).min(1), np.nan)

    out = pd.DataFrame({
        "sig_tissues": [", ".join(tissues[s]) for s in sig],
        "n_sig_tissues": sig.sum(1),
        "lymphoid_fdr": lym_fdr,
    }, index=mat.index)
    out.index.name = "gene"
    out["tier"] = [classify(s) for s in out["sig_tissues"]]
    return out


def classify(sig_tissues):
    """Assign a vaccine-relevant tier from a gene's outlier tissue set."""
    t = {x.strip() for x in str(sig_tissues).split(",") if x.strip()}
    if not t:
        return "4_broad"
    if t & (NEURAL | VITAL):
        return "3_vital"
    if t & LYMPHOID:
        return "1a_lymphoid"
    if t <= IMMUNE_PRIVILEGED:
        return "1b_testis_CT"
    return "2_restricted"


if __name__ == "__main__":
    if not HPA_CONSENSUS.exists():
        raise SystemExit(f"missing HPA consensus table: {HPA_CONSENSUS}\n"
                         "Run 00_fetch_data.py first.")
    print("Outlier-testing normal tissue expression (HPA consensus)...")
    out = outlier_tissues(HPA_CONSENSUS)
    out.to_csv(f"{DATA}/tissue_tiers.csv")
    print(f"  {len(out)} genes tested across the consensus tissue panel")
    print("\n  tier distribution, genome-wide:")
    for t, n in out["tier"].value_counts().sort_index().items():
        print(f"    {t:<14} {n:>6}")
    print(f"\nSaved {DATA}/tissue_tiers.csv")
