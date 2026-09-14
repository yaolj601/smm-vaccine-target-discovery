"""
Step 5b: novelty_confidence - how attractive is this antigen to *this* program?

Combines the two categorical judgements step 5 returns into one 0-10 score:

  CONFIDENCE - is this plausibly a myeloma antigen at all?  (A / B / C / D)
  NOVELTY    - is the target still available to us? After Yao et al., Cancer
               Res 2023 (Fig 3D), with tier 1 split by indication:
                 1a  vaccine already in clinical study in ANOTHER cancer
                 1b  vaccine already in clinical study IN myeloma (crowded)
                 2   clinical in myeloma as CAR / antibody / ADC / bispecific
                 3   published myeloma relevance, no clinical program
                 4   not previously described as a myeloma marker

They are combined by a lookup over the grid rather than by a weighted sum,
because the two axes do not point the same way. Confidence is monotonic - A
always beats D. Novelty is not: an antigen already in a myeloma vaccine trial
is validated but taken. Summing two monotonic axes can never produce the
inversion the grid does - DKK1 (A, already a vaccine in myeloma) scores 3 while
TBXAS1 (C, unclaimed) scores 4 - and that inversion is the whole point of
tracking novelty.

This replaces an earlier `literature_score`: a 0-10 integer with no formula
behind it, carrying a quarter of the composite, and largely redundant once the
tissue tier existed (Spearman -0.65 against the development tier, -0.40 against
the tissue tier). A category is auditable; that number was not.

Input:  $SMM_DATA_DIR/llm_scores.json  (or llm_scores_insession.json)
Output: $SMM_DATA_DIR/novelty_confidence.csv
"""

import json
import os

import pandas as pd

from config import DATA

CONFIDENCE_LEVELS = {
    "A": "established myeloma antigen or oncogene",
    "B": "plausible - related pathway, family or lineage rationale",
    "C": "no known myeloma role",
    "D": "known housekeeping / broadly essential",
}

# Development-status tiers, after Yao et al., Cancer Res 2023 (Fig 3D), with
# tier 1 split by indication because that distinction is the one that matters
# most for a new program: a vaccine proven in another cancer is de-risked AND
# still unclaimed here, while one already in myeloma trials is simply taken.
NOVELTY_TIERS = {
    "1a_vaccine_elsewhere": "clinical vaccine program in another cancer",
    "1b_vaccine_in_mm": "clinical vaccine program in myeloma (crowded)",
    "2_clinical_other": "clinical in myeloma as CAR / antibody / ADC / bispecific",
    "3_literature": "published myeloma relevance, no clinical program",
    "4_uncharacterised": "not previously described as a myeloma marker",
}

# confidence x novelty -> 0-10 "attractiveness to this program".
#
# The novelty axis is deliberately NON-monotonic in development status. An
# antigen with published myeloma relevance and no clinical program (tier 3)
# scores ABOVE one already in myeloma vaccine trials (tier 1b): in trials it is
# validated but taken, and the value for a new program sits in
# validated-but-unclaimed. A vaccine proven in ANOTHER cancer (1a) is the best
# cell on the board - de-risked modality, open indication. No weighted sum of
# two monotonic axes can produce that inversion, which is why this is a lookup.
#
# None marks a self-contradictory pairing (an established antigen cannot be
# uncharacterised); `score()` falls back to the nearest defensible cell.
GRID = {
    "A": {"1a_vaccine_elsewhere": 10, "1b_vaccine_in_mm": 3, "2_clinical_other": 8,
          "3_literature": 9, "4_uncharacterised": None},
    "B": {"1a_vaccine_elsewhere": 8, "1b_vaccine_in_mm": 2, "2_clinical_other": 7,
          "3_literature": 7, "4_uncharacterised": 6},
    "C": {"1a_vaccine_elsewhere": 5, "1b_vaccine_in_mm": None, "2_clinical_other": 5,
          "3_literature": 4, "4_uncharacterised": 3},
    # D is disqualifying whatever its development status. A known housekeeping
    # gene does not belong on a vaccine shortlist at any rank, however good its
    # expression data looks - the same logic that rejects vital-tissue outliers.
    "D": {k: 0 for k in ("1a_vaccine_elsewhere", "1b_vaccine_in_mm",
                         "2_clinical_other", "3_literature", "4_uncharacterised")},
}

REJECT_CONFIDENCE = {"D"}

# legacy tier string -> current one, so older cache files still resolve
ALIASES = {"1_clinical_vaccine": "1b_vaccine_in_mm",
           "2_clinical_antibody": "2_clinical_other"}


def score(confidence, novelty):
    """Look up the 0-10 score for a (confidence, novelty) pair."""
    novelty = ALIASES.get(novelty, novelty)
    if confidence not in GRID or novelty not in NOVELTY_TIERS:
        return None
    v = GRID[confidence][novelty]
    if v is None:
        # contradictory pairing: treat as the adjacent, defensible cell
        v = GRID[confidence]["3_literature"]
    return v


def rejected(confidence):
    return confidence in REJECT_CONFIDENCE


if __name__ == "__main__":
    # Precedence: API grades beat the in-session fallback.
    #
    # MIND THE ORDER. dict.update() lets whichever source loads LAST win, so
    # the fallback must be loaded FIRST and the API grades second. Writing the
    # list the other way round reads like a preference order ("use the API
    # file, else the fallback") but does the opposite - and that was a real
    # bug: with llm_scores.json listed first, every gene the two files
    # disagreed on silently took the in-session grade. It affected 6 of 32
    # genes and left SCYL2 (API grade D = reject) sitting at rank 19.
    llm, provenance = {}, {}
    for source, path in [("in-session", f"{DATA}/llm_scores_insession.json"),
                         ("api", f"{DATA}/llm_scores.json")]:
        if os.path.exists(path):
            with open(path) as f:
                loaded = json.load(f)
            llm.update(loaded)
            provenance.update({g: source for g in loaded})
            print(f"  loaded {len(loaded):>3} grades from "
                  f"{os.path.basename(path)}  [{source}]")
    if llm:
        n_api = sum(v == "api" for v in provenance.values())
        print(f"  effective provenance: {n_api} api, "
              f"{len(provenance) - n_api} in-session fallback")
    # Verified overrides: primary-source checks that supersede a model grade.
    # llm_scores.json is never edited - it is the evidence trail of what the
    # model said, and rewriting it would destroy the record. Corrections live
    # in verified_overrides.json, each carrying its source, and are logged
    # here so a changed grade is never silent.
    ov_path = f"{DATA}/verified_overrides.json"
    if os.path.exists(ov_path):
        with open(ov_path) as f:
            overrides = json.load(f)
        for gene, o in overrides.items():
            if gene.startswith("_") or gene not in llm:
                continue
            for field in ("mm_confidence", "novelty_tier"):
                if field in o and llm[gene].get(field) != o[field]:
                    print(f"  OVERRIDE  {gene}.{field}: "
                          f"{llm[gene].get(field)} -> {o[field]}")
                    llm[gene][field] = o[field]
                    provenance[gene] = "verified"

    if not llm:
        raise SystemExit(f"no LLM scores found in {DATA}. Run 05_llm_scoring.py, "
                         "or provide llm_scores_insession.json.")

    rows = []
    for gene, v in llm.items():
        if not isinstance(v, dict):
            continue
        conf, nov = v.get("mm_confidence"), v.get("novelty_tier")
        rows.append({"gene": gene, "mm_confidence": conf, "novelty_tier": nov,
                     "novelty_confidence": score(conf, nov),
                     "rejected": rejected(conf),
                     "existing_programs": v.get("existing_programs", ""),
                     "llm_source": provenance.get(gene, "unknown")})
    out = pd.DataFrame(rows).set_index("gene").sort_index()
    out.to_csv(f"{DATA}/novelty_confidence.csv")

    print(f"scored {len(out)} genes")
    print("\n  grid occupancy (confidence x novelty):")
    print(pd.crosstab(out["mm_confidence"], out["novelty_tier"]).to_string())
    print(f"\n  confidence D -> rejected downstream: {int(out['rejected'].sum())}")
    print("  grades by source: "
          + out["llm_source"].value_counts().to_dict().__repr__())
    print(f"\nSaved {DATA}/novelty_confidence.csv")
