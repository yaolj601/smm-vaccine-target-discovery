#!/bin/bash
# Run the full SMM vaccine-target pipeline end-to-end.
#
#   pip install -r ../requirements.txt
#   mhcflurry-downloads fetch models_class1_pan models_class1_presentation
#   export ANTHROPIC_API_KEY="sk-ant-..."   # optional, see README 'LLM step'
#   ./run_all.sh
#
# Paths default to <repo>/{data,reference,results} and can be overridden with
# SMM_DATA_DIR / SMM_REF_DIR / SMM_OUT_DIR. See ../README.md.
set -euo pipefail
cd "$(dirname "$0")"

python 00_fetch_data.py           # GEO + HPA download (both automatic)
python 01_preprocess.py
python 01b_cohort_overview.py   # cohort structure figure
python 02_differential.py
python 02b_tissue_restriction.py  # HPA normal-tissue outlier test
python 03_specificity_safety.py
python 04_epitopes.py

# Step 5 needs a workspace-scoped Anthropic key. If it is missing or the API
# rejects the call, fall through to the provenance-labelled in-session scores
# rather than failing the whole run - step 6 prefers API scores when present.
if ! python 05_llm_scoring.py; then
  echo "WARNING: LLM API step did not complete; step 5b will fall back to"
  echo "         llm_scores_insession.json (see README 'LLM step')."
fi

python 05b_novelty_confidence.py   # combine confidence x novelty into one score
python 06_rank_and_report.py
python 07_slides.py                # build the .pptx deck from results/

echo
echo "Done. Outputs in ${SMM_OUT_DIR:-$(cd .. && pwd)/results}/"
