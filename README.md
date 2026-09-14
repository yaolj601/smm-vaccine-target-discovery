# Ranked vaccine-target shortlist for smouldering multiple myeloma (SMM)

A six-step pipeline that takes the GSE193531 single-cell RNA-seq dataset
(MGUS -> SMM -> MM plus healthy bone marrow) and produces a ranked, evidence-
backed shortlist of candidate peptide/HLA-I vaccine targets.

The reasoning behind the scoring criteria, what the LLM step added, next steps
and caveats are in **[WRITEUP.md](WRITEUP.md)**. This file is just how to run it.

## Quickstart

```bash
python -m venv .venv && source .venv/bin/activate     # Python >= 3.10
pip install -r requirements.txt
mhcflurry-downloads fetch models_class1_pan models_class1_presentation

export ANTHROPIC_API_KEY="sk-ant-..."                 # see 'LLM step'
export ANTHROPIC_WORKSPACE_ID="wrkspc_..."             # only if key is not workspace-scoped

cd pipeline && ./run_all.sh
```

Runtime is dominated by two steps: loading the count matrix (step 1, ~1-2 min)
and MHCflurry prediction over 38 proteins x 27 alleles (step 4, tens of minutes).

## Data

`00_fetch_data.py` downloads the two public inputs automatically:

| Source | File | How |
|---|---|---|
| [GSE193531](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193531) | `GSE193531_umi-count-matrix.csv.gz` -> `data/counts.csv.gz` | automatic |
| GSE193531 | `GSE193531_cell-level-metadata.csv.gz` -> `data/meta.csv.gz` | automatic |
| [Human Protein Atlas](https://www.proteinatlas.org/about/download) | `proteinatlas.tsv` -> `reference/` | automatic |
| Human Protein Atlas | `rna_tissue_consensus.tsv` -> `reference/` (per-tissue nTPM, needed for the outlier test) | automatic |

Both inputs download unattended, so `run_all.sh` works from a clean clone with
no manual steps. DepMap CRISPR essentiality was evaluated and deliberately
dropped - rationale in [WRITEUP.md](WRITEUP.md) section 1.

## Layout and path overrides

```
pipeline/        00_fetch_data.py ... 06_rank_and_report.py, config.py, run_all.sh
data/            downloaded inputs + intermediates (created on first run)
reference/       Human Protein Atlas table (created on first run)
results/         ranked_targets.{csv,md}, all_candidates_scored.csv,
                 evidence/, figures/
```

Nothing is hard-coded to an absolute path. All of it lives in
[pipeline/config.py](pipeline/config.py) and is overridable:

| Variable | Default | Purpose |
|---|---|---|
| `SMM_DATA_DIR` | `<repo>/data` | inputs + intermediates |
| `SMM_REF_DIR` | `<repo>/reference` | Human Protein Atlas table |
| `SMM_OUT_DIR` | `<repo>/results` | final deliverables |
| `SMM_HPA_TSV` | `<ref>/proteinatlas.tsv` | explicit HPA annotation path |
| `SMM_HPA_CONSENSUS` | `<ref>/rna_tissue_consensus.tsv` | explicit consensus-tissue path |
| `SMM_LLM_MODEL` | `claude-opus-5` | model used in step 5 |
| `ANTHROPIC_API_KEY` | — | required for step 5 |
| `ANTHROPIC_WORKSPACE_ID` | — | sent as the `anthropic-workspace-id` header when the key is not workspace-scoped |

## LLM step

Step 5 sends each candidate's evidence block to Claude and asks for structured
JSON: known MM relevance, existing vaccine/antibody/CAR programs, safety
signals, and a 0-10 literature score. Responses are cached to
`data/llm_scores.json`, so re-runs are free and only missing genes are queried.

The key is read from `ANTHROPIC_API_KEY` by the SDK. If it was created *outside*
a workspace the API returns:

> *"This API key is not scoped to a workspace, so this request must include the
> `anthropic-workspace-id` header with the ID of the workspace to use. Add the
> header, or use an API key that is scoped to a workspace."*

Either remedy works. To use the header, export the workspace id:

```bash
export ANTHROPIC_WORKSPACE_ID="wrkspc_..."   # your own workspace ID
```

If the step cannot run, `run_all.sh` continues and step 5b falls back to
`llm_scores_insession.json` - assessments made in-session against the identical
rubric, labelled as such by the `llm_available` column in the output. Cost for
a full run is a few cents (43 short calls).

## Outputs

| File | What |
|---|---|
| `results/ranked_targets.csv` / `.md` | the top-15 shortlist with per-gene rationale |
| `results/all_candidates_scored.csv` | all 43 scored candidates |
| `results/figures/` | composite score breakdown, specificity vs HLA coverage, evidence heatmap |
| `results/evidence/` | intermediate tables (DE results, annotations, Wilcoxon sensitivity, LLM scores) so every number in the shortlist can be traced back |

## Pipeline steps

| Step | Does |
|---|---|
| `00_fetch_data.py` | download GEO + HPA |
| `01_preprocess.py` | align counts/metadata, recover SMM risk stratum, light QC, write AnnData |
| `01b_cohort_overview.py` | cohort structure figure |
| `02_differential.py` | pseudobulk DESeq2 across 3 contrasts, within-patient recurrence, per-stage coverage |
| `02b_tissue_restriction.py` | HPA consensus outlier test → normal-tissue tiers |
| `03_specificity_safety.py` | apply tissue tiers, reject neural/vital outliers, HPA protein annotation |
| `04_epitopes.py` | MHCflurry presentation over 27 common HLA-A/B alleles |
| `05_llm_scoring.py` | Claude returns two categoricals: confidence (A–D) and novelty tier |
| `05b_novelty_confidence.py` | combine the pair into one 0–10 score via the lookup grid |
| `06_rank_and_report.py` | weighted composite, shortlist, figures, evidence trail |

Composite weights live in one place, `WEIGHTS` in
[pipeline/06_rank_and_report.py](pipeline/06_rank_and_report.py) - the formula
printed into `ranked_targets.md` is derived from it, so the documented weights
cannot drift from the applied ones.
