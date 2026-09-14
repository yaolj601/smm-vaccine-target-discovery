"""
Step 5: LLM-assisted scoring via the Anthropic API.

For each candidate gene, we send the bioinformatics evidence (DE stats,
epitope counts, normal-tissue restriction tier) to Claude and
ask for structured knowledge that traditional pipelines cannot get from this
dataset alone:
  - known role in myeloma biology
  - existing vaccine / antibody / CAR programs against it
  - safety signals in normal tissue
  - two categorical judgements - mm_confidence (is it plausibly a myeloma
    antigen?) and novelty_tier (is it still unclaimed?) - which are combined
    into a single `novelty_confidence` score by the lookup table in
    step 5b. Development tiers follow Yao et al., Cancer Res 2023 (Fig 3D),
    with tier 1 split by indication since a vaccine proven in another cancer
    is a very different proposition from one already in myeloma trials.

Responses are cached to JSON, so reruns cost nothing.

Requires ANTHROPIC_API_KEY (see README section 'LLM step').

Input:  $SMM_DATA_DIR/candidates_annotated.csv
        $SMM_DATA_DIR/epitope_results.csv
        $SMM_DATA_DIR/de_results.csv
Output: $SMM_DATA_DIR/llm_scores.json
"""

import json
import os
import time

import anthropic
import pandas as pd
from anthropic import Anthropic

from config import DATA, LLM_MODEL

MODEL = LLM_MODEL
CACHE_PATH = DATA / "llm_scores.json"

# Credentials come from the environment - never hard-code, never commit:
#
#     export ANTHROPIC_API_KEY="sk-ant-..."
#     export ANTHROPIC_WORKSPACE_ID="wrkspc_..."      # only if the key is
#                                                     # not workspace-scoped
#
# A key created OUTSIDE a workspace fails with:
#   "This API key is not scoped to a workspace, so this request must include
#    the anthropic-workspace-id header with the ID of the workspace to use."
# Either create the key inside the workspace in the Console, or set
# ANTHROPIC_WORKSPACE_ID and the header below supplies it per request.
if not (os.environ.get("ANTHROPIC_API_KEY") or os.environ.get("ANTHROPIC_AUTH_TOKEN")):
    print("NOTE: ANTHROPIC_API_KEY is not set; falling back to an `ant auth login` "
          "profile if one exists. See README section 'LLM step'.")

_workspace = os.environ.get("ANTHROPIC_WORKSPACE_ID")
_headers = {"anthropic-workspace-id": _workspace} if _workspace else None
if _workspace:
    print(f"using workspace {_workspace}")
client = Anthropic(default_headers=_headers)

# ------------------------------------------------- build the evidence table
annot = pd.read_csv(f"{DATA}/candidates_annotated.csv")
epi = pd.read_csv(f"{DATA}/epitope_results.csv")
de = pd.read_csv(f"{DATA}/de_results.csv")

# DE stats from the SMM-vs-healthy contrast (the vaccine indication)
de_b = de[de["contrast_name"] == "B_SMM_vs_NBM"].set_index("gene")

merged = annot.merge(epi, on="gene", how="inner")

# bring the per-cell expression stats over from the DE table
# (they are identical across contrasts, so take them from contrast A)
de_a = de[de["contrast_name"] == "A_neoplastic_vs_normal"].set_index("gene")
for col in ["pct_cells_neoplastic", "pct_cells_normal", "patient_prevalence"]:
    merged[col] = merged["gene"].map(de_a[col])
print(f"{len(merged)} candidates with full evidence")

# ------------------------------------------------- the prompt
SYSTEM_PROMPT = """You are an immunologist advising an off-the-shelf preventive
vaccine program for smouldering multiple myeloma (SMM). The goal is peptide/HLA-I
vaccine targets: tumor-associated antigens expressed by neoplastic plasma cells,
presented on HLA-I, and safe to break tolerance against.

You will receive computational evidence from a single-cell RNA-seq dataset
(neoplastic vs normal plasma cells across MGUS/SMM/MM) plus immunoinformatics.
Judge each candidate on:
1. Biological plausibility as an MM/SMM antigen (is the gene known to be
   expressed/dysfunctional in myeloma? is it a cancer-testis or oncofetal antigen?)
2. Existing development (vaccines, antibody-drug conjugates, CAR-T, bispecifics -
   existing programs validate immunogenicity and feasibility)
3. Safety risk (broad normal-tissue expression, essential housekeeping function,
   risk of autoimmune attack on normal plasma cells or other tissues)

Return TWO categorical judgements rather than a numeric score. They are
combined downstream by a published lookup grid (step 5b), so a
category is auditable in a way a 0-10 number is not.

CONFIDENCE - is this plausibly a myeloma antigen at all?
  "A"  established myeloma antigen or oncogene
  "B"  plausible - related pathway, family member, or lineage rationale
  "C"  no known myeloma role
  "D"  known housekeeping or broadly essential gene

D is disqualifying: a gene marked D is dropped from the shortlist regardless
of how good its expression data looks. Use it for ribosomal proteins, tubulins,
core metabolic and mitochondrial machinery, and the like.

NOVELTY - how far has this target already been developed? After Yao et al.
(Cancer Res 2023, Fig 3D), with tier 1 split by indication:

  "1a_vaccine_elsewhere"  vaccine already in clinical study in ANOTHER cancer
  "1b_vaccine_in_mm"      vaccine already in clinical study IN myeloma
                          (validated, but crowded)
  "2_clinical_other"      clinical in myeloma as CAR / antibody / ADC / bispecific
  "3_literature"          published myeloma relevance, no clinical program
  "4_uncharacterised"     not previously described as a myeloma marker

The 1a/1b distinction matters: a vaccine proven elsewhere but unclaimed in
myeloma is the most attractive position there is, while one already in myeloma
trials is the least. If you assert 1a, 1b or 2, name the agent or trial.

"4_uncharacterised" is the correct answer for a gene you do not recognise. Do
not invent a program to justify a higher tier; an unfounded clinical claim is
far more damaging than an honest 4.

Be honest about uncertainty. If you do not know a gene, say so."""

def make_user_message(row):
    """One compact evidence block per gene."""
    return f"""Gene: {row['gene']}
UniProt: {row['hpa_uniprot_id']}

Computational evidence from the SMM scRNA-seq pipeline:
- SMM neoplastic vs healthy plasma cells: log2FC={row_de(row, 'log2FoldChange')}, padj={row_de(row, 'padj')}
- Expressed in {row['pct_cells_neoplastic']:.0%} of neoplastic cells vs {row['pct_cells_normal']:.0%} of normal plasma cells
- Detected in {row['patient_prevalence']:.0%} of patients with neoplastic cells
- Normal-tissue restriction (HPA consensus outlier test): tier {row['tissue_tier']}; significant outlier in: {row['hpa_outlier_tissues'] or 'no tissue'}
- HPA subcellular location: {row['hpa_subcellular']}
- HLA-I epitopes: {row['n_strong_binders']} strong binders (<50 nM) across {row['n_alleles_with_binder']}/27 common alleles; best affinity {row['best_affinity_nM']:.0f} nM

Respond with ONLY a JSON object:
{{"mm_confidence": "A|B|C|D",
  "novelty_tier": "1a_vaccine_elsewhere|1b_vaccine_in_mm|2_clinical_other|3_literature|4_uncharacterised",
  "existing_programs": "name the specific agent/trial if you claim 1a, 1b or 2, else 'none known'",
  "safety_concerns": "...", "rationale": "2-3 sentences"}}"""


def row_de(row, col):
    """Look up a DE stat for this gene (NaN-safe)."""
    if row["gene"] in de_b.index:
        val = de_b.loc[row["gene"], col]
        return f"{val:.2f}" if pd.notna(val) else "n/a"
    return "n/a"


# ------------------------------------------------- call the API (with cache)
cache = {}
if os.path.exists(CACHE_PATH):
    with open(CACHE_PATH) as f:
        cache = json.load(f)
    print(f"loaded {len(cache)} cached LLM scores")

for i, row in merged.iterrows():
    gene = row["gene"]
    if gene in cache:
        continue
    print(f"querying LLM for {gene} ({i + 1}/{len(merged)})...")
    try:
        resp = client.messages.create(
            model=MODEL,
            max_tokens=4000,
            system=SYSTEM_PROMPT,
            messages=[{"role": "user", "content": make_user_message(row)}],
        )
        # response.content is a list of blocks. Current models return a
        # ThinkingBlock first, so content[0].text raises AttributeError -
        # select the text block by type rather than by position.
        text = "".join(b.text for b in resp.content if b.type == "text")
        if resp.stop_reason == "max_tokens":
            raise ValueError(f"response truncated at max_tokens for {gene}")
        # strip markdown code fences if the model added them
        text = text.strip().removeprefix("```json").removeprefix("```").removesuffix("```").strip()
        cache[gene] = json.loads(text)
        with open(CACHE_PATH, "w") as f:
            json.dump(cache, f, indent=2)
    except TypeError as e:
        # The SDK raises a bare TypeError (not an anthropic.* error) when it
        # cannot resolve ANY credential - the most likely failure for someone
        # running this fresh. Translate it into something actionable.
        raise SystemExit(
            f"\nFATAL: no Anthropic credentials found.\n  {e}\n\n"
            "Set a key created inside your Anthropic workspace:\n"
            "    export ANTHROPIC_API_KEY=\"sk-ant-...\"\n\n"
            "Or skip this step - step 6 falls back to the provenance-labelled "
            "in-session scores in data/llm_scores_insession.json."
        )
    except (anthropic.AuthenticationError, anthropic.PermissionDeniedError,
            anthropic.BadRequestError) as e:
        # Credential / request-shape problems affect every gene equally - stop
        # now rather than writing 43 identical error records to the cache.
        raise SystemExit(
            f"\nFATAL: the API rejected the request for {gene}:\n  {e}\n\n"
            "If this mentions workspace scoping, the key was created outside the "
            "workspace this project runs in. Either create a key inside that "
            "workspace in the Console, or set ANTHROPIC_WORKSPACE_ID to it and "
            "re-export ANTHROPIC_API_KEY."
        )
    except (anthropic.RateLimitError, anthropic.APIStatusError,
            anthropic.APIConnectionError) as e:
        # Transient / per-gene failures: record and carry on. Re-running the
        # script retries only the genes that are missing from the cache.
        print(f"  ERROR for {gene}: {e}")
        cache[gene] = {"error": str(e)}
        with open(CACHE_PATH, "w") as f:
            json.dump(cache, f, indent=2)
    time.sleep(1)  # be polite to the API

print(f"\nLLM scoring complete: {len(cache)} genes scored")
