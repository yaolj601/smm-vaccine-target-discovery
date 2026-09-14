# Vaccine Target Shortlist for Smouldering Multiple Myeloma

Ranked peptide/HLA-I targets for an off-the-shelf preventive vaccine in SMM,
from GSE193531 scRNA-seq (35 patients, 29,387 CD138+ plasma cells, NBM through
MM). 22,273 genes → 32 scored → 15 shortlisted. Every number traces back
through `results/evidence/`.

**Shortlist, in order:** CCND1, FCRLA, NLGN4X, CD1D, KIT, SPN, POU2F2, MLLT3, IGF1, MYEOV, TBXAS1, IFITM1, HGF, HERC5, HEY2.

## 1. Scoring criteria, and why

SMM patients are not yet sick and many never progress, so the product sets the
bar. Four layers, weighted by how much I trust each:

| Layer | Weight | Why that weight |
|---|---|---|
| Differential expression (pseudobulk DESeq2) | 0.30 | Most grounded. Effect size × recurrence × tumour specificity |
| Normal-tissue restriction (HPA outlier test) | 0.25 | A preventive vaccine cannot tolerate vital-tissue targets |
| Knowledge layer (LLM-derived) | 0.25 | Crowding, mechanism and known toxicity are not in the data |
| HLA-I presentation (MHCflurry, 27 alleles) | 0.20 | Nothing fails it — 38/38 carry a sub-50 nM binder — and peptide choice stays revisable after target selection |

Four choices worth defending:

- **Three contrasts, because "up in tumour" is ambiguous until you say
  *compared with whose normal cells*.** A pools all tumour vs all normal; B is
  12 SMM tumours vs 9 healthy donors — the clinical question, and the effect
  size used; C is the 11 paired patients against themselves. Entry is a union,
  but pooled-signal admissions must still be up in SMM.
- **Specificity is measured on the indication.** The pooled neoplastic set is
  57% MM cells, so apparent tumour specificity can be an MM phenomenon. The
  term uses SMM tumour cells against healthy-donor plasma cells.
- **Tissue restriction asks *which* tissue, not *how many*.** A leave-one-out
  outlier test over ~50 HPA tissues, FDR ≤ 0.05. Lymphoid and testis score
  1.00; neural or vital-organ outliers are **rejected**. It moves the shortlist
  most, cutting nine otherwise-ranking genes.
- **Knowledge is a 4×5 grid, not a weighted sum.** Plausibility is monotonic;
  novelty is not — an antigen already in a myeloma vaccine trial is validated
  *but taken*, and no sum of two monotonic axes produces that inversion.

**The largest assumption:** all seven clinically validated myeloma targets
(BCMA, GPRC5D, CD38, SLAMF7, FCRL5, CD138, XBP1) are plasma-cell *lineage*
antigens, and all seven fail this pipeline's DE gate. For high-risk SMM — 9 of
the 12 SMM samples here — eliminating the lineage compartment is the strategy
with clinical precedent; tumour-versus-normal selection has none.

## 2. Where the LLM was used, and what it added

The model returns two **categoricals** — plausibility (A–D) and development
status (1a–4) — never a score. A fixed grid converts the pair to 0–10, so the
scoring policy lives in twenty numbers you can argue with rather than inside a
model's answer.

DESeq2 + HPA + MHCflurry alone give a ranked list with no notion of whether a
target is already taken, no mechanism, and no prior on known toxicity. The LLM
added four things:

- **Existing clinical programmes with trial identifiers** — checkable, unlike a
  recollection.
- **A mechanism the data could not supply.** NLGN4X ranks 3 because its
  Y-paralogue is an H-Y minor histocompatibility antigen driving
  graft-versus-myeloma — unimplied by any expression value.
- **Safety signals invisible to the dataset** — interferon-stimulated genes
  look tumour-specific here but are induced in all cells during inflammation.
- **Rejection of housekeeping genes** strong DE alone would have promoted.

**AI assistance with the code, as distinct from the pipeline's LLM step.**
Everything above concerns step 05, where the Claude API grades candidate genes —
an LLM used as a pipeline component. Separately, the code itself was written
with Claude Code (Claude Opus 5): the scripts in `pipeline/`, the figures, and
these documents are AI-authored under my direction. The methods are mine:
**three of the four scoring layers adapt my own published approach** — tissue
restriction by leave-one-out outlier test, within-patient recurrence with its
0.25 log2FC floor, and the novelty tiering, all from [Yao et al., *Cancer Res*
2023](https://pubmed.ncbi.nlm.nih.gov/36779841/) and
[`tumor_specific_marker`](https://github.com/Lijunyao8/tumor_specific_marker).
I also rejected scoring by adjusted p-value — significance qualifies a target,
it does not rank it — and caught two things reported to me as done that were
not: tumour specificity was still stage-pooled, and the LLM step had never run
via the API. The assistant caught defects I would not have — a gene filter silently
discarding 38 non-immunoglobulin genes, a dict-precedence bug that made the
whole API grading run a no-op, and candidates lost to network timeouts. The
failure modes differ: mine were judgement calls about method, its were silent
defects in working code.

## 3. What I would do next

1. **Replicate on an independent SMM/MM cohort.** The shortlist rests on one
   dataset of 35 patients, 12 of them SMM. FCRLA is the one to watch — it has
   the tightest lymphoid restriction in the set and the same tissue profile as
   FCRL5, already in the clinic.
2. **Tooling:** cross-check top epitopes against netMHCpan-4.1, which disagrees
   with MHCflurry on 10–20% of peptides; score allele coverage against a
   length-matched null, since it still correlates 0.89 with protein length; and
   add the HPA protein (IHC) layer, the half of the tissue method not yet
   implemented.

## 4. Caveats

- **The core contrast may be wrong for the enrolment population** (§1 above).
  This is an assumption, not a finding.
- **Everything is RNA**, including the tissue tiers; where HPA protein data
  exists it mostly agrees, but CPQ stains highest in kidney and liver. And
  recurrence rests on the 11 of 35 samples that are paired.
