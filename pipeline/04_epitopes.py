"""
Step 4: Epitope-level assessment with MHCflurry.

WHY THIS STEP EXISTS
--------------------
Steps 2-3 asked "is this gene turned on in tumour and off in normal tissue?"
That is necessary but not sufficient. A CD8 T-cell never sees a gene, or even
a protein - it sees a 8-11 aa peptide sitting in the groove of an HLA class I
molecule on the cell surface. A gene can be beautifully tumour-specific and
still be a useless vaccine target if its protein yields no peptide that binds
the patient's HLA alleles.

So this step asks a different question: IF this protein is made, how much of
it is visible to the immune system, and in how many people?

That last clause is what makes this an off-the-shelf-vaccine question rather
than a personalised one. A personalised vaccine only needs peptides matching
one patient's six HLA alleles. An off-the-shelf vaccine needs peptides that
work across a population, so we score every candidate against a fixed panel of
common alleles and ask how many of them it can reach.

PIPELINE POSITION
-----------------
  step 3 -> candidates_annotated.csv  (genes that survived DE + tissue safety)
  step 4 -> epitope_results.csv       (this file: how presentable is each one)
  step 6 -> composite score, weight 0.20 on the epitope component

Step 6 combines two columns from here into that 0.20:
    epitope_density      = n_strong_binders / protein_length_aa * 100
    allele_coverage_frac = n_alleles_with_binder / 27
and multiplies the min-max-scaled versions together, so a gene must be BOTH
epitope-dense AND broadly presentable to score well. Neither alone is enough.

Be aware that neither input is clean. n_strong_binders is almost perfectly
collinear with protein length (r = 0.99), which is why step 6 divides it by
length; allele_coverage_frac is still strongly length-driven even though it
looks normalised (spearman 0.89). See the notes in step 4.3.

WHAT THIS STEP DOES NOT DO
--------------------------
Read the numbers as "how visible could this protein be", not "will this work":

  * Binding is not immunogenicity. MHCflurry predicts that a peptide will sit
    in the groove. It says nothing about whether a T-cell exists that can
    recognise it. For self-antigens - which every gene here is - high-avidity
    clones are often deleted in the thymus. This is the single largest gap
    between this score and clinical reality, and nothing in this pipeline
    closes it.
  * No proteome-uniqueness check. A 9-mer from our candidate that also occurs
    verbatim in some other human protein is both a tolerance problem and an
    off-target-toxicity risk. Should be BLASTed against the human proteome
    before any peptide is synthesised. Not done here.
  * Canonical isoform only (see fetch_uniprot_sequence).
  * No cross-check against a second predictor (netMHCpan-4.1).

Input:  $SMM_DATA_DIR/candidates_annotated.csv
Output: $SMM_DATA_DIR/epitope_results.csv     (one row per gene - the summary)
        $SMM_DATA_DIR/epitopes_detail.csv.gz  (one row per peptide x allele)
"""

import time
import urllib.request

import numpy as np
import pandas as pd
from mhcflurry import Class1PresentationPredictor

from config import DATA

# ---------------------------------------------------------------- the panel
# 27 common HLA-I alleles: 15 HLA-A, 12 HLA-B. These are the alleles that
# recur in off-the-shelf vaccine design because between them they cover a
# large fraction of most populations - most individuals carry at least one.
#
# Two deliberate limits, both of which bias the coverage number upward:
#
#   1. No HLA-C. HLA-C reaches the cell surface at roughly a tenth the density
#      of HLA-A/B and is less well predicted, so it is conventionally omitted.
#      But it is genuinely presented, and omitting it means a gene that is only
#      visible through HLA-C scores zero here when it should score something.
#
#   2. Allele frequencies are ancestry-dependent, and this panel leans toward
#      alleles common in European and East Asian populations. "Covers 24 of 27
#      alleles" is therefore NOT the same as "covers 89% of patients" - real
#      population coverage has to be computed against the HLA distribution of
#      the intended trial population, which we have not done.
ALLELE_PANEL = [
    "HLA-A*01:01", "HLA-A*02:01", "HLA-A*02:03", "HLA-A*02:06", "HLA-A*03:01",
    "HLA-A*11:01", "HLA-A*23:01", "HLA-A*24:02", "HLA-A*26:01", "HLA-A*30:01",
    "HLA-A*31:01", "HLA-A*32:01", "HLA-A*33:01", "HLA-A*68:01", "HLA-A*68:02",
    "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01", "HLA-B*18:01", "HLA-B*27:05",
    "HLA-B*35:01", "HLA-B*38:01", "HLA-B*39:01", "HLA-B*40:01", "HLA-B*44:02",
    "HLA-B*51:01", "HLA-B*57:01",
]

# ------------------------------------------------------------- the threshold
# 50 nM is the long-standing IEDB/netMHC convention for a "strong binder".
#
# Known weakness, and it matters for n_alleles_with_binder specifically: a
# fixed nM cut is not equally strict across alleles. Promiscuous alleles like
# A*02:01 bind many peptides tightly and clear 50 nM easily; fussier alleles
# like A*01:01 and B*57:01 rarely do. So this threshold systematically
# over-counts binders for some alleles and under-counts for others, which
# means allele_coverage_frac is partly measuring "how many promiscuous alleles
# are in the panel" rather than purely a property of the gene.
#
# MHCflurry's own documentation prefers a percentile rank cut (e.g. <0.5% or
# <2%), which is allele-normalised by construction and would remove this bias.
# Switching would re-rank the epitope component, so it is left as a documented
# limitation rather than changed silently.
STRONG_NM = 50.0


def fetch_uniprot_sequence(accession, retries=4):
    """Fetch a protein sequence from UniProt. Returns None after all retries fail.

    Retries with exponential backoff (1s, 2s, 4s) because this loop is the
    pipeline's one hard dependency on a flaky external service, and a silent
    failure here is worse than it looks: a gene whose fetch times out is
    dropped from `sequences`, never scored, and never reaches the shortlist -
    with nothing in the final output recording that it was ever a candidate.
    That makes the result depend on network luck rather than on the data.
    An earlier run lost MFAP3L and RPL36A exactly this way.

    The caller prints a loud banner listing anything lost, so a degraded run
    is visible rather than merely quieter.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                lines = resp.read().decode().splitlines()
            # FASTA: drop the '>' header line, concatenate the sequence lines.
            # This returns the CANONICAL (reviewed Swiss-Prot) isoform only.
            # Alternative splice products are not fetched, so a junction-
            # spanning neoepitope from a tumour-specific isoform is invisible
            # to this step.
            return "".join(line for line in lines if not line.startswith(">"))
        except Exception as e:
            if attempt == retries - 1:
                print(f"  WARNING: could not fetch {accession} "
                      f"after {retries} attempts: {e}")
                return None
            time.sleep(2 ** attempt)


def make_peptides(seq, lengths=(9, 10)):
    """Every contiguous 9-mer and 10-mer in the protein (a sliding window).

    We generate ALL of them rather than pre-filtering, because which peptides
    are worth looking at is exactly what the predictor decides.

    Two things follow from this that matter downstream:

      * Peptide count is linear in protein length: (L-8) + (L-9) = 2L - 17.
        A 1000 aa protein yields ~1983 peptides, a 200 aa protein ~383. So
        n_strong_binders is roughly proportional to length and CANNOT be
        compared across genes raw - a big protein wins automatically. Step 6
        divides it by protein_length_aa for exactly this reason.

      * Only 9- and 10-mers. HLA-I also presents 8-mers and 11-mers, but
        9-mers dominate real immunopeptidomes (~60-70%) with 10-mers next, so
        this captures most of the repertoire at half the compute.

    Note also that a sliding window ignores whether the proteasome would ever
    cut there. MHCflurry's presentation model scores processing separately
    (processing_score), which partly compensates.
    """
    peps = []
    for length in lengths:
        for i in range(len(seq) - length + 1):
            peps.append(seq[i:i + length])
    return peps


# ================================================================= STEP 4.1
# Load the surviving candidates and fetch one protein sequence per gene.
# ==========================================================================
annot = pd.read_csv(f"{DATA}/candidates_annotated.csv")

# A gene with no UniProt accession cannot be scored here at all. This drops
# lncRNAs, pseudogenes and anything HPA could not map to a protein - which is
# correct for a peptide vaccine (no protein, no epitope), but note it is a
# silent exclusion criterion that never appears in the ranked output.
annot = annot[annot["hpa_uniprot_id"].notna()]
print(f"{len(annot)} candidates with a UniProt protein")

sequences = {}
for gene, acc in zip(annot["gene"], annot["hpa_uniprot_id"]):
    # ".split('.')[0]" strips a version suffix (P12345.2 -> P12345); the REST
    # API wants the bare accession.
    seq = fetch_uniprot_sequence(str(acc).split(".")[0])
    if seq:
        sequences[gene] = seq
        print(f"  {gene} ({acc}): {len(seq)} aa")

# Make any loss impossible to miss. Without this, a gene lost to a network
# error looks identical to a gene that was never a candidate.
lost = sorted(set(annot["gene"]) - set(sequences))
if lost:
    print("\n" + "!" * 70)
    print(f"!! {len(lost)} CANDIDATE(S) LOST TO FAILED UNIPROT FETCHES:")
    print(f"!! {', '.join(lost)}")
    print("!! These genes will be absent from epitope_results.csv and will")
    print("!! therefore never appear in the final ranking. Re-run before")
    print("!! treating this shortlist as complete.")
    print("!" * 70 + "\n")


# ================================================================= STEP 4.2
# Predict presentation: every peptide against every allele in the panel.
# ==========================================================================
print("\nRunning MHCflurry (this takes several minutes)...")
predictor = Class1PresentationPredictor.load()

t0 = time.time()
detail_rows = []
for gene, seq in sequences.items():
    peps = make_peptides(seq)

    # THE PER-ALLELE LOOP - the subtlety in this whole script.
    #
    # Class1PresentationPredictor.predict(alleles=[...]) interprets the list
    # as ONE PERSON'S GENOTYPE (up to 6 alleles), and returns a single row per
    # peptide reporting only the best-scoring allele of that genotype. So
    # passing all 27 at once would collapse to one row per peptide and throw
    # away precisely what we need: WHICH alleles can present this protein.
    #
    # Looping one allele at a time instead gives 27 rows per peptide, one per
    # allele, so allele coverage becomes countable. It also means that within
    # each call `best_allele` is trivially the single allele passed in - which
    # is what makes the nunique() count in step 4.3 correct.
    #
    # Cost: 27 predictor calls per gene instead of 1. This is the slow part.
    per_allele = [predictor.predict(peptides=peps, alleles=[a]) for a in ALLELE_PANEL]

    df = pd.concat(per_allele, ignore_index=True)
    df["gene"] = gene
    detail_rows.append(df)
    print(f"  {gene}: {len(peps)} peptides x {len(ALLELE_PANEL)} alleles "
          f"({time.time() - t0:.0f}s elapsed)")

# MHCflurry returns, per peptide x allele:
#   affinity           - predicted binding affinity in nM, LOWER IS TIGHTER
#   processing_score   - likelihood the peptide is produced and transported
#   presentation_score - 0-1, combines affinity + processing; the model's own
#                        preferred single number
# We threshold on raw `affinity` (the conventional 50 nM cut) and carry the
# presentation score along as a summary only.
detail = pd.concat(detail_rows, ignore_index=True)

# Gzipped because this is large: ~38 genes x ~1000 peptides x 27 alleles is
# on the order of a million rows. Kept as the audit trail - it is what lets
# you pull the actual peptide sequences for any gene later.
detail.to_csv(f"{DATA}/epitopes_detail.csv.gz", index=False, compression="gzip")


# ================================================================= STEP 4.3
# Collapse a million peptide-allele rows into one row per gene.
# ==========================================================================
print("\nSummarizing per gene...")
summary = []
for gene, g in detail.groupby("gene"):
    strong = g[g["affinity"] < STRONG_NM]

    # How many DISTINCT alleles in the panel have at least one strong binder.
    # Works because of the per-allele loop above: in a single-allele call,
    # best_allele is always that allele, so the distinct values here are
    # exactly the alleles that produced a sub-50 nM peptide.
    alleles_covered = strong["best_allele"].nunique()

    summary.append({
        "gene": gene,
        # Carried forward so step 6 can length-normalise n_strong_binders.
        "protein_length_aa": len(sequences[gene]),
        "n_peptides_screened": g["peptide"].nunique(),

        # RAW COUNT - length-confounded, do not compare across genes directly.
        # Step 6 converts this to epitope_density (binders per 100 aa).
        "n_strong_binders": len(strong),

        # THE POPULATION-COVERAGE PROXY, 0-27 and its fraction. Step 6 uses
        # allele_coverage_frac RAW, and on this candidate set that turns out
        # to be mostly a protein-length measurement in disguise:
        #
        #     spearman(protein_length_aa, allele_coverage_frac) = 0.89
        #
        # A longer protein gets more independent chances to hit each allele.
        # Concretely: DST (7570 aa) reaches 26/27 alleles, RPL36A (106 aa)
        # reaches 6/27 - a ranking that says more about size than about
        # antigenicity. Observed range is 6 to 26 of 27, median 20, and
        # nothing hits 27/27, so the metric is NOT saturated - it discriminates
        # fine, it just substantially discriminates on the wrong thing.
        #
        # Partial accidental correction: step 6 forms the epitope component as
        # scale(allele_coverage_frac) * scale(epitope_density), and
        # epitope_density (binders per 100 aa) favours SHORT proteins while
        # coverage favours long ones, so the product cancels some of the bias.
        # That is luck, not design. The principled fix is to score coverage
        # against a length-matched null - listed as an open item, not done.
        "n_alleles_with_binder": alleles_covered,
        "allele_coverage_frac": alleles_covered / len(ALLELE_PANEL),

        # The single tightest peptide and the allele that binds it. Useful as
        # a lead peptide to take to the bench, but noisy as a score: it is one
        # extreme value out of ~1000 predictions, so it rewards outliers.
        # Not used in the composite for that reason.
        "best_affinity_nM": g["affinity"].min(),
        "best_peptide": g.loc[g["affinity"].idxmin(), "peptide"],
        "best_allele": g.loc[g["affinity"].idxmin(), "best_allele"],

        # Recorded for the evidence trail; not used in the composite.
        "max_presentation_score": g["presentation_score"].max(),
    })

# Sorted by raw binder count for readability only - this is NOT the ranking.
# The actual epitope component is computed in step 6, length-normalised.
epi = pd.DataFrame(summary).sort_values("n_strong_binders", ascending=False)
epi.to_csv(f"{DATA}/epitope_results.csv", index=False)
print(epi.to_string(index=False))
