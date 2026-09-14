"""
Shared configuration for the SMM vaccine-target pipeline.

Every path is resolved relative to the repository root so the pipeline runs
anywhere, and every path can be overridden with an environment variable if
your data lives elsewhere:

    SMM_DATA_DIR    intermediate + input data      (default <repo>/data)
    SMM_OUT_DIR     final deliverables             (default <repo>/results)
    SMM_REF_DIR     reference datasets (HPA)        (default <repo>/reference)
    SMM_HPA_TSV     explicit path to the HPA annotation table
    SMM_HPA_CONSENSUS  explicit path to HPA consensus tissue nTPM
    SMM_LLM_MODEL   Claude model id used in step 5

Step 5 also reads ANTHROPIC_API_KEY, and ANTHROPIC_WORKSPACE_ID when the key
is not itself scoped to a workspace.

This module is also the single home for the gene-exclusion filter, which used
to be copy-pasted (and subtly wrong) in steps 1 and 3.
"""

import os
import re
from pathlib import Path

_HERE = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get("SMM_PROJECT_ROOT", _HERE.parent))

DATA = Path(os.environ.get("SMM_DATA_DIR", PROJECT_ROOT / "data"))
OUT = Path(os.environ.get("SMM_OUT_DIR", PROJECT_ROOT / "results"))
REF = Path(os.environ.get("SMM_REF_DIR", PROJECT_ROOT / "reference"))

for _d in (DATA, OUT, REF):
    _d.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------ raw input
COUNTS_CSV = DATA / "counts.csv.gz"   # GSE193531_umi-count-matrix.csv.gz
META_CSV = DATA / "meta.csv.gz"       # GSE193531_cell-level-metadata.csv.gz

# ------------------------------------------------- reference datasets (step 3)
HPA_CONSENSUS = Path(os.environ.get(
    "SMM_HPA_CONSENSUS", REF / "rna_tissue_consensus.tsv"))


def hpa_tsv():
    """Locate the Human Protein Atlas table.

    Accepts either the full proteinatlas.tsv download or the column subset
    used during development; both carry the same column names.
    """
    explicit = os.environ.get("SMM_HPA_TSV")
    if explicit:
        return Path(explicit)
    for name in ("proteinatlas.tsv", "proteinatlas_subset.tsv"):
        p = REF / name
        if p.exists():
            return p
    return REF / "proteinatlas.tsv"  # not found: let the caller raise


# --------------------------------------------------------------- step 5 (LLM)
# Default to the current flagship. Override with SMM_LLM_MODEL to use a
# cheaper model (e.g. claude-sonnet-5) - 43 short calls either way.
LLM_MODEL = os.environ.get("SMM_LLM_MODEL", "claude-opus-5")


# ------------------------------------------------------------ gene exclusions
# Genes we never want as vaccine targets: immunoglobulin loci (the tumour's own
# clonal product, and the dominant signal in any plasma-cell DE test), HLA
# genes, and mitochondrial genes.
#
# NOTE: an earlier version of this filter used gene.startswith("IG"), which
# silently discarded IGF1R, IGF2, IGF2BP3, IGSF8, IGFBP7 and friends - none of
# them immunoglobulins, and IGF1R/IGF2BP3 are genuine myeloma target candidates.
# The patterns below match the actual immunoglobulin loci and nothing else.
_EXCLUDE_PATTERNS = (
    r"^IG[HKL][VDJ]\d",    # V/D/J segments: IGHV3-23, IGKV1D-39, IGLJ3
    r"^IG[HKL][VDJ]$",     # bare segment-locus names
    r"^IGH[GAMDE]\d*$",    # heavy constant: IGHG1-4, IGHA1-2, IGHM, IGHD, IGHE
    r"^IGHGP$",            # heavy constant pseudogene
    r"^IGKC$",             # kappa constant
    r"^IGLC\d+$",          # lambda constant: IGLC1-7
    r"^IGLL\d+$",          # surrogate light chain: IGLL1, IGLL5
    r"^HLA-",              # HLA class I/II
    r"^MT-",               # mitochondrially encoded
    r"^MTRNR2L\d+$",       # nuclear mitochondrial rRNA-like pseudogenes
)
_EXCLUDE_RE = re.compile("|".join(_EXCLUDE_PATTERNS))


def is_excluded(gene):
    """True if `gene` is an immunoglobulin, HLA or mitochondrial gene."""
    return bool(_EXCLUDE_RE.match(str(gene).upper()))
