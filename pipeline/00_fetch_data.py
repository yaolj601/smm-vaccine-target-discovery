"""
Step 0: Fetch the public input data.

Downloads both public inputs unattended:
  - GSE193531 supplementary files (UMI count matrix + cell-level metadata)
  - Human Protein Atlas proteinatlas.tsv

Already-present files are skipped, so this is safe to re-run.

Output: $SMM_DATA_DIR/counts.csv.gz, $SMM_DATA_DIR/meta.csv.gz
        $SMM_REF_DIR/proteinatlas.tsv
"""

import shutil
import sys
import urllib.request
import zipfile

from config import COUNTS_CSV, DATA, META_CSV, REF

GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193531/suppl"
DOWNLOADS = [
    (f"{GEO_BASE}/GSE193531_umi-count-matrix.csv.gz", COUNTS_CSV),
    (f"{GEO_BASE}/GSE193531_cell-level-metadata.csv.gz", META_CSV),
]
HPA_ZIP_URL = "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"
HPA_TSV = REF / "proteinatlas.tsv"
HPA_CONSENSUS_URL = "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip"
HPA_CONSENSUS_TSV = REF / "rna_tissue_consensus.tsv"


def download(url, dest):
    """Stream `url` to `dest`, skipping if it already exists."""
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  skip (exists): {dest.name}")
        return
    print(f"  downloading {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        with urllib.request.urlopen(url, timeout=120) as resp, open(tmp, "wb") as fh:
            shutil.copyfileobj(resp, fh)
    except Exception as e:
        tmp.unlink(missing_ok=True)
        raise SystemExit(f"FAILED to download {url}\n  {e}")
    tmp.replace(dest)
    print(f"  -> {dest} ({dest.stat().st_size / 1e6:.1f} MB)")


print(f"GEO supplementary files -> {DATA}")
for url, dest in DOWNLOADS:
    download(url, dest)

print(f"\nHuman Protein Atlas -> {REF}")
if HPA_TSV.exists() and HPA_TSV.stat().st_size > 0:
    print(f"  skip (exists): {HPA_TSV.name}")
else:
    zip_path = REF / "proteinatlas.tsv.zip"
    download(HPA_ZIP_URL, zip_path)
    print("  extracting...")
    with zipfile.ZipFile(zip_path) as zf:
        member = next((n for n in zf.namelist() if n.endswith(".tsv")), None)
        if member is None:
            raise SystemExit(f"no .tsv inside {zip_path}")
        with zf.open(member) as src, open(HPA_TSV, "wb") as dst:
            shutil.copyfileobj(src, dst)
    zip_path.unlink()
    print(f"  -> {HPA_TSV} ({HPA_TSV.stat().st_size / 1e6:.1f} MB)")

# per-tissue nTPM for every tissue - needed for the outlier test in step 3.
# proteinatlas.tsv only lists a gene's ENRICHED tissues, which is not enough to
# ask "is this gene a significant outlier in bone marrow?".
print(f"\nHPA consensus tissue nTPM -> {REF}")
if HPA_CONSENSUS_TSV.exists() and HPA_CONSENSUS_TSV.stat().st_size > 0:
    print(f"  skip (exists): {HPA_CONSENSUS_TSV.name}")
else:
    zp = REF / "rna_tissue_consensus.tsv.zip"
    download(HPA_CONSENSUS_URL, zp)
    print("  extracting...")
    with zipfile.ZipFile(zp) as zf:
        member = next((x for x in zf.namelist() if x.endswith(".tsv")), None)
        if member is None:
            raise SystemExit(f"no .tsv inside {zp}")
        with zf.open(member) as src, open(HPA_CONSENSUS_TSV, "wb") as dst:
            shutil.copyfileobj(src, dst)
    zp.unlink()
    print(f"  -> {HPA_CONSENSUS_TSV} ({HPA_CONSENSUS_TSV.stat().st_size / 1e6:.1f} MB)")

print("\nAll inputs present.")
