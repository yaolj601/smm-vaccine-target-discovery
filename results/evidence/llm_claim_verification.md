# Verification of LLM-asserted external claims

The novelty/confidence grades in `llm_scores.json` are model knowledge, not
retrieved citations. This file records which of those assertions have actually
been checked against a primary source, and what the check returned.

`llm_scores.json` is **not** edited to match — it is the evidence trail of what
the model said, and rewriting it would destroy the record. Corrections live
here and in the derived prose (WRITEUP.md, the slides).

---

## DKK1 — NCT03591614

Checked 2026-09-13 against the ClinicalTrials.gov API v2.

**Model asserted** (`existing_programs`):
> pUMVC3-hDKK1 DNA vaccine in smouldering myeloma (Mayo Clinic, NCT03591614);
> also anti-DKK1 mAb BHQ880 (Novartis) in MM phase 1b/2 trials, and DKN-01 in
> other solid tumours

**Registry returns:**

| Field | Value |
|---|---|
| Title | Dendritic Cell DKK1 Vaccine for Monoclonal Gammopathy and Stable or Smoldering Myeloma |
| Sponsor | Case Comprehensive Cancer Center |
| Phase | Early Phase 1 |
| Status | **WITHDRAWN** — "PI left institution" |
| Enrollment | **0** |
| Start date | 2023-12-01 |
| Intervention | DKK1 peptide-loaded autologous dendritic cell vaccine, 5–10×10⁶ in 0.5 mL Plasma-Lyte A + 5% HSA |

**Verdict: the grade itself is wrong, not just the details.**

A full DKK1 query on ClinicalTrials.gov returns **exactly one vaccine ever
registered** — NCT03591614 above — and it never dosed a patient. DKK1's real
clinical footprint is antibodies:

| Trial | Agent | Indication | Phase | Status | n |
|---|---|---|---|---|---|
| NCT01302886 | BHQ880 | **smouldering myeloma** | 2 | COMPLETED | 41 |
| NCT01337752 | BHQ880 + bortezomib/dex | MM, renal insufficiency | 2 | COMPLETED | 9 |
| NCT01711671 | DKN-01 + len/dex | MM | 1 | COMPLETED | 8 |

So `1b_vaccine_in_mm` is unsupported — no vaccine has been *in clinical study*,
only registered and abandoned. The correct tier is `2_clinical_other`: an
antibody programme exists in myeloma, and the vaccine space is open. That is
precisely the "validated but not yet vaccinated against" case the grid rewards.

**Effect:** novelty_confidence 3 → 8; DKK1 rank 28 → 17. Top 15 unchanged.
Applied via `data/verified_overrides.json`; `llm_scores.json` untouched.

What the model got right and wrong, separated:
- **Right** — DKK1 has prior clinical development in myeloma, and there is a
  real NCT number in our exact indication. No expression data could supply that.
- **Wrong** — the modality (dendritic-cell, not the `pUMVC3-hDKK1` DNA plasmid
  named), the sponsor (Case Comprehensive, not Mayo), the trial's status
  (withdrawn, zero enrolled), and consequently the tier.

**The generalisable lesson:** a fabricated identifier was never the failure
mode. The identifier was real and the indication was right; what drifted was
the detail around it — and the detail was what set the score. Treat an NCT
number as a lookup key, never as the claim, and check the *status* field, not
just that the trial exists.

---

## Not yet checked

- **NLGN4X** — the NLGN4Y H-Y minor-histocompatibility rationale (rank 3, and
  the single most load-bearing unverified claim left, since it is the entire
  reason NLGN4X scores as it does).
- **KIT** — asserted as an established MM diagnostic/flow marker (rank 6).
- **CCND1, MYEOV, IGF1, HGF** — asserted MM literature; high prior, unchecked.
- Every other `existing_programs` string in `llm_scores.json`.
