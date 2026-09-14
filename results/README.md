Outputs from the run described in ../WRITEUP.md. `run_all.sh` overwrites this
directory in place, so these files are the reference result to diff against.

`evidence/` holds the complete intermediate trail behind the shortlist: the DE
results for all four contrasts, the per-cell Wilcoxon sensitivity check, the
HPA annotations, the per-gene epitope summary, and the LLM scores. Every number
in `ranked_targets.csv` can be traced back through these without re-running.
