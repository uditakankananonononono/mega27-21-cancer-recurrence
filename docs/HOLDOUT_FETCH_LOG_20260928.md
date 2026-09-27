# Prospective holdout fetch log (lane 21, item 5)

Date: 2026-09-28. Per docs/PREREG_PROSPECTIVE_HOLDOUT_20260927.md: the first
qualifying cohort fetched is the holdout.

## Locked model file at execution
- results/risk_model.json SHA-256: 8ae3992542e60468b8506286b3fd8fdd4b8d5ed880e6381f868c525aea257dfc
  (unchanged since the prereg; verified at fetch time)

## Candidates examined, in order
1. GSE31519 - INELIGIBLE. Series matrix holds only 67 samples (< 150 locked
   minimum; the 579 cases in the publication are assembled across many
   series). No recurrence endpoint fields in the series matrix.
   Matrix: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE31nnn/GSE31519/matrix/GSE31519_series_matrix.txt.gz
2. GSE3494 - INELIGIBLE. Per-platform series matrices (GPL96: 251 samples)
   contain no clinical characteristics at all (only tissue), so no declared
   recurrence endpoint is recorded with the series.
   Matrix: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE3nnn/GSE3494/matrix/GSE3494-GPL96_series_matrix.txt.gz
3. GSE21653 - QUALIFIES. DESIGNATED HOLDOUT.
   - Not previously analyzed here (not METABRIC, TCGA-BRCA, GSE7390,
     GSE2990, GSE11121, GSE25066, GSE20685, GSE2034).
   - Breast cancer primary tumors (invasive early breast adenocarcinomas,
     Institut Paoli-Calmettes; tissue: breast cancer tumor), 266 samples
     (>= 150 locked minimum).
   - Recurrence-type endpoint declared natively in the series matrix:
     "dfs evt" / "dfs time (months)"; 83 events (>= 40 locked minimum),
     169 non-events. NOTE (interpretation, flagged): the endpoint is
     disease-free survival rather than a field literally named RFS/DRFS/
     DMFS; it is the cohort's declared recurrence-type endpoint, not
     harmonized away.
   - Platform GPL570 (HG-U133_Plus_2); the endpoint-harmonization work
     mapped 69/70 panel genes on GPL570, satisfying the >= 60 locked
     minimum (to be re-verified on this file during evaluation).
   - Core clinical fields present: age at diagnosis, sbr grade, pn (nodal),
     pt (tumor), er/pr/erbb2/ki67/p53 IHC, molecular subtype.
   - Matrix: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE21nnn/GSE21653/matrix/GSE21653_series_matrix.txt.gz
   - SHA-256: f13e3bf972bd92a6fd53a1fbe42ad7c84caf4651df3a5e836f70dd3b22b5e76f
   - Local copy (git-ignored): data_cache/geo/GSE21653_series_matrix.txt.gz

## Locked evaluation reminder (from the prereg)
Primary: concordance of the frozen score on the holdout endpoint, bootstrap
95% CI (B=1000, patient-level). Secondary: calibration slope, IPCW 5-year
Brier vs cohort-null if horizon allows. Decision rule: transport supported
if the concordance CI excludes 0.5 AND the point estimate exceeds the
cohort's clinical-only comparator (NPI/grade/nodal covariates available in
that cohort). No refitting, no pooling, no threshold tuning, no gene
substitution.
