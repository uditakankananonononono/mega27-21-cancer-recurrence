# Pre-registration: fully locked prospective holdout with pre-registered genes

Date: 2026-09-27 (committed before the holdout cohort is fetched or scored).
Lane 21 verdict item: "Add a fully locked prospective holdout with pre-registered genes."

## Locked model and genes (pre-registered)
- Score: the frozen penalized-Cox risk model committed as results/risk_model.json
  (87 features: 17 clinical + 70-gene panel; trained on 1,580 METABRIC patients).
  SHA of the file at lock time: recorded in the execution commit.
- No refitting, no recalibration, no feature changes on the holdout.

## Holdout cohort eligibility (locked)
A cohort qualifies only if ALL hold:
1. Not previously analyzed in this repository (excludes METABRIC, TCGA-BRCA,
   GSE7390, GSE2990, GSE11121, GSE25066, GSE20685, GSE2034 - all already used).
2. Breast cancer, primary tumors, with a recurrence-type endpoint (RFS/DRFS/DMFS
   declared and recorded, not harmonized away), >= 150 patients, >= 40 events.
3. Expression for >= 60 of the 70 panel genes and the core clinical fields.
Candidate sources: SCAN-B (GEO), METABRIC-like public series on GEO/ArrayExpress
matching the above. The first qualifying cohort fetched is the holdout; the
fetch log is committed with the result.

## Locked evaluation
1. Primary metric: concordance of the frozen score on the holdout endpoint,
   bootstrap 95% CI (B=1000, patient-level).
2. Secondary: calibration slope, and IPCW 5-year Brier vs cohort-null
   (run_brier_dca.py unchanged), if the endpoint horizon allows.
3. Decision rule: the transport is supported if the concordance CI excludes 0.5
   AND exceeds the cohort's clinical-only comparator fit (NPI/grade/nodal
   covariates available in that cohort) on point estimate; otherwise reported
   as a non-beat. alpha is descriptive; no post-hoc endpoint switching.

## Out of scope
No pooling with existing cohorts, no threshold tuning, no gene substitution.
Any deviation requires a new dated commit before outcomes.
