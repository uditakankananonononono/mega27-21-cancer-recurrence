# Pre-registration: fresh prospective transport test on an untouched cohort (GSE17705 primary)

Date: 2026-09-28 IST (committed before the designated cohort's expression data is
fetched or scored; only GEO series-matrix headers - field names, sample counts,
series summaries - were inspected to establish eligibility, no expression values
and no outcome quantities).
Lane 21 verdict item 5 follow-up. The 2026-09-27 prospective-holdout execution
(GSE21653) was ruled invalid for the locked transport claim because GSE21653's
native endpoint is DFS (parent ruling 2026-09-28 04:51 IST); the locked claim
(recurrence-endpoint transport of the frozen METABRIC model) remains UNTESTED on
an unused cohort. The ruling requires a NEW dated prereg - this document - before
any outcome, with either a DFS-aware claim or a cohort with native RFS/DRFS/DMFS.
This prereg takes the second path: DFS-only cohorts are ineligible here.

## Locked model and genes (unchanged)
- Score: the frozen penalized-Cox risk model committed as results/risk_model.json
  (87 features: 17 clinical + 70-gene panel; trained on 1,580 METABRIC patients).
  SHA of the file at lock time: recorded in the execution commit.
- No refitting, no recalibration, no feature changes on the holdout.

## Holdout cohort eligibility (locked)
A cohort qualifies only if ALL hold:
1. Not previously analyzed in this repository. Excluded: METABRIC, TCGA-BRCA,
   GSE7390, GSE2990, GSE11121, GSE25066, GSE20685, GSE2034, GSE21653, GSE3494,
   GSE31519.
2. Breast cancer, primary tumors, with a NATIVE recurrence-type endpoint
   (RFS/DRFS/DMFS declared in the cohort's own sample records, not harmonized
   away). DFS-only cohorts are ineligible (per the 2026-09-28 04:51 ruling).
   >= 150 patients, >= 40 events on that endpoint.
3. Expression for >= 60 of the 70 panel genes and the core clinical fields.

## Candidate order and scouting evidence (locked)
1. GSE17705 (PRIMARY). GEO series-matrix header inspected 2026-09-28: 298
   ER-positive primary breast tumors, uniformly tamoxifen-treated 5 years,
   platform GPL96; native endpoint fields "distant relapse (1=dr, 0 censored)"
   with "event time (years)" = DRFS; clinical fields include nodal status
   (0/1). Two profiling labs (MD Anderson, Jules Bordet) - recorded in the
   fetch log; the lab breakdown is reported descriptively with no gate.
2. GSE6532 (BACKUP 1, only if GSE17705 fails eligibility at fetch). Native
   e.rfs/t.rfs AND e.dmfs/t.dmfs fields, grade/node/size/age/er/pgr present.
   Multi-platform (GPL96/97/570): the GPL570 submatrix alone has 87 samples
   (< 150, ineligible alone); eligible only as the GPL96+GPL570 combined fetch
   if the combined n >= 150 and events >= 40, with the platform merge recorded
   in the fetch log.
3. GSE19615 is INELIGIBLE (115 samples < 150), recorded here so it is not
   re-scouted.
The first qualifying cohort fetched is the holdout; the fetch log
(docs/TRANSPORT_FRESH_FETCH_LOG_20260928.md) is committed with the result.

## Locked evaluation
1. Primary metric: concordance of the frozen score on the holdout's native
   endpoint, bootstrap 95% CI (B=1000, patient-level), using the committed
   transport convention (point C-index uses concordance_index(tv, -eta, ev)).
2. Secondary: calibration slope, and IPCW 5-year Brier vs cohort-null
   (run_brier_dca.py unchanged), if the endpoint horizon allows. (The frozen
   model carries no baseline hazard; if Brier cannot be computed honestly that
   is recorded, as in the GSE21653 execution.)
3. Decision rule: the transport is supported if the concordance CI excludes
   0.5 AND exceeds the cohort's clinical-only comparator fit on point estimate.
   The comparator uses the NPI/grade/nodal covariates AVAILABLE in that cohort,
   locked here per candidate: GSE17705 -> nodal status only (grade and age are
   absent from its sample records; ER is constant ER+); GSE6532 -> grade +
   node + size + age. Otherwise reported as a non-beat. alpha is descriptive;
   no post-hoc endpoint switching.

## Out of scope
No pooling with existing cohorts, no threshold tuning, no gene substitution,
no endpoint switching, no model changes. Any deviation requires a new dated
commit before outcomes. Execution script (run_transport17705.py for the
primary candidate) is committed before any outcome-producing run.
