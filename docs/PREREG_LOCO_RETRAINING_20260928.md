# Preregistration: leave-cohort-out retraining across the harmonized external cohorts (2026-09-28)

Queue item 11 (docs/QUEUE_20260927.md): "Leave-cohort-out retraining - PARTIAL:
leave-one-cohort-out evaluation exists (masked coefficients); full retraining per
held-out cohort NOT DONE." Verdict weakness (docs/USER_VERDICT_20260927.md): "No
leave-cohort-out retraining." Written 2026-09-28 before any leave-cohort-out model
was fit. Thresholds and the decision rule below are locked; they will not move
after outcomes.

## Question
Is the 70-gene panel's external discrimination stable to FULL retraining, or is
the single frozen METABRIC fit itself load-bearing? The existing
leave-one-cohort-out artifact (results/transport_stability.json) only masks the
frozen METABRIC coefficients per cohort; it never refits. This branch refits.

## Scope (fixed)
- Strata and cohorts exactly as the committed harmonization record
  (results/endpoint_harmonization.json): HARM-RFS = {GSE2990, GSE7390, GSE20685};
  HARM-DMFS = {GSE2990, GSE7390, GSE11121, GSE20685, GSE25066}. No other cohorts;
  METABRIC excluded (discovery cohort with its own benchmark; this branch concerns
  the external transport pool only).
- Endpoints: the committed key-verbatim derivations in
  run_endpoint_harmonization.py, unchanged. Per-cohort n, events and genes_mapped
  must reproduce the committed record exactly or the run aborts before any
  leave-cohort-out fit (integrity gate, below).
- Features: the 70-gene expression panel only (expr-only, matching the transported
  model; no clinical covariates). Per-cohort probe-to-gene mapping and per-gene
  z-scoring within each cohort's analysis set exactly as the committed extraction.
  Genes unmapped in a cohort are all-zero columns there and contribute zero to any
  risk score, matching the committed frozen-transport convention.

## Model and procedure (fixed)
- Frozen comparator: the committed METABRIC expr-only penalized Cox (penalizer
  0.05), refit with the identical fitting code as run_endpoint_harmonization.py.
  Integrity gate: recomputed per-cohort frozen C-indices must match the committed
  endpoint_harmonization.json values to 4 decimals; any mismatch aborts the run
  before any leave-cohort-out fit and is reported as a blocker, not a verdict.
- Leave-cohort-out arm: for each cohort C in a stratum, fit lifelines
  CoxPHFitter (penalizer 0.05, no strata term) on the pooled patients of the
  stratum's OTHER cohorts (each cohort z-scored within itself first, as
  committed). Zero-variance genes in the training pool are dropped from that
  arm's fit; their coefficients are effectively zero on the held-out cohort,
  mirroring the frozen convention. Evaluate the fitted model on C with the
  committed lane statistic: concordance_index(t, -eta, e).
- Primary readout per cohort: delta_C = LOCO C-index minus the committed frozen
  per-cohort C-index for the same cohort and stratum (values in
  results/endpoint_harmonization.json).
- Uncertainty: paired patient bootstrap within the held-out cohort - resample
  patients of C, recompute both the LOCO and the frozen C-index on the resample,
  delta per replicate; B = 1000, percentile 95% CI, one shared rng stream seeded
  0, cohorts processed in the fixed order GSE2990, GSE7390, GSE20685 (HARM-RFS)
  then GSE2990, GSE7390, GSE11121, GSE20685, GSE25066 (HARM-DMFS).

## Locked decision rule
- A stratum is RETRAINING-STABLE if every one of its cohorts has a paired
  bootstrap 95% CI of delta_C that includes 0.
- If any cohort's CI excludes 0, the stratum is NOT RETRAINING-STABLE; the
  direction (retraining helps or harms, per cohort) is reported plainly.
- The two strata are verdicted separately and never pooled into one verdict.
  A NOT-STABLE result is reported as measured; it is not reframed, and the
  frozen-transport numbers already in the paper are not altered.
- Abort conditions (no verdict, reported as execution failure): integrity-gate
  mismatch; any Cox fit failing to converge; any cohort/stratum combination
  missing from the committed record.

## Secondary (descriptive only, no gate)
- The LOCO C-indices themselves, with paired-bootstrap CIs, per cohort.
- Per-stratum n-weighted mean of per-cohort LOCO C-indices with a
  cohort-stratified patient bootstrap (500 reps, same rng stream), compared
  descriptively against the committed pooled frozen values (0.6165 HARM-RFS,
  0.6477 HARM-DMFS).

## Declared asymmetry
GSE20685 maps 69/70 panel genes (GPL570); the GPL96 cohorts map 63/70. When
GSE20685 is held out, its 6 GPL570-only genes get no training-pool coefficient
and contribute zero to its LOCO score, while the committed frozen comparator
scored them. This mirrors the unmapped-gene convention and is disclosed, not
adjusted away; no extra arms are run.

## Deliverables (either verdict)
- run_loco_retraining.py (committed before any outcome), results/loco_retraining.json.
- Queue item-11 addendum and a paper subsection addendum reporting the verdict
  as measured.
