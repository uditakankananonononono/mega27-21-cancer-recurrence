# NEW dated preregistration: cohort-stratified joint training (pivot 11A) (2026-09-28, 07:37 IST)

Pivot 11A under the 2026-09-28 posture; parent adjudication order #2
(2026-09-28 07:22 IST). The claim being pivoted AWAY from: item 11's LOCO
retraining-stability claim - NOT retraining-stable in either stratum
(e1b5ba2); the frozen METABRIC fit is load-bearing. That record stands and is
never re-tuned. This prereg asks a NEW method question: instead of leaving a
cohort OUT, does training JOINTLY on cohort-stratified folds - every training
fold contains patients from every cohort in the stratum - beat the frozen
panel on held-out patients from the same external pool?

## Locked question
Within a harmonized stratum, is the pooled out-of-fold C-index of a
cohort-stratified JOINTLY trained penalized Cox (expr-only 70-gene panel)
higher than the frozen METABRIC panel's C-index on the same held-out
patients?

## Scope and data (fixed)
- Strata and cohorts exactly as the committed harmonization record:
  HARM-DMFS = {GSE2990, GSE7390, GSE11121, GSE20685, GSE25066} (PRIMARY);
  HARM-RFS = {GSE2990, GSE7390, GSE20685} (SECONDARY). Never pooled together.
- Endpoints, probe-to-gene mapping, per-cohort z-scoring, zero-fill of
  unmapped genes, frozen comparator: exactly the committed LOCO loading path.
  The script reuses run_loco_retraining.py's committed loading code by
  source extraction up to its analysis marker (no edit to the locked
  artifact), so the data build is byte-identical to e1b5ba2's.
- INTEGRITY GATE: the committed per-cohort n / events / genes_mapped and
  frozen per-cohort C-indices must reproduce to 4 decimals before any joint
  fit; any mismatch aborts and is reported as a blocker, not a verdict.

## Locked procedure
1. Per stratum, pool the cohorts' patients. One fixed 5-fold partition:
   sklearn StratifiedKFold(n_splits=5, shuffle=True, random_state=20260928)
   on cohort-x-event labels, so every fold contains patients from every
   cohort with preserved event ratios. One partition per stratum; declared
   limit: no repeated splitting.
2. Per fold: fit CoxPHFitter(penalizer 0.05) on the pooled training patients
   (all cohorts jointly). Zero-variance genes in the training pool are
   dropped from that fold's fit (committed convention). Predict the linear
   predictor on the held-out fold. Feasibility check: every fold's training
   pool must retain >= 5 events per cohort (guaranteed by the stratification
   given committed event counts); any violation aborts - never re-split
   silently.
3. After all folds, every patient has exactly one out-of-fold joint score.
   Readout: pooled OOF C-index, joint minus frozen, per stratum, using the
   committed lane statistic concordance_index(t, -eta, e).
4. Uncertainty: paired patient bootstrap over the stratum's pooled patients
   on the OOF scores (models are NOT refit inside the bootstrap; this is the
   out-of-fold evaluation frame), B = 2000, percentile 95% CI, one rng stream
   seeded 20260928, strata processed in the fixed order HARM-DMFS then
   HARM-RFS. Replicates with a degenerate resample (undefined C-index) are
   skipped and counted.
5. Descriptive only (no gate): per-cohort OOF deltas with B=500 bootstrap
   from the same stream; per-fold coefficient sign-agreement counts.

## Locked decision rule
- HARM-DMFS (primary): POSITIVE iff the pooled OOF delta 95% CI excludes
  zero in the positive direction. Otherwise the result is a measured
  negative for the joint-training question, reported plainly.
- HARM-RFS (secondary): same rule, verdicted separately; the two strata are
  never combined into one verdict.
- A POSITIVE is explicitly a pooled-joint-method result within the external
  pool: joint training saw every cohort during training, so it says NOTHING
  about transport to an unseen cohort, and it does not repair item 11's LOCO
  negative (different question). Any new-cohort transport claim for the
  jointly trained model would require its own dated prereg.
- A negative does not weaken any committed positive (frozen benchmark
  comparisons stand on their own preregs).

## Deliverables (either outcome)
run_joint_training.py, results/joint_training.json, queue item-11 addendum
and a paper paragraph reporting the strata verdicts as measured.
