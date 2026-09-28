# NEW dated prereg: cohort-anchored slope recalibration pivot (2026-09-28, 07:23 IST)

Pivot 6A under the 2026-09-28 user posture (parent adjudication 07:22 IST, order #1).
The failed claims being pivoted AWAY from: gene-level sign stability of the frozen
score (items 5/6 record: strict-CI gate failed, all-four sign gate failed, fresh-cohort
transport NOT SUPPORTED). Those results stay recorded as measured and are never
re-tuned. This prereg asks a NEW question about deployability, not stability:
can each cohort's own data repair the frozen score's calibration out of sample?

## Locked question
For each of the five committed external cohorts, does a ONE-PARAMETER slope
recalibration of the frozen METABRIC 70-gene expression-only risk score - the
single coefficient alpha fit on that cohort's own training half - produce a
held-out calibration slope inside [0.8, 1.25]?

## Locked data and conventions (integrity gate first)
- Cohorts, endpoints, probe mapping, per-cohort z-scoring, missing-gene zero-fill:
  exactly run_transport_stability.load_cohort / COHORTS, unchanged.
- Frozen score: METABRIC expr-only penalized Cox (penalizer 0.05) coefficients,
  refit with the identical fitting code as the committed calibration audit.
- Slope convention: univariate penalized Cox (penalizer 0.05) of the linear
  predictor against the cohort endpoint, matching results/transport_calibration.json.
- INTEGRITY GATE: the recomputed full-cohort frozen slopes must match the
  committed results/transport_calibration.json values to 4 decimals for all five
  cohorts BEFORE any split is scored. Any mismatch aborts; no verdict.

## Locked procedure per cohort
1. Stratified 50/50 split by event indicator, sklearn train_test_split
   random_state=20260928, one fixed split per cohort, declared limitation: no
   repeated splitting.
2. Training half: fit alpha = slope of the frozen eta (penalizer 0.05). If the
   training half has fewer than 5 events the cohort is recorded as skipped
   (never silently re-split).
3. Held-out half: recalibrated score eta' = alpha * eta. Primary readout =
   held-out slope of eta' (penalizer 0.05, same convention). Context readout =
   held-out slope of the un-recalibrated eta on the same patients.
4. Uncertainty: B=1000 patient bootstrap over the held-out half ONLY (alpha
   fixed from training; this is the deployability frame: recalibrate once,
   evaluate out of sample), percentile 95% CI, seed 20260928, one rng stream,
   cohorts processed in the committed order (GSE7390, GSE11121, GSE2990,
   GSE20685, GSE25066). Replicates with fewer than 2 distinct event values are
   skipped and counted.

## Locked decision rule
- Per-cohort PASS: held-out recalibrated slope POINT estimate inside
  [0.8, 1.25]. The band is locked now, before any outcome; it is a
  calibration-adequacy band, not a significance test. Held-out halves carry
  roughly 23-55 events, so 95% CIs will be wide; they are reported for honesty
  but are NOT the gate (declared: a CI-based gate would be structurally
  unpassable at these event counts).
- Item-level verdict: POSITIVE if at least 4 of 5 cohorts PASS. Otherwise the
  per-cohort table is reported as measured with no reframing. Endpoint
  heterogeneity across cohorts is inherited from the committed record and
  declared, not re-harmonized; GSE7390's committed full-cohort slope is
  shallow (0.464), so its alpha is expected to be far below 1 - the test is
  whether that shrinkage replicates out of sample.
- A PASS does not claim clinical utility, endpoint comparability, or transport
  of discrimination; it says only that local one-parameter recalibration was
  stable out of sample. A FAIL leaves the recorded negatives untouched.

## Deliverables (either outcome)
run_recalibration_pivot.py, results/recalibration_pivot.json, queue item-6
addendum and a paper paragraph reporting the table as measured.
