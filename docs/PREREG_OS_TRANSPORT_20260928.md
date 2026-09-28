# Pre-registration: OS-endpoint transport of the frozen METABRIC panel in the designated cohorts (GSE7390 + GSE20685)

Date: 2026-09-28 IST (committed before any OS-endpoint scoring; the only OS
facts inspected are the eligibility counts in the committed scouting record
d2064b0/d4f2cd8 - field presence and death counts, no outcome quantities).

Lane 21 verdict item 5, pivot A second step per the adjudicated order
2026-09-28 07:22 IST ("5A OS-endpoint scouting then the designated-cohort
question if a cohort clears scouting"). Scouting (POSITIVE) designated two
cohorts. This is a NEW question - endpoint transport - not a re-tuning of the
failed DMFS fresh-transport claim (GSE17705, measured NOT SUPPORTED).

## The question (new)
Does the frozen METABRIC risk score carry prognostic concordance for
NATIVE OVERALL SURVIVAL in platform-matched cohorts whose OS endpoint has
never been evaluated in this repository?

## Honest scope statement (locked)
- The frozen score (results/risk_model.json; 87 features: 17 clinical +
  70-gene panel; trained on 1,580 METABRIC patients) was never fitted,
  tuned, or selected on GSE7390 or GSE20685. Evaluating it there is unbiased
  with respect to model selection.
- Both cohorts were previously used in this repository for RFS/DMFS-endpoint
  evaluation of the same frozen score; their OS endpoints are UNTOUCHED
  (every prior locked run used RFS/DRFS/DMFS endpoints only). This prereg
  therefore claims ENDPOINT transport, not fresh-cohort transport.
- The 6A recalibration fitted intercepts on these cohorts' DMFS endpoint.
  Recalibration parameters are NOT part of this evaluation: the frozen score
  only, no recalibration, no refitting, no feature changes.
- GSE7390 time-unit note from scouting: stored t.os median 4561.5 suggests
  days. Units are verified verbatim at fetch and recorded in the fetch log;
  concordance is unit-invariant, so the primary metric is unaffected.

## Designated cohorts (locked, from the scouting gate)
1. GSE7390 (GPL96): 198 patients, native t.os/e.os, 56 deaths.
2. GSE20685 (GPL570): 327 patients, native follow_up_duration/event_death,
   83 deaths.
No other cohort is added after this designation.

## Locked evaluation
1. Primary metric per cohort: concordance of the frozen score on native OS,
   bootstrap 95% CI (B=1000, patient-level, seed 0), committed transport
   convention (point C-index uses concordance_index(tv, -eta, ev)).
2. Comparator per cohort: clinical-only Cox fit in-cohort on the same OS
   endpoint, using the covariates available in that cohort's committed
   clinical extract (GSE7390: age, size, node, grade, er; GSE20685: age,
   t_stage, n_stage). In-cohort fitting favors the comparator, so the
   comparison is conservative for the panel. Feature imputation, if any,
   declared in the fetch log before scoring.
3. Brier: not computed (frozen model carries no baseline hazard; recorded as
   before).
4. Decision rule (locked): OS-endpoint transport is SUPPORTED only if the
   concordance CI excludes 0.5 AND the point estimate exceeds the
   clinical-only comparator point estimate in BOTH designated cohorts. A
   mixed result is reported as measured (per-cohort), never pooled into a
   PASS.
5. Missingness: panel genes map per platform as in the committed transport
   machinery; cohorts with < 60 of 70 panel genes mapped would have failed
   eligibility (both cohorts already meet the expression gate from prior
   locked runs: gene mapping identical to the DMFS/RFS evaluations).

## Negative handling (user directive 2026-09-28 07:18 IST)
A measured null is reported honestly as supporting rigor; if SUPPORTED, the
positive leads the queue line and paper paragraph; if NOT SUPPORTED, the
pivot family records the result and this question closes - it is not
re-tuned.
