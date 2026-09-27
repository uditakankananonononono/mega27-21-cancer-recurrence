# Pre-registration: rerun of the leave-cohort-out stability test

Date: 2026-09-27 (committed before any rerun outcome is computed).
Lane: 21 (RecurScan). Verdict item: "Pre-register the stability test before re-running".

## Why a new pre-registration
The original stability test (run_transport_stability.py, committed pre-outcome)
failed its own gate: 3 of 5 folds selected zero genes, so the planned five-fold
mean delta and its permutation comparisons were undefined. Per the lane rule,
any rerun requires this new dated pre-registration before outcomes.

## Locked changes relative to the failed run (all decided now, before outcomes)
1. Selection rule per training fold set: a gene is "stable" if its METABRIC
   penalized-Cox coefficient direction matches its sign in at least 4 of 5
   training cohorts (was: bootstrap CI excluding zero in all four). Rationale:
   the zero-gene failure came from CI width, not sign inconsistency; sign
   agreement is the weaker, pre-declared criterion.
2. Minimum panel: if fewer than 3 genes pass in a fold, that fold is recorded
   as gate-failed (undefined reduced score), exactly as before. No mean is
   computed over passing folds only.
3. Gate: the rerun is a success only if all 5 folds yield >= 3 stable genes AND
   the five-fold mean delta (reduced minus full concordance) has a bootstrap
   95% CI excluding 0 in either direction, reported whichever way it lands.
4. The reduced score masks unstable METABRIC coefficients; no retraining on the
   held-out cohort. Same cohorts, same preprocessing, same alpha (0.05 for this
   descriptive stability gate, fixed here).

## Locked interpretation
- Success means only: a sign-stable subpanel transports without losing
  concordance. It does not license pathway, mechanism, or clinical claims.
- Failure is reported as a second bounded negative; the full-panel comparator
  audit is unaffected either way.
- The result lands in results/transport_stability_rerun.json with this file's
  commit hash recorded inside the JSON.

## Out of scope
No new cohorts, no endpoint redefinition, no penalizer changes, no gene
substitutions after outcomes are seen.
