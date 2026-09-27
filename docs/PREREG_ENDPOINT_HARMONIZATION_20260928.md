# Preregistration: patient-level endpoint harmonization, lane 21 (2026-09-28)

Queue item 3 (docs/QUEUE_20260927.md): "Harmonize endpoints before pooling -
PARTIAL ... patient-level harmonization NOT DONE." Written before any
harmonized endpoint is derived or evaluated. Rules below are locked and will
not move after outcomes.

## Question
What is the transport discrimination of the committed 70-gene expr-only Cox
model when every cohort is evaluated only on patient-level endpoint
definitions that are derivable from its raw data, and no pooling ever mixes
endpoint types?

## Locked definitions
- Stratum HARM-RFS: event = any recurrence (locoregional or distant),
  time = the accompanying time field, converted to months where needed.
  A cohort's handling of death without recurrence inside its own RFS field is
  not altered; each cohort's verbatim source definition is recorded in the
  derivation audit table.
- Stratum HARM-DMFS: event = distant metastasis / distant relapse only,
  time = time to that event. Locoregional-only recurrences are not events.
- A cohort enters a stratum ONLY if its raw patient-level fields express that
  stratum's event. Otherwise it is excluded from that stratum with the reason
  recorded. No imputation, no proxy mapping across endpoint types, ever.
- Pooling: per-cohort results are pooled only within a stratum (pooled
  C-index, patient-level bootstrap stratified by cohort). RFS and DMFS
  results are never pooled together or compared as if one endpoint. TCGA is
  excluded: the committed TCGA analysis is a case-control early-recurrence
  comparison, not a patient-level survival transport.

## Model and data
- Model: identical to the committed transport protocol (run_transport.py):
  penalized Cox (penalizer 0.05) on the METABRIC 70-gene expression panel,
  fit once on the committed discovery data. No re-tuning after seeing any
  harmonized outcome.
- Discovery endpoint remains METABRIC RFS (committed data_cache/clinical.json
  RFS_MONTHS/RFS_STATUS).
- External cohort clinical data are re-pulled with the same documented
  GEOparse pattern as run_transport*.py into data_cache/external/geo/; the
  pulled clinical files are committed before any outcome is computed, so the
  derivation inputs are frozen and inspectable.
- Candidate cohorts and their committed-script endpoint fields:
  GSE2990 (event.rfs/time.rfs), GSE7390 (RFS-comparable per
  results/transport_calibration.json), GSE11121 (t.dmfs/e.dmfs),
  GSE20685 (event_metastasis), GSE25066 (drfs). Eligibility of each cohort
  for each stratum is decided from the pulled raw files under the inclusion
  rule above and recorded either way.

## Report (no performance gate; methodological item)
- Per-stratum, per-cohort: n, events, C-index with 95% bootstrap CI.
- Within-stratum pooled C-index with cohort-stratified bootstrap CI.
- Derivation audit table: cohort, raw fields used, verbatim source
  definition, unit conversion, exclusion reasons.
- Everything is labeled a descriptive audit. No clinical-threshold claim.
- This study does not reopen any committed transport result and changes no
  existing "do not pool" language; it adds the patient-level harmonized
  analysis the verdict asked for.
