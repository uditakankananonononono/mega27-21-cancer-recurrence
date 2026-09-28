# OS transport fetch/designation log - 2026-09-28 (PREREG_OS_TRANSPORT_20260928.md, 48b72e4)

All facts below verified verbatim in the committed family SOFT extracts
(data_cache/external/geo/GSE7390_family.soft.gz, GSE20685_family.soft.gz)
BEFORE any OS scoring. No outcome quantities computed here.

## GSE7390 (TRANSBIG, GPL96, 198 patients)
- OS fields: `t.os` / `e.os` present per sample (first sample t.os=937, e.os=1;
  stored values are DAYS - t.rfs=723/t.os=937 on the same scale; cohort median
  t.os 4561.5). Unit recorded as stored-days; concordance is unit-invariant.
- Series summary verbatim: "198 N- systemically untreated patients".
  Therefore HORMONE_THERAPY=YES=0 and CHEMOTHERAPY=YES=0 for ALL samples
  (series-recorded constants, set not imputed).
- Frozen-model clinical mappings (verbatim from characteristics):
  AGE_AT_DIAGNOSIS<-age (years); TUMOR_SIZE<-size (stored in CM; converted
  x10 to the model's millimetre scale - declared conversion);
  GRADE<-grade; LYMPH_NODES_EXAMINED_POSITIVE<-node (0 for all: series-recorded
  node-negative cohort); ER_STATUS=Positive<-er (1->1.0, 0->0.0);
  NPI<-NPI (present verbatim).
- Imputed to training mean (absent from series records): PR_STATUS=Positive,
  HER2_STATUS=Positive, the six CLAUDIN_SUBTYPE dummies, RADIO_THERAPY=YES.
- Comparator covariates (prereg-locked set age/size/node/grade/er): node is
  CONSTANT 0 (node-negative cohort) and is dropped as degenerate; comparator
  = in-cohort Cox on age, size(cm, in-cohort scale), grade, er.

## GSE20685 (Kao 2011, GPL570, 327 patients)
- OS fields: `event_death` / `follow_up_duration (years)` present per sample
  (first sample event_death=0, follow-up 7.0 years).
- Frozen-model clinical mappings (verbatim): AGE_AT_DIAGNOSIS<-age at
  diagnosis; CHEMOTHERAPY=YES<-adjuvant_chemotherapy (yes->1.0, no->0.0).
- Imputed to training mean (absent or non-mappable: t_stage/n_stage are
  ordinal stage categories, not size-in-mm or node counts): TUMOR_SIZE,
  GRADE, NPI, LYMPH_NODES_EXAMINED_POSITIVE, ER_STATUS=Positive,
  PR_STATUS=Positive, HER2_STATUS=Positive, the six CLAUDIN_SUBTYPE dummies,
  HORMONE_THERAPY=YES, RADIO_THERAPY=YES.
- Comparator covariates (prereg-locked): age at diagnosis, t_stage, n_stage
  (ordinal numeric, in-cohort scale).

## Common
- Expression: streamed from the family SOFTs; probes mapped to the 70 panel
  genes via the committed GPL96/GPL570 tables (first symbol before " ///",
  per-gene mean across probes), identical machinery to the prior locked
  cohort runs. The >=60/70 mapped-gene gate from prior runs is re-asserted.
- Brier not computed (frozen model carries no baseline hazard).
