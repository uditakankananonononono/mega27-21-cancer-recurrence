# Fetch log: fresh transport holdout (per PREREG_TRANSPORT_GSE17705_20260928.md)

2026-09-28 IST. Candidate order per the locked prereg: GSE17705 primary,
GSE6532 backup, GSE19615 ineligible (n=115 < 150, recorded pre-fetch).

## GSE17705 - DESIGNATED (first qualifying cohort)

- Source: GEO series matrix, official NCBI FTP
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE17nnn/GSE17705/matrix/GSE17705_series_matrix.txt.gz
  pulled 2026-09-28 ~05:48 IST into data_cache/geo/GSE17705_series_matrix.txt.gz
  (51,521,008 bytes, md5 c68334acd17a3904ea37d928179b28eb).
- Platform: GPL96 (Affymetrix HG-U133A); annotation already cached at
  data_cache/external/geo/GPL96.txt (ID + Gene Symbol columns verified).
- Eligibility against the locked rules, from the series matrix itself:
  1. Untouched: GSE17705 appears nowhere in this repository's history (checked
     by full-repo grep before designation); not on the prereg exclusion list.
  2. Breast cancer primary tumors, 298 samples (>= 150). NATIVE recurrence
     endpoint: "distant relapse (1=dr, 0 censored)" + "event time (years)" =
     DRFS, declared in the cohort's own sample records (71 events >= 40;
     times numeric for all 298, 0.50-16.27 years). Not DFS.
  3. Panel-gene coverage: asserted >= 60/70 inside run_transport17705.py at
     execution (recorded in the result JSON).
- Cohort facts for the record: all 298 ER+, uniformly tamoxifen-treated 5
  years; profiled by two labs (MD Anderson, Jules Bordet) - lab breakdown is
  descriptive only per the prereg; nodal status recorded for 287 samples
  (175 negative / 112 positive / 11 NA) and is the locked clinical-only
  comparator field (grade and age are absent; ER constant).
- Backup GSE6532 not fetched (primary qualified).
