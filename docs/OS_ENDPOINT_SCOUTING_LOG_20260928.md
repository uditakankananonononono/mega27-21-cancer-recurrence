# OS-endpoint scouting log - 2026-09-28 (PREREG_OS_ENDPOINT_SCOUTING_20260928.md)

| Cohort | Native OS fields | n | OS-complete | Deaths | Gate |
|---|---|---|---|---|---|
| GSE7390 | t.os / e.os (TRANSBIG extract; time units to be verified at designation - stored values suggest days, median 4561.5) | 198 | 198 | 56 | QUALIFIES (>=30 deaths) |
| GSE20685 | follow_up_duration (years) / event_death | 327 | 327 | 83 | QUALIFIES (>=30 deaths) |
| GSE11121 | none (t.dmfs/e.dmfs only) | 200 | - | - | no OS |
| GSE2990 | none (RFS/DMFS only) | 189 | - | - | no OS |

Facts read from the committed clinical extracts (data_cache/external/geo/
clinical_*.csv, assembled from the GEO SOFT characteristics during cohort
building; "key: value" prefixes stripped for counting). Neutral eligibility
facts only; no test statistics computed.

**SCOUTING OUTCOME: POSITIVE.** Two cohorts (GSE7390, GSE20685) carry a
native OS endpoint with >= 30 deaths. The locked sufficiency gate (>= 2
cohorts) is met. The OS-endpoint fresh-transport question proceeds under its
own NEW designated-cohort prereg.
