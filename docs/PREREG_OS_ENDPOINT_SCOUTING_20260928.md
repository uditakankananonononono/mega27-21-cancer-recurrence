# Locked preregistration - Item 5 pivot A: OS-endpoint scouting (2026-09-28)

**Standing:** MEGA27 lane 21 item 5 (fresh-cohort transport), pivot per the
adjudicated order 2026-09-28 07:22 IST: "(3) 5A OS-endpoint scouting then the
designated-cohort question if a cohort clears scouting." New question, new
dated prereg; the failed claim (DMFS-endpoint fresh transport, measured NOT
SUPPORTED in GSE17705) is never re-tuned.

## Scouting question
Which of the available platform-matched GEO breast-cancer cohorts carry a
native overall-survival (OS) endpoint with enough deaths to support a locked
OS-endpoint fresh-cohort question?

## Candidate set (fixed before inspection)
- GSE7390 (GPL96), GSE11121 (GPL96), GSE2990 (GPL96), GSE20685 (GPL570).
- GSE17705 is excluded: measured to carry DMFS only (the completed
  fresh-transport run), so it cannot serve an OS-endpoint question.

## Facts to verify per cohort (neutral eligibility facts, no test statistics)
1. Sample count in the series matrix.
2. Presence of OS fields in `!Sample_characteristics_ch1` (OS time and OS
   event / vital status keys).
3. Native OS event count (deaths) when parseable from the header block.

## Sufficiency gate (locked before inspection)
- **>= 2 cohorts** with a native OS endpoint and **>= 30 deaths** each:
  scouting SUCCEEDS; the OS-endpoint question proceeds under its own NEW
  designated-cohort prereg (cohorts designated there, not here).
- **0-1 qualifying cohorts:** the OS-endpoint question is exhausted on
  available platform-matched cohorts; the scouting null is the documented
  exhaustion proof and item 5's pivot family closes.

## Scouting null handling
A null is recorded as the scouting outcome in the queue and the paper as the
documented exhaustion proof for this pivot family. Success is not
pre-defined beyond the eligibility gate above.
