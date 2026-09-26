# MEGA27-21 Revival Preregistration (locked 2026-09-26, before any new outcome)

Current audited state (paper/main.tex @ delegated HEAD): five-model METABRIC
benchmark complete, penalized Cox clinical+expression C-index 0.680 best;
deep models trail (honest negative). Cross-cohort transport beat: frozen
70-gene Cox beats GGI in GSE7390 (0.6009 vs 0.5366) and GSE25066
(0.640 vs 0.6025), Veridex-76 (0.6009 vs 0.5757), tumor grade in GSE11121
(0.6979 vs 0.6294); honest non-beats in GSE2990 and GSE20685. R3 mitotic
signature enriched (BH q to 8.8e-5) but failed nested-CV (AUC 0.50) and
TCGA-BRCA replication (p=0.36): bounded negative, pivot required (rule 4).

## Locked gates (declared before outcomes)
- G1 discovery pivot: the next discovery candidate must replicate on an
  untouched external cohort at pre-declared alpha before any claim. Ranked
  candidates (to be ordered by a rule-6 ChatGPT redirection round before
  testing): (a) treatment-era-matched subtype-specific risk re-calibration;
  (b) interaction of panel score with nodal status in GSE20685 (where
  clinical covariates won); (c) a transport-stability map: which genes drive
  cross-platform transport beats vs losses, with a locked permutation
  control.
- G2 benchmark: existing transport beats stand as the benchmark-beat
  evidence; any new model claim must beat the frozen Cox on the same locked
  cohorts with paired bootstrap CIs.
- Judge: >= 10 ChatGPT rounds (weaknesses + additions), verbatim
  docs/JUDGE_ROUNDS.md; paper expanded toward 50+ text-body pages (methods,
  per-cohort audits, negative-round analyses) without padding.
- ISEF archetype: biomarker identify-then-VERIFY end-to-end (the
  Kulviwat-style pattern the user named 4:14:37); replication-first claim
  discipline is the explicit judge-facing strength.
