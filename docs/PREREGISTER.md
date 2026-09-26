# MEGA27-21 Revival Preregistration (locked 2026-09-26, before any new outcome)

Current audited state (paper/main.tex @ delegated HEAD): five-model METABRIC
benchmark complete, penalized Cox clinical+expression C-index 0.680 best;
deep models trail (honest negative). Cross-cohort transport beat: frozen
70-gene Cox beats GGI in GSE7390 (0.6009 vs 0.5366) and GSE25066
(0.640 vs 0.6025), Veridex-76 (0.6009 vs 0.5757), tumor grade in GSE11121
(0.6979 vs 0.6294); honest non-beats in GSE2990 and GSE20685. R3 mitotic
signature enriched (BH q to 8.8e-5) but failed nested-CV (AUC 0.50) and
TCGA-BRCA replication (p=0.36): bounded negative, pivot required (rule 4).

## Provenance of this document (honest attribution)
- USER STANDING RULES (verbatim, WhatsApp channel history): 4:11:18 (complete
  all projects except deleted ones; ask CHATGPT for ideas/redirection; minimum
  10 judging rounds on weaknesses/additions; never count a negative as a
  result), 4:11:49 (each project beats benchmarks - improve until it does -
  and produces an actual new discovery), 4:12:25 (ask ChatGPT how to redirect
  when a negative is not moving forward), 4:14:37 (take inspiration from
  previous ISEF winners, e.g. Natasha Kulviwat).
- RESEARCHER-LOCKED METHODOLOGICAL CHOICES (this agent, 2026-09-26, locked
  before inspecting new outcomes): every numeric threshold, alpha level,
  seed count, pivot ladder, gate name and scope framing below. These are the
  lane's own preregistration decisions, NOT user-specified values; they exist
  so results cannot be fished past moving goalposts. Pre-existing gates
  declared by earlier builders in repo history (e.g. the RMSD < 2.0 A redock
  gate already in this repo's README) are inherited, not invented here.

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

## Judge round 1 transport-stability correction (2026-09-26 17:25 IST, pre-outcome)
Any reduced-panel assessment on the same cohorts used to select its stable
genes would be outcome leakage. Use five leave-one-cohort-out folds: on each,
stability is fit on METABRIC and the other four cohorts; the held-out cohort
contributes only its independent C-index. Keep all five required for the
primary verdict. An all-five-cohort stability map is descriptive, not a
validation result. The decision threshold and permutation count recorded in
JUDGE_ROUNDS.md remain unchanged.

## Judge round 3 secondary biological program falsifier (locked before stability result)
The previous judge proposed program-level coherence, not another trained score.
The new `run_program_coherence.py` uses five already-locked stable-gene sets,
**only after** the primary result; it cannot change the primary transport gate.
Within the original 70-gene panel (not the human genome), use MSigDB Hallmark
v2025.1.Hs (official Broad symbols GMT; SHA logged at run time), filter sets
to >=3 panel members, and test >=3 overlap via hypergeometric upper tail,
BH q<.05 per fold. A secondary coherence claim requires one Hallmark program
significant in >=4/5 folds AND the maximum repeated-fold statistic to exceed
1000 common gene-label permutations stratified by absolute METABRIC beta
quartile (seed 3), empirical p<.01. A common permutation preserves fold
intersection and beta-scale structure; using the *maximum* protects against
choosing a pathway after inspection. This is conditional on a panel already
enriched for breast-cancer programs, so a negative remains a negative and a
positive is only program association, not causal mechanism. Check the probe
multiplicity audit separately; biology vs technical confounding remains open.
Source: https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/
