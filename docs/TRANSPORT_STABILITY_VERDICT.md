# Locked transport-stability verdict, run 2026-09-27

The five-cohort leave-one-out gate **fails**, exactly under the pre-outcome
zero-gene rule in `run_transport_stability.py`. Three held-out folds have zero
stable genes: GSE11121, GSE2990, GSE25066. GSE7390 has 3, with stable
C-index 0.5772 versus the full panel's 0.5993; GSE20685 has 2, with stable
C-index 0.5706 versus full 0.6264. Consequently no mean delta, paired
bootstrap p, or gene-label permutation p is defined; none is imputed. This
cannot be reframed as a biological mechanism or discovery. Five-cohort
endpoint/population differences remain a key interpretability limit.

The independently preregistered technical audit finds no gene stable in
three folds, so binary technical LOO AUROC is undefined; probe-count
Spearman tests are non-significant (GPL96 rho -0.127, p .296; GPL570 rho
-0.181, p .134; mapped-cohort count rho .092, p .447), which does **not**
exclude technical confounding. Cohort z-score preprocessing destroys raw
intensity/noise information. The Hallmark program-coherence test fails by
its own zero-fold rule, with no enrichment/permutation result or biological
claim. Exact JSON outputs live in `results/transport_stability.json`,
`results/probe_confound.json`, and `results/program_coherence.json`; all
scripts were committed before the primary result. The prior frozen-panel
external comparator point-estimate beats are independent of this failed
*new* discovery test and are not upgraded to statistical significance.

Next work is a newly locked, independent discovery pivot, not retroactive
relaxation of the stability threshold. Existing G1 remains failed.
