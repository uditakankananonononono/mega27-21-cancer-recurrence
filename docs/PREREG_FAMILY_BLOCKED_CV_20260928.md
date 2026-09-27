# Preregistration: gene-family-blocked evaluation of the 70-gene panel (2026-09-28)

Queue item 10 (docs/QUEUE_20260927.md): "Use protein-family-blocked CV -
gene-family-blocked splits for the 70-gene panel." Written before any family
arm was fit. Thresholds and the decision rule below are locked; they will not
move after outcomes.

## Question
Which functional gene families carry the penalized-Cox panel's held-out
discrimination, and which are dispensable or harmful? Context: the fixed
seed-42 holdout already showed removing a 7-gene proliferation subset did
not hurt (results/beyond_proliferation.json, full 0.6464 vs reduced 0.6544).

## Families (fixed a priori, PAM50-standard annotation for the panel core;
pathway annotation for the remaining genes; all 70 genes covered exactly once)
- proliferation_mitotic (21): MKI67 BIRC5 CCNB1 CDC20 CENPF MYBL2 UBE2C CEP55
  KIF2C MELK NDC80 NUF2 ORC6 PTTG1 EXO1 CDC6 RRM2 TYMS UBE2T CCNE1 ANLN
- er_luminal (14): ESR1 PGR FOXA1 GATA3 BCL2 MAPT NAT1 SLC39A6 MLPH CXXC5
  GPR160 TMEM45B BLVRA BAG1
- her2_amplicon (2): ERBB2 GRB7
- basal_myoepithelial (10): KRT5 KRT14 KRT17 EGFR CDH3 FOXC1 MIA SFRP1 PHGDH MMP11
- pi3k_akt_mtor (4): AKT1 PIK3CA MTOR PTEN
- dna_repair_genome_stability (8): BRCA1 BRCA2 ATM CHEK2 PALB2 TP53 RB1 MDM2
- cellcycle_apoptosis_signaling (5): CCND1 BCL2L1 MCL1 STAT3 MYC
- hypoxia_adhesion_detox_other (6): HIF1A VEGFA FGFR4 CDH1 ACTR3B GSTM1
  (GSTM1 is here as detox/metabolism, not in cellcycle_apoptosis_signaling.)

Erratum (post-outcome, cosmetic only): the original text showed "(7)" for
cellcycle_apoptosis_signaling with an inline note routing GSTM1 to the "other"
family; the executed family membership in run_family_blocked_cv.py was always
the 5 genes listed above and the runtime assertion verified all 70 panel genes
covered exactly once. Family membership and the locked rule are unchanged.

## Design
- Data: rebuilt public cBioPortal METABRIC cache (1,975 patients, 800 RFS
  events), same assemble() path as run_beyond_proliferation.py.
- Features: the 4 declared clinical covariates (grade, ER, HER2, positive
  nodes) PLUS the 70-gene expression matrix, exactly as beyond_proliferation's
  "full" arm. Clinical covariates are present in every arm; blocking applies
  to genes only.
- Model: lifelines CoxPHFitter, penalizer 0.05, same as beyond_proliferation.
- Evaluation A (primary): 5-fold patient CV, folds fixed by
  sklearn KFold(n_splits=5, shuffle=True, random_state=42), stratified on
  event. For each family: drop its genes, fit on 4 folds, predict risk on the
  held-out fold, pool all out-of-fold predictions, one C-index per arm.
  Baseline = full 70-gene arm on identical folds. Reported per family:
  pooled out-of-fold C-index and delta vs baseline.
- Evaluation B (control, fixed seed-42 80/20 holdout identical to
  beyond_proliferation): for each family, 100 random size-matched gene sets
  drawn without replacement from the 70 genes (rng seed 20260928), refit,
  held-out C-index; the family's own drop is placed in this null and its
  empirical percentile reported.

## Locked decision rule
- A family is declared LOAD-BEARING if (i) its 5-fold pooled out-of-fold
  C-index drops by > 0.005 vs baseline AND (ii) its fixed-holdout drop exceeds
  the 95th percentile of its size-matched random null.
- A family whose removal IMPROVES either evaluation by > 0.005 is declared
  DISPENSABLE-OR-HARMFUL and reported as such, honestly (consistent with the
  existing proliferation negative).
- Everything else is reported descriptively; no family-level superiority or
  clinical claim is made; this audit does not reopen any prior verdict.
- Failures of fit (nonconvergence) are recorded as NA, never imputed.
