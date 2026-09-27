# Feasibility check: direct same-patient MammaPrint / Oncotype DX comparison

Date: 2026-09-27. Lane 21 verdict item: "Compare against MammaPrint/Oncotype DX on same patients."

## Data constraint (checked against committed cache)
The committed METABRIC expression cache (data_cache/genes.json) contains exactly
the 70 panel genes, not the full transcriptome. Any signature outside those 70
genes cannot be reconstructed without a new data pull.

## Oncotype DX (21-gene Recurrence Score)
- Cancer genes in cache: 12 of 16 (with aliases: AURKA, CTSV, SCUBE2, CD68 missing).
- Reference genes (ACTB, GAPDH, RPLPO, GUS, TFRC): 0 of 5 in cache.
- Algorithm: proprietary weights and reference normalization; not public.
Verdict: exact reconstruction impossible from the committed cache. Any attempt
would be an unlabeled approximation of a marketed test - not done.

## MammaPrint (70-gene)
- The committed panel overlaps MammaPrint's proliferation cluster, but the full
  published 70-gene list has not been verified against the panel in this repo;
  overlap is partial by construction (panel = PAM50 + core drivers).
- Algorithm: published as correlation to a centroid template, but the exact
  template coefficients are not fully public.
Verdict: at best a labeled approximation; not the marketed test.

## Consequence for the verdict item
A direct same-patient comparison with the marketed scores is NOT executable on
the committed data. What exists instead (already committed): indirect
comparators - GGI (continuous and class), Veridex-76, DLDA-30, tumor grade,
nodal stage - across five external cohorts, with honest non-beats reported.
If this item is to be executed properly it needs (a) a full-transcriptome
METABRIC pull and (b) published approximation implementations labeled as
approximations. Both are out of scope for today's compute budget.
