# User verdict archive, 2026-09-27 (lane 21 = "main (12): Breast Cancer Recurrence")

Provenance: inbound WhatsApp wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMEFGRDY4MzY4OTkxNzFEQURGRAA=,
received 2026-09-27 12:02:56 IST, author = user (+918134098571). Full-message body SHA256:
0795f5ff4a86b44e7edd301caa8c9c4ce17abbda034f0a0b29ae7b40ff76a6a6.
Header directive (verbatim, first line of the original message): "IGNORE ABOUT ISEF DELIVERABLES,
IMPROVE PAGE COUNT". Effect: page-length weaknesses fleet-wide are superseded; grow pages with
substantive content; 12-slide storyboard deliverables dropped.

## Verbatim section for this lane
2. main (12): Breast Cancer Recurrence — METABRIC/TCGA Cox
Weaknesses (20) — computational only:

Ceiling is a penalized linear Cox model — no deep-learning win.

DeepSurv, CoxCNN, CoxGCN all underperform clinical+expression Cox.

GCN has the highest seed variance (±0.015).

Locked stability test failed: 3 folds had no stable genes.

Multicohort transport wins are point-estimate only; CIs overlap.

Two cohorts lose (GGI, nodal stage).

Comparator encodings had sign-convention bugs (fixed, but shows fragility).

DLDA-30 metadata column dropped mid-analysis.

Endpoint mismatch: RFS vs DFS across cohorts.

Cross-platform (Affymetrix vs Illumina) without harmonization.

Treatment-era differences between METABRIC and TCGA.

Panel is expression-only; clinical covariates still carry signal (GSE20685).

No competing-risks handling.

No survival calibration (Brier, decision curves).

Tool count below 40-gate.

Dataset count is 9 studies, not 120 (nested records).

Failed discovery arms (low-NPI mitotic contrast) only in appendix.

"Benchmark beat" framing overstates point-estimate wins.

GCN graph is k-NN on expression — not biological graph.

No leave-cohort-out retraining.

Additions (computational):

Add survival calibration and decision-curve analysis.

Add competing-risks models (Fine-Gray).

Harmonize endpoints before pooling.

Add a fully locked prospective holdout with pre-registered genes.

Compare against MammaPrint/Oncotype DX on same patients.

Report per-cohort calibration slopes.

Add "when does deep learning help?" n-threshold analysis.

Pre-register the stability test before re-running.

Add gene-level biological interpretation of top Cox coefficients.

Use protein-family-blocked CV.


## Verbatim cross-cutting themes

Cross-Cutting Computational Themes
Recurring weaknesses:

Dataset/tool count inflation (nested records counted as independent).

Long papers (47-58 pages) — not ISEF-ready.

Negative-heavy narratives that obscure positive contributions.

Single-seed headline numbers.

Ad hoc gates/thresholds rather than theoretically derived.

Homology leakage in random-split benchmarks.

No leave-family-out CV in most projects.

CIs often overlapping — point-estimate wins only.

Universal computational additions:

One primary question per paper.

One locked primary endpoint.

Cluster-level bootstrap CIs everywhere.

Leave-family-out CV as the primary protocol.

12-slide storyboard as the ISEF deliverable.

One-page summary card.
