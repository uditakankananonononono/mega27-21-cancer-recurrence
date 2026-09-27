# Locked descriptive competing-risk CIF audit, 2026-09-28

Before running: use the already committed 1,975 METABRIC patient-level recurrence/death coding from `run_competing_risks.code_competing`, with 800 recurrences, 431 deaths before recurrence and 744 censored expected. Compute the Aalen-Johansen cumulative incidence of recurrence and of death before recurrence at 60 months. Compare the recurrence CIF with a naive Kaplan-Meier event probability treating competing death as ordinary censoring. No new risk model is fitted and no Fine-Gray claim is made; this is a descriptive estimate of competing-event impact.

Report estimates with Aalen-Johansen variance-derived confidence bounds if available; record raw counts and the exact coding/horizon. If input counts differ from the earlier committed table, stop and investigate. Avoid translating the cause-specific Cox into Fine-Gray. This audit is not endpoint harmonization across cohorts or treatment decision guidance.
