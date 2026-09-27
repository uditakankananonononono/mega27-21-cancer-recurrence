# Supplementary Gemini direction consult, 2026-09-27 09:35 IST

Conversation: https://gemini.google.com/app/705379adce23b48a. This is an
external, untrusted suggestion, NOT empirical evidence and NOT a counted
ChatGPT judge round. The old locked stability gate remains failed. The paper
text was pasted, not uploaded.

## Verbatim prompt
```
Supplementary Gemini direction consult for a breast-cancer recurrence computational study. The current paper text and locked-test failure are pasted below. Please propose ONE new, testable discovery direction *within breast-cancer recurrence* based on already-available METABRIC+five GEO cohorts that can lead to an independent verified result, not rebranding the failed stability test. Conditions: no clinical decisions, no causal claim from retrospective cohorts, do not pool RFS, DMFS, DRFS as equivalent for absolute calibration; do not move prior preregistered thresholds. Suggest a public independent benchmark or untouched validation split, specific biology hypothesis, explicit negative control/falsifier and minimal 2-core CPU plan. Be alert to leakage and sample-size limits. Keep the successful benchmark point estimates separate from discovery. This is NOT a ChatGPT round and its answer will not count toward the mandated ten ChatGPT judge rounds.

CURRENT PAPER TEXT:
\subsection{Cross-cohort transportability: benchmark beat}
A penalized Cox model (70-gene panel only), trained once on METABRIC, was
transported raw to five independent Affymetrix cohorts and compared with each
cohort's own published signature calls or clinical comparators. It
outperforms the Genomic Grade Index in two cohorts
($0.6009$ vs $0.5366$ in GSE7390; $0.640$ vs $0.6025$ in GSE25066),
Veridex-76 in GSE7390 ($0.6009$ vs $0.5757$), and tumor grade in the
node-negative GSE11121 cohort ($0.6979$ vs $0.6294$), with zero retraining
across platforms. Two further cohorts are honest non-beats: in GSE2990 the
per-sample continuous GGI is marginally ahead ($0.6651$ vs $0.6559$), and in
GSE20685 nodal stage leads ($0.6996$ vs $0.6264$), though our score still
exceeds $0.6$ with its bootstrap interval excluding $0.5$. Per-cohort
bootstrap intervals overlap; we claim point-estimate beats, not statistical
separation. Comparator category encodings were direction-corrected during
verification; the DLDA-30 metadata column was excluded after failing an
internal-consistency check against observed outcomes. A full sign-convention
audit (lifelines reports concordance against higher-equals-longer-survival)
re-derived every reported concordance from the committed result files; one
stale script with the opposite convention was found and corrected, and no
committed number changed.
\begin{figure}
\centering
\includegraphics[width=0.95\linewidth]{figures/fig5_transport.png}
\caption{Transported model versus each cohort's own published comparator
across five independent cohorts; numbers from \texttt{results/transport*.json}.}
\label{fig:transport}
\end{figure}

\subsection{Locked stability test: negative boundary}
A five-cohort leave-one-cohort-out test asked whether coefficient-sign-stable
genes explain cross-platform transport without selecting on the test cohort.
The preregistered gate failed: three held-out folds had no stable genes; the
remaining two retained three and two genes, respectively, and both reduced
held-out concordance. The prescribed paired bootstrap and label-permutation
comparisons are undefined under that zero-gene rule, not silently replaced
with a smaller analysis. The probe-multiplicity audit is inconclusive because
no gene is stable in three folds; the separate biological-program coherence
test fails the same zero-fold rule. Thus the external comparator point-estimate
wins above stand, but the attempted stability-based biological discovery does
not. Full fold values and technical audit outputs are in the repository.



LOCKED FAILURE/AUDIT:
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

```

## Verbatim visible response
```
Gemini said
Independent Discovery Pivot Proposal

Rather than attempting to rescue or reframe the failed cross-platform coefficient stability test, this proposal introduces a distinct hypothesis within breast cancer recurrence dynamics.

It leverages the existing METABRIC dataset (discovery/training) and five GEO validation cohorts without altering preregistered thresholds, pooling non-equivalent survival endpoints, or asserting causal clinical utility.

1. Biological Hypothesis & Specific Target Mechanism
Mechanistic pathways influencing breast cancer recurrence and dormancy. Source: PMC - NIH
Hypothesis: Early vs. Late Recurrence Subtype-Specific Disruption Ratio (SDR)

Recurrence timing (early <5 years vs. late ≥5 years) in estrogen receptor-positive (ER
+
) breast cancer is governed by distinct microenvironmental and proliferation programs:

Early recurrence is primarily driven by high baseline cell-cycle activity (e.g., MKI67, CCNE1, AURKA).

Late recurrence (dormancy escape) is mediated by stromal remodeling, extracellular matrix turnover, and low-grade inflammatory signaling (e.g., TGFB1, COL1A1, STAT3).

Instead of scoring continuous linear risks across all recurrence events simultaneously, we hypothesize that the ratio of proliferation signature score to stromal/dormancy score (SDR=
Stromal
Proliferation
	​

) stratifies patients into mutually exclusive hazards for early versus late relapse within ER
+
 patients, irrespective of platform noise.

2. Independent Benchmark / Untouched Validation Plan

Discovery Training Set: METABRIC ER
+
 subset (N≈1,500).

Validation Benchmark Set: GEO cohorts with explicit, unpooled recurrence timing metrics (specifically GSE11121 and GSE7390 for recurrence-free survival [RFS], keeping distant metastasis-free survival [DMFS] cohorts strictly separate to avoid absolute calibration bias).

Untouched Validation Split: Lock a 30% stratified random holdout within METABRIC and hold out GSE25066 entirely as an un-accessed final verification cohort.

3. Explicit Negative Control & Falsifier
Negative Controls:

Permuted Gene Sets: Generate 1,000 random gene-set pairs matched for size and expression abundance from the GPL96/GPL570 platforms to calculate pseudo-SDR ratios.

ER
−
 Cohort Control: Apply the SDR model to ER
−
 patients within METABRIC and GEO cohorts. Because late recurrence is rare in ER
−
 disease, the late-hazard stratification must fail.

Explicit Falsification Rule:

The discovery hypothesis is declared falsified if:

The SDR ratio fails to demonstrate a non-zero time-varying hazard ratio (HR(t) crossing 1.0 between years 0–5 vs. years 5–10) in both METABRIC holdout and held-out validation cohorts.

The real SDR empirical performance falls within the 95th percentile distribution of the permuted gene-set negative controls.

4. Minimal 2-Core CPU Execution Plan

To respect hardware limitations and avoid runtime issues across large microarray datasets:

Pre-computation & Subsampling:

Extract pre-calculated signature scores for Proliferation (S
prolif
	​

) and Stromal (S
strom
	​

) modules using z-score medians per platform rather than full matrix inversions.

Parallelized Cox Time-Varying Models:

Use Python (lifelines) or R (survival) on 2 CPU cores using standard process pools (multiprocessing.Pool(processes=2)).

Fit piecewise Cox proportional hazard models split at t=5 years.

Execution Time:

Feature reduction & scoring: ∼45 seconds.

Piecewise Cox fitting + 1,000 permutations across 2 cores: ∼4.5 minutes.

5. Separation of Benchmark Point Estimates vs. Discovery

Frozen Transport Benchmark (Prior Results): Remains untouched as point-estimate wins without statistical separation (e.g., GSE11121 C-index 0.6979 vs 0.6294). The stability test failure (G1) stays logged as failed.

New Discovery Metric: Evaluated purely on time-dependent hazard ratios (HR
<5yr
	​

 vs HR
≥5yr
	​

) within ER
+
 subgroups, completely separated from previous frozen-panel C-index comparisons.

Would you like the minimal Python 2-core script structure for fitting the time-varying SDR Cox model?
Yes




Flash
```

## Independent assessment, before any new analysis
Gemini suggests an ER-positive early-versus-late recurrence timing test using
proliferation and stromal/dormancy modules. It incorrectly calls GSE11121 a
recurrence-free endpoint (our dataset uses DMFS), and calls GSE25066 untouched
although it participated in the earlier stability test. Those claims must not
be used. Its specific stromal genes TGFB1/COL1A1/STAT3 may not be among the
70 mapped panel genes, so source checks are required. Neither its proposed
hazard crossing nor its 95th-percentile random-set cutoff is adopted as a
post-hoc claim. Candidate next step: independently verify gene coverage, ER
stratification and late-event counts; then date-lock a feasible primary test
before running. No biological discovery yet.
