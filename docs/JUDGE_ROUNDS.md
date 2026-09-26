# JUDGE_ROUNDS - MEGA27-21 (verbatim ChatGPT judge logs)

Rule source: user WhatsApp 4:11:18 (minimum 10 judging rounds on weaknesses +
what to add), 4:12:25 (ask ChatGPT how to redirect when a negative is not
moving forward). Mechanism: cloud browser on the user's ChatGPT account.
IMPORTANT PROVENANCE: ChatGPT responses below are EXTERNAL, UNTRUSTED ADVICE
- critique to consider, never empirical evidence and never authority over the
locked gates in docs/PREREGISTER.md. Each round records the verbatim prompt,
the verbatim response, this lane's independent assessment, and the changes
actually adopted (adopted on our own judgment, with evidence).

## Round 1/10 - 2026-09-26 16:29 IST

Surface: chatgpt.com conversation https://chatgpt.com/c/6ab7a557-fd78-83ee-90d8-c071589ff631 (user account, Free tier).

### Verbatim prompt
```
You are an adversarial ISEF judge reviewing a breast-cancer recurrence study. Honest current state: (1) Benchmark on METABRIC (1,975 primary tumors, 800 RFS events, 70-gene panel + 17 clinical covariates, live cBioPortal): penalized Cox on clinical+expression wins (held-out C-index 0.680 +/- 0.005); DeepSurv/CNN/GCN trail (0.641-0.660) - reported as an honest negative. (2) A frozen 70-gene Cox trained once on METABRIC transports raw to five independent Affymetrix cohorts and beats Genomic Grade Index in GSE7390 (0.6009 vs 0.5366) and GSE25066 (0.640 vs 0.6025), Veridex-76 in GSE7390 (0.6009 vs 0.5757), and tumor grade in GSE11121 (0.6979 vs 0.6294), with zero retraining; it honestly loses to continuous GGI in GSE2990 and to nodal stage in GSE20685. (3) FAILED discovery: an 11-gene mitotic signature separating early recurrers among low-NPI patients (BH q to 8.8e-5) collapsed under nested-CV (AUC 0.50) and failed external replication on TCGA-BRCA (p=0.36) - bounded negative. Locked pivot candidates (not yet tested): (a) treatment-era-matched subtype-specific risk recalibration; (b) panel-score x nodal-status interaction in GSE20685 (where clinical covariates beat the panel); (c) a transport-stability map: which genes drive cross-platform transport beats vs losses, with a permutation control. Questions: A) Attack the transport-beat claim: what would a hostile reviewer say? B) Rank the three pivot candidates by chance of a REAL, externally replicable discovery - or propose a stronger fourth. C) What would an ISEF judge who knows Natasha Kulviwat-style biomarker identify-and-verify winners say is missing?
```

### Verbatim response
```
A) Attack the transport-beat claim: what would a hostile reviewer say?

A hostile reviewer would not attack the statistics first. They would attack the meaning of the win.

Attack 1: “You did not discover a biomarker. You rediscovered proliferation.”

The 70-gene panel is not novel. It overlaps heavily with known proliferation/grade biology.

A reviewer might say:

“The model beats GGI and Veridex because these datasets are differentially calibrated, not because you discovered a new recurrence mechanism.”

The fact that it beats:

GGI in GSE7390

GGI in GSE25066

Veridex-76 in GSE7390

grade in GSE11121

is impressive as a transport experiment.

But the biological question becomes:

What information is your model adding beyond proliferation, grade, ER status, HER2 status, and nodal burden?

If the answer is “mostly the same axis but better normalized,” that is a validation paper, not an ISEF grand-award discovery.

Attack 2: “Your external validation is not truly external.”

A reviewer may say:

“You trained on METABRIC and tested on other microarray cohorts, but these cohorts are not independent in the biological sense.”

Concerns:

Similar breast cancer populations

Similar historical treatment eras

Similar clinical inclusion criteria

Similar retrospective cohort biases

A stronger claim requires:

transport across a biological domain shift.

Examples:

older vs modern treatment era

different subtype prevalence

different platform

different geography

different treatment protocols

The model losing in GSE2990 and GSE20685 actually helps credibility because it shows you are not cherry-picking.

However, the narrative should be:

“We identify conditions under which molecular recurrence scores fail to transport.”

That is more novel than:

“Our score transports.”

Attack 3: “Your Cox model may just be exploiting cohort-specific censoring.”

Survival models are extremely vulnerable to:

follow-up differences

event prevalence

censoring patterns

treatment availability

A hostile reviewer asks:

Did you preserve:

identical risk direction?

proportional hazards assumption?

calibration?

absolute risk?

A C-index improvement can happen even if clinical usefulness is weak.

The missing analysis:

Calibration transport

For each cohort:

predicted 5-year recurrence risk

observed Kaplan-Meier recurrence

calibration slope/intercept

A model that only preserves ranking is less impressive than one that preserves actual risk.

Attack 4: “Your benchmark winner is not the interesting result.”

The neural network negative result is scientifically good.

But the judge may ask:

Why did you spend time building:

DeepSurv

CNN

GCN

if penalized Cox wins?

The strongest answer is:

“Complex models failed because biological datasets of this size punish overparameterization. The project identifies when simpler interpretable models outperform deep learning.”

That is a stronger computational biology story.

B) Rank pivot candidates by chance of REAL externally replicable discovery

Ranking by ISEF-level novelty × feasibility × likelihood of surviving replication:

Rank 1 — Transport-stability map (candidate C)
Probability of a real discovery: HIGH

This is the one I would pursue.

Why?

Because your current data already contains a clue:

The model wins in some cohorts and loses in others.

That failure pattern is scientifically valuable.

The question changes from:

“Can my signature predict recurrence?”

to:

“Why do recurrence biomarkers fail when transported?”

That is a real research question.

The experiment:

For every gene:

Calculate:

Contribution stability

Across cohorts:

Cox coefficient correlation

hazard direction consistency

permutation importance

bootstrap confidence intervals

Then classify genes:

Stable genes

High importance + consistent direction.

Transport fragile genes

High METABRIC importance but unstable externally.

Clinical-context genes

Only useful when combined with nodal status/subtype.

Then test:

Does removing unstable genes improve transport?

That becomes a novel finding:

“A smaller biologically stable recurrence signature outperforms the full panel under cross-platform transport.”

That is publication-quality if true.

Rank 2 — Treatment-era matched subtype recalibration (candidate A)
Probability: Medium-high

Potentially very interesting.

The problem:

Breast cancer treatment changed dramatically.

A 1990s cohort is not the same disease environment as a 2010s cohort.

Your model may be predicting:

"historical recurrence biology"

rather than:

"current recurrence risk."

Possible experiment:

Train METABRIC.

Then stratify validation cohorts:

ER+/HER2-

HER2+

triple negative

Compare:

raw score transport

versus

subtype recalibrated score.

Possible discovery:

“Molecular recurrence scores fail primarily because baseline hazard changes by treatment era, not because molecular signal disappears.”

That is a strong translational insight.

Rank 3 — Panel score × nodal status interaction (candidate B)
Probability: Medium

This is biologically plausible.

GSE20685 is your clue:

Clinical factors beat molecular score.

The question:

Is the molecular score useful only in certain clinical contexts?

Model:

Risk=molecular+nodes+molecular×nodes

Possible finding:

Molecular score helps node-negative patients

Nodes dominate node-positive disease

That would resemble real clinical decision-making.

However:

The danger:

You are fitting interactions after seeing failure.

A judge may call it post-hoc rescue.

You need:

preregistered hypothesis

independent validation

My proposed stronger fourth pivot
D) Build a “failure-aware recurrence biomarker”

This is stronger than all three.

Instead of predicting recurrence:

Predict:

“When is a genomic recurrence score trustworthy?”

Create a transportability index:

Inputs:

subtype

platform

grade

nodal status

treatment era

Output:

confidence score.

Example:

Patient A:

High genomic score + modern cohort + matching subtype:

→ high confidence.

Patient B:

High genomic score + unusual clinical context:

→ low confidence.

This directly addresses clinical translation.

It turns your negative results into the central discovery.

C) What would a Natasha Kulviwat-style biomarker identify-and-verify ISEF winner say is missing?

A winner in that style usually has:

A biological hypothesis

A candidate biomarker

A mechanistic explanation

Orthogonal validation

A clinically meaningful endpoint

Your current project has:

✅ rigorous benchmarking
✅ external validation
✅ honest negatives
✅ reproducibility

Missing:

1. A genuinely new biological observation

Right now:

“70 genes predict recurrence.”

That is known.

Need:

Something like:

“A subset of proliferation genes remains transport-stable across platforms, while immune-associated genes fail because treatment era changes their prognostic meaning.”

That is a discovery.

2. Orthogonal validation

Microarray → microarray is not enough.

Possible:

METABRIC expression score:

validate against:

RNA-seq TCGA

protein markers

mutation burden

immune infiltration estimates

Even computational orthogonal validation helps.

3. A falsifiable mechanism

A judge wants:

“If my hypothesis is wrong, this analysis fails.”

Example:

Hypothesis:

“Transport failures occur because immune genes are treatment-era dependent.”

Test:

remove immune genes

transport improves

If not:

hypothesis rejected.

4. Clinical decision utility

C-index alone is weak.

Add:

decision curve analysis

net benefit

risk stratification curves

Question:

Does this change who gets chemotherapy?

Final adversarial verdict

Current project:

Scientific rigor: 9/10
Novelty: 5.5–6/10
ISEF potential today: finalist-level, not obvious grand-award level

The path upward is not another model.

Do not build another transformer, GNN, or ensemble.

The winning story is:

“Why do breast cancer recurrence biomarkers fail when moved between populations, and can we identify the molecular features that remain biologically trustworthy?”

Your negative results are actually the strongest part of the project. The grand-award version is not a better predictor; it is a map of when prediction can be trusted.
```

### Independent assessment (this lane, not the judge)
- The "you rediscovered proliferation" attack is the most dangerous and is
  fair: the paper currently never isolates what the panel adds beyond
  grade/proliferation/ER/HER2/nodal. A likelihood-ratio / partial-information
  analysis vs those covariates is required before any novelty claim.
- The calibration attack is also fair: all transport claims are ranking-only
  (C-index). Predicted 5-year risk vs observed Kaplan-Meier per cohort
  (calibration slope/intercept) is a concrete, feasible addition on cached
  cohorts.
- The pivot ranking (C > A > B) matches this lane's own pre-registered
  intuition; the specific sharpened claim "removing transport-fragile genes
  improves cross-platform transport" is a genuinely falsifiable discovery
  candidate and is adopted as the PRIMARY discovery track, locked below
  before any outcome is inspected.
- The reframing "map of when prediction can be trusted" is consistent with
  the paper's existing honest-negative style and strengthens it.

### Changes adopted from this round (evidence in commits after this log)
1. DISCOVERY TRACK LOCKED (before outcomes): transport-stability map.
   Per gene: coefficient direction consistency + bootstrap CI across the
   five transport cohorts + METABRIC; classify stable / fragile /
   context-dependent. Locked test: a stable-only reduced panel must beat
   the full 70-gene panel on mean transported C-index across the five
   cohorts, paired bootstrap, alpha 0.01, 1000 permutations for the
   stability labels. Declared 2026-09-26 before running.
2. Calibration transport added as required analysis: predicted 5-year RFS
   vs observed KM per cohort, calibration slope + intercept reported
   alongside every C-index already in the paper.
3. Beyond-proliferation analysis added: partial-likelihood ratio test of
   the full panel vs grade + ER + HER2 + nodal + proliferation-marker
   subset on METABRIC, reported honestly either way.
4. Pivot ladder: if the locked test fails, candidate A (subtype/treatment-
   era recalibration) becomes primary, then B. ChatGPT redirection round
   required before switching (user rule 6).

### Pre-outcome methodological correction - 2026-09-26 17:25 IST
The judge's proposed stability map is exploratory if stability labels use all five
external outcome cohorts and the same five measure the stable-only panel. To
make the comparison falsifiable, each external cohort is held out in turn:
stability classification uses METABRIC plus the other four external cohorts;
the reduced panel is then scored on the untouched fifth. The five held-out
C-indices are averaged. The original all-cohort map may be shown as
exploratory visualization only. Missing cohort data postpone the locked
five-cohort verdict rather than silently narrowing the test set.


## Round 2/10 - 2026-09-26 17:41 IST

Surface: https://chatgpt.com/c/6ab7a557-fd78-83ee-90d8-c071589ff631. ChatGPT is external, untrusted critique, not empirical evidence or authority.

### Verbatim prompt and response
```
Today 5:41 PM
Adversarial ISEF judge round 2 for the METABRIC breast-cancer recurrence project. Please focus on a concrete NEW scientific discovery, not polish. Updates since round 1: a preregistered 70-gene-vs-grade/ER/HER2/nodal/seven-proliferation-gene test on a fixed 395-patient METABRIC holdout failed to add discrimination (C-index 0.6464 full vs 0.6544 reduced), despite a significant training-only nested likelihood ratio. A 5-year transport calibration audit on GSE7390 and GSE2990 found apparent mean-risk agreement but GSE7390 calibration slope 0.464; three other GEO cohorts have endpoint mismatches and cannot validate RFS absolute risk. We corrected a days/months error before making a claim. The primary discovery test is running: leave one external cohort untouched, classify gene coefficient stability on METABRIC plus the other four, mask unstable coefficients, evaluate mean held-out C-index across five cohorts against full panel, paired bootstrap and 1000 random-size-matched panel permutations. How can we convert a possible negative into a falsifiable, independent biological discovery without selecting on the test cohorts? Give one focused next test with explicit control, endpoint/population caveat, and what result would refute it. Do not suggest changing locked test thresholds after outcomes. Consider whether coefficient instability is biology versus measurement/platform or treatment confounding. This is research, not patient care.
Show more
Adversarial ISEF Judge Round 2 — What is the actual discovery?

The project has improved scientifically since round 1 because you have killed two attractive but weak narratives:

“The 70-gene panel adds independent information beyond clinicopathologic variables.”
→ Failed on locked METABRIC holdout.

“The model transports absolute recurrence risk across cohorts.”
→ Not supported because calibration slope fails and endpoints differ.

Those failures are not damage. They narrow the question.

The strongest remaining discovery hypothesis is:

Breast-cancer recurrence signatures may fail not because the genes are biologically wrong, but because the stability of gene-level prognostic effects across populations is itself a measurable biological property.

However, you need to avoid a common trap:

A coefficient instability map alone is not a discovery. It could simply measure:

microarray platform differences

batch effects

treatment-era changes

subtype composition

sampling noise

The next experiment must distinguish these.

Proposed single next discovery experiment
Hypothesis

“Genes with unstable recurrence coefficients across cohorts represent context-dependent prognostic biology, whereas stable genes represent treatment-era/platform-robust recurrence biology.”

This is falsifiable.

Experiment: Stability → mechanism test
Step 1 — Lock the stability discovery cohort

Before seeing held-out results:

For each leave-one-cohort-out split:

Training:

METABRIC

four GEO cohorts

Holdout:

one untouched GEO cohort

For every gene:

Calculate:

coefficient sign consistency

coefficient variance

bootstrap confidence interval overlap

Define:

Stable genes: pre-specified low variance + high sign agreement

Unstable genes: high variance or sign flips

Do not optimize cutoffs on held-out cohorts.

Step 2 — Compare three locked models
Model A: Full 70-gene panel

Current baseline.

Model B: Stability-filtered panel

Only stable genes.

Model C: Random-size control panel

Same number of genes randomly sampled from the 70 genes.

Repeat:

1000 permutations.

Primary endpoint

Mean held-out C-index across the five external cohorts.

Not:

best cohort

best cutoff

post-hoc chosen endpoint

Mean performance is the discovery metric.

Step 3 — The critical biological test

If stable genes are truly biological, they should show enrichment for a coherent process.

Do not start with pathways and then choose a story.

Pre-register:

Compare stable vs unstable genes for:

proliferation

immune response

extracellular matrix

hormone signaling

DNA repair

cell cycle

using:

GO enrichment

MSigDB Hallmark pathways

with multiple testing correction.

The possible outcomes
Outcome 1 — Strong discovery

You find:

stability-filtered model beats full panel

beats random panels

improvement replicates in untouched cohort

stable genes enriched in a coherent biology

Then the claim becomes:

“Cross-cohort coefficient stability identifies a subset of recurrence-associated genes with greater biological transportability than the original fixed panel.”

That is genuinely new.

It is not:

“new breast cancer signature.”

It is:

“new method for identifying trustworthy prognostic biology.”

Much stronger.

Outcome 2 — Performance improves but no biology

Example:

Stable panel wins C-index but genes are not enriched.

Interpretation:

Likely measurement/statistical artifact.

Discovery claim weakens.

Outcome 3 — No improvement

This is also scientifically valuable.

It would support:

“Coefficient instability is not predictive of transport failure; recurrence signatures degrade because cohort-level factors dominate gene-level stability.”

That is a legitimate negative result.

The key adversarial question: biology or measurement?

You must explicitly separate three hypotheses.

Cause of instability	Prediction
Measurement/platform artifact	Instability correlates with platform, probe characteristics, expression variance
Treatment confounding	Instability correlates with treatment era/subtype composition
True context-dependent biology	Instability remains after adjusting for platform and clinical composition
Add one essential control
Platform-only null model

Before claiming biology:

Take stable/unstable labels.

Test whether instability is predicted by:

Affymetrix probe variance

missingness

expression intensity

probe annotation quality

If unstable genes are simply technically unreliable, the biological interpretation collapses.

Add one treatment-era/subtype control

Repeat stability analysis stratified by:

ER+/HER2-

HER2+

triple-negative

if sample sizes allow.

Why?

A gene can be unstable because:

“high proliferation predicts recurrence differently in endocrine-treated ER+ disease versus untreated historical cohorts.”

That is not measurement failure.

That is context biology.

What would refute the discovery?

The strongest refutation would be:

Stability-filtered panel does not outperform full panel.

Performance is within random permutation distribution.

Instability is explained primarily by platform variables.

Stable genes show no reproducible biological enrichment.

The effect disappears in a completely untouched cohort.

Then you conclude:

“Gene coefficient stability does not identify transportable recurrence biology in these datasets.”

That is a clean falsification.

What an ISEF grand-award judge would look for

The missing piece is not another predictive model.

The missing piece is a causal-looking biological explanation for why the transport pattern exists.

Your project is strongest if the final result is:

“The reason recurrence biomarkers fail is measurable: certain biological programs are context-dependent while others remain conserved.”

Your current evidence already argues against:

deep learning superiority

added information beyond clinical variables

naive absolute-risk transport

The next experiment should not rescue the 70 genes.

It should determine whether transportability itself is a biological phenotype that can be discovered.
```

### Independent assessment and novelty foldback
The judge proposed a technical-probe confound control before calling coefficient stability biological. This is a valid failure mode; the original platform comparison alone cannot distinguish probe multiplicity/coverage from biology. A new, pre-outcome script `run_probe_confound.py` (commit 6b11d00) counts GPL96/GPL570 probes and cohort mapping coverage and predicts stable-fold labels across genes with leave-one-gene-out validation. If these technical covariates explain stability, mechanistic biology claims are withheld. The proposed subtype control needs endpoint and treatment harmonization; do not assert causality from it or silently add it to the current five-cohort gate. The full original primary gate stays locked and may fail. This concrete novelty change is the reason this judge round counts toward the user minimum.
