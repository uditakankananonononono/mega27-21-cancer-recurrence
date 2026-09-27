# Locked sample-size sensitivity, 2026-09-28 (lane 21 item 8)

Before running: ask whether DeepSurv (MLP) gains relative to a penalized full Cox model as the training sample grows on METABRIC. This is a within-cohort sensitivity study, not a claim of biological benefit or a new independent validation. No threshold is tuned to the result.

- Use the original stratified 25% test split for seeds 0-4 from `benchmark.stratified_event_split`, fixed across all training sizes within each seed. Select nested stratified subsets of the remaining 1,481 training patients at n=200, 400, 800, 1,200, and all 1,481. Selection is deterministic per seed and without replacement; no test patient enters training or preprocessing.
- Same 17 clinical and 70 expression features and same DeepSurv architecture, optimizer, 300 epochs, learning rate 0.001, seed and model training routines as the original benchmark. Cox full model uses penalizer 0.1, as in original. Standardization fits training subset only. No model selection on the test endpoint.
- Metric: test Harrell C-index per seed and size, DeepSurv minus Cox on the same test patients. Report all 25 paired comparisons, mean and seed range by size. Because seeds share patients, do not treat the 5 points as independent trials or infer a causal scaling law.
- Decision description: empirical advantage at a size only if the mean paired DeepSurv-minus-Cox C-index is positive at that size. This is descriptive, not a p-value gate or a claim of superior transport. No interpolation to an unseen size, no post-hoc picking the best size.
- Result checkpoints written after each seed/size; push small JSON regularly. Do not overwrite the original five-seed benchmark.
