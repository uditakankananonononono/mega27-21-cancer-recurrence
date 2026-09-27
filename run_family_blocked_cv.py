"""Queue item 10: gene-family-blocked evaluation of the 70-gene panel.
Preregistered in docs/PREREG_FAMILY_BLOCKED_CV_20260928.md before any family
arm was fit. Resumable: checkpoints each finished unit to
results/family_blocked_cv_ckpt.json and stops before the wall-clock budget;
re-invoke until it writes results/family_blocked_cv.json. No clinical claim.
"""
import json, os, sys, time
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.model_selection import StratifiedKFold, train_test_split
sys.path.insert(0, 'src')
from recurscan.data.dataset import assemble

BUDGET = float(os.environ.get('FBCV_BUDGET', '90'))
CKPT = 'results/family_blocked_cv_ckpt.json'
FINAL = 'results/family_blocked_cv.json'

FAM = {
 'proliferation_mitotic': ['MKI67','BIRC5','CCNB1','CDC20','CENPF','MYBL2','UBE2C','CEP55','KIF2C','MELK','NDC80','NUF2','ORC6','PTTG1','EXO1','CDC6','RRM2','TYMS','UBE2T','CCNE1','ANLN'],
 'er_luminal': ['ESR1','PGR','FOXA1','GATA3','BCL2','MAPT','NAT1','SLC39A6','MLPH','CXXC5','GPR160','TMEM45B','BLVRA','BAG1'],
 'her2_amplicon': ['ERBB2','GRB7'],
 'basal_myoepithelial': ['KRT5','KRT14','KRT17','EGFR','CDH3','FOXC1','MIA','SFRP1','PHGDH','MMP11'],
 'pi3k_akt_mtor': ['AKT1','PIK3CA','MTOR','PTEN'],
 'dna_repair_genome_stability': ['BRCA1','BRCA2','ATM','CHEK2','PALB2','TP53','RB1','MDM2'],
 'cellcycle_apoptosis_signaling': ['CCND1','BCL2L1','MCL1','STAT3','MYC'],
 'hypoxia_adhesion_detox_other': ['HIF1A','VEGFA','FGFR4','CDH1','ACTR3B','GSTM1'],
}
CLIN = ['GRADE','ER_STATUS=Positive','HER2_STATUS=Positive',
        'LYMPH_NODES_EXAMINED_POSITIVE']

t0 = time.time()
d = assemble()
assert sorted(g for v in FAM.values() for g in v) == sorted(d.gene_names)
clin_idx = [d.clin_names.index(x) for x in CLIN]
gpos = {g: j for j, g in enumerate(d.gene_names)}
T = d.time.astype(float); E = d.event.astype(int)

def design(drop=()):
    keep = [j for g, j in gpos.items() if g not in drop]
    return np.column_stack([d.X_clin[:, clin_idx], d.X_expr[:, keep]])

def frame(arr, idx):
    f = pd.DataFrame(arr[idx], columns=[f'x{i}' for i in range(arr.shape[1])])
    f['t'] = T[idx]; f['e'] = E[idx]
    return f

def fit_predict(arr, tr, te):
    m = CoxPHFitter(penalizer=.05)
    m.fit(frame(arr, tr), 't', 'e')
    return m.predict_partial_hazard(frame(arr, te).drop(columns=['t','e'])).values.reshape(-1)

ck = json.load(open(CKPT)) if os.path.exists(CKPT) else {'done': {}}
def save():
    json.dump(ck, open(CKPT, 'w'))

def unit_key(arm, kind, extra=''):
    return f'{arm}|{kind}|{extra}'

arms = {'full_baseline': []} | {k: v for k, v in sorted(FAM.items())}
skf_splits = list(StratifiedKFold(n_splits=5, shuffle=True, random_state=42).split(np.zeros(len(T)), E))
tr, te = train_test_split(np.arange(len(T)), test_size=.2, random_state=42, stratify=E)

def draw_set(k, i):
    rng = np.random.default_rng(20260928 + i)  # fixed per-draw seed, order-independent
    return list(rng.choice(list(d.gene_names), size=k, replace=False))

for arm, genes in arms.items():
    if time.time() - t0 > BUDGET: break
    drop = set(genes)
    arr = design(drop)
    # unit 1: 5-fold pooled out-of-fold C-index
    k = unit_key(arm, 'cv5')
    if k not in ck['done']:
        risk = np.full(len(T), np.nan)
        for ftr, fte in skf_splits:
            risk[fte] = fit_predict(arr, ftr, fte)
        ck['done'][k] = float(concordance_index(T, -risk, E))
        save(); print(f'{k} -> {ck["done"][k]:.4f}', flush=True)
        if time.time() - t0 > BUDGET: break
    # unit 2: fixed holdout C-index
    k = unit_key(arm, 'holdout')
    if k not in ck['done']:
        ck['done'][k] = float(concordance_index(T[te], -fit_predict(arr, tr, te), E[te]))
        save(); print(f'{k} -> {ck["done"][k]:.4f}', flush=True)
        if time.time() - t0 > BUDGET: break
    # units 3..102: permutation draws (family arms only), batches of 10
    if drop:
        for batch in range(10):
            k = unit_key(arm, 'perm', str(batch))
            if k not in ck['done']:
                vals = []
                for i in range(batch * 10, batch * 10 + 10):
                    rarr = design(set(draw_set(len(drop), i)))
                    vals.append(float(concordance_index(T[te], -fit_predict(rarr, tr, te), E[te])))
                ck['done'][k] = vals
                save(); print(f'{k} done', flush=True)
            if time.time() - t0 > BUDGET: break

# finalize when all units done
need = []
for arm, genes in arms.items():
    need += [unit_key(arm, 'cv5'), unit_key(arm, 'holdout')]
    if genes:
        need += [unit_key(arm, 'perm', str(b)) for b in range(10)]
if all(k in ck['done'] for k in need):
    base_cv = ck['done'][unit_key('full_baseline', 'cv5')]
    base_ho = ck['done'][unit_key('full_baseline', 'holdout')]
    out = {'prereg': 'docs/PREREG_FAMILY_BLOCKED_CV_20260928.md',
           'n': int(len(T)), 'events': int(E.sum()),
           'cv_folds': 5, 'cv_seed': 42, 'holdout_seed': 42,
           'perm_seed_base': 20260928, 'n_perm': 100, 'penalizer': 0.05,
           'families': FAM,
           'baseline': {'cv5_pooled_cindex': base_cv, 'holdout_cindex': base_ho},
           'arms': {}}
    for arm, genes in sorted(FAM.items()):
        cv = ck['done'][unit_key(arm, 'cv5')]
        ho = ck['done'][unit_key(arm, 'holdout')]
        null = [v for b in range(10) for v in ck['done'][unit_key(arm, 'perm', str(b))]]
        null_drops = [base_ho - v for v in null]
        obs_drop = base_ho - ho
        out['arms'][arm] = {
            'n_genes_dropped': len(genes),
            'cv5_pooled_cindex': cv,
            'cv5_drop_vs_baseline': base_cv - cv,
            'holdout_cindex': ho,
            'holdout_drop_vs_baseline': obs_drop,
            'null_drop_mean': float(np.mean(null_drops)),
            'null_drop_p95': float(np.percentile(null_drops, 95)),
            'observed_drop_percentile_vs_null': float(100 * np.mean([nd <= obs_drop for nd in null_drops])),
            'load_bearing_locked_rule': bool((base_cv - cv) > 0.005 and obs_drop > float(np.percentile(null_drops, 95))),
            'dispensable_or_harmful_locked_rule': bool((base_cv - cv) < -0.005 or obs_drop < -0.005)}
    json.dump(out, open(FINAL, 'w'), indent=1)
    print('FINAL written:', FINAL)
else:
    missing = [k for k in need if k not in ck['done']]
    print(f'checkpointed; {len(missing)} units remain')
