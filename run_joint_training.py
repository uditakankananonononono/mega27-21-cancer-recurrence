"""Pivot 11A: cohort-stratified joint training, prereg 2026-09-28 07:37 IST.
Reuses the committed LOCO loading path by source extraction up to its
analysis marker (no edit to the locked artifact); integrity gate runs there.
Joint 5-fold stratified training; pooled OOF delta vs frozen, B=2000."""
import json
import numpy as np, pandas as pd
from pathlib import Path
from sklearn.model_selection import StratifiedKFold
from lifelines import CoxPHFitter

# --- committed loading path, byte-identical to e1b5ba2 (runs its integrity gate) ---
src = Path('run_loco_retraining.py').read_text()
head = src.split('# ---- LOCO retrains')[0]
exec(compile(head, 'run_loco_retraining.py', 'exec'), globals())

SEED = 20260928
rng = np.random.default_rng(SEED)
out_strata = {}
for stratum in ['HARM-DMFS', 'HARM-RFS']:  # primary first, fixed order
    cohorts = STRATA[stratum]
    X = np.vstack([data[(c, stratum)]['X'] for c in cohorts])
    t = np.concatenate([data[(c, stratum)]['t'] for c in cohorts])
    e = np.concatenate([data[(c, stratum)]['e'] for c in cohorts])
    etaf = np.concatenate([data[(c, stratum)]['eta_frozen'] for c in cohorts])
    clab = np.concatenate([[c] * data[(c, stratum)]['n'] for c in cohorts])
    ylab = np.array([f'{c}_{ev}' for c, ev in zip(clab, e)])
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    eta_joint = np.full(len(t), np.nan)
    sign_maps, coef_folds = [], []
    for k, (tr, te) in enumerate(skf.split(X, ylab)):
        # feasibility: >=5 events per cohort in training pool; never re-split
        for c in cohorts:
            ev_tr = e[tr][clab[tr] == c].sum()
            assert ev_tr >= 5, (stratum, c, k, ev_tr)
        keep = X[tr].var(0) > 0
        df = {f'x{j}': X[tr][:, j] for j in np.where(keep)[0]}
        df['t'], df['e'] = t[tr], e[tr]
        b = CoxPHFitter(penalizer=0.05).fit(pd.DataFrame(df), 't', 'e').params_.values
        eta_joint[te] = X[te][:, keep] @ b
        coef_folds.append((keep, b))
        sign_maps.append(pd.Series(b, index=np.where(keep)[0]))
    # descriptive sign agreement across folds (genes fit in ALL folds)
    common = set(sign_maps[0].index)
    for sm in sign_maps[1:]:
        common &= set(sm.index)
    n_sign_agree = int(sum(len({int(np.sign(sm.loc[j])) for sm in sign_maps}) == 1 for j in common))
    # pooled OOF readout
    c_joint = float(cidx(t, -eta_joint, e))
    c_frozen = float(cidx(t, -etaf, e))
    delta = c_joint - c_frozen
    n = len(t)
    dlt, skipped = [], 0
    for _ in range(2000):
        bi = rng.integers(0, n, n)
        cj = cidx(t[bi], -eta_joint[bi], e[bi])
        cf = cidx(t[bi], -etaf[bi], e[bi])
        if np.isnan(cj) or np.isnan(cf):
            skipped += 1; continue
        dlt.append(cj - cf)
    dlo, dhi = [float(x) for x in np.percentile(dlt, [2.5, 97.5])]
    # descriptive per-cohort OOF deltas, B=500, same stream
    per_cohort = []
    for c in cohorts:
        m = clab == c
        cj_c = float(cidx(t[m], -eta_joint[m], e[m]))
        cf_c = float(cidx(t[m], -etaf[m], e[m]))
        dd = []
        nn = int(m.sum())
        tm, em, jm, fm = t[m], e[m], eta_joint[m], etaf[m]
        for _ in range(500):
            bi = rng.integers(0, nn, nn)
            a = cidx(tm[bi], -jm[bi], em[bi]); b_ = cidx(tm[bi], -fm[bi], em[bi])
            if not (np.isnan(a) or np.isnan(b_)):
                dd.append(a - b_)
        lo, hi = [float(x) for x in np.percentile(dd, [2.5, 97.5])]
        per_cohort.append({'cohort': c, 'n': nn, 'c_joint': round(cj_c, 4),
                           'c_frozen': round(cf_c, 4), 'delta': round(cj_c - cf_c, 4),
                           'delta_ci95': [round(lo, 4), round(hi, 4)]})
    out_strata[stratum] = {
        'cohorts': cohorts, 'n_patients': n, 'events': int(e.sum()),
        'pooled_oof_c_joint': round(c_joint, 4), 'pooled_oof_c_frozen': round(c_frozen, 4),
        'pooled_oof_delta': round(delta, 4), 'delta_ci95': [round(dlo, 4), round(dhi, 4)],
        'boot_skipped': skipped, 'positive': bool(dlo > 0),
        'per_cohort_descriptive': per_cohort,
        'sign_agreement_descriptive': {'genes_fit_all_folds': len(common),
                                       'sign_consistent_across_folds': n_sign_agree}}
    print(stratum, 'joint', round(c_joint, 4), 'frozen', round(c_frozen, 4),
          'delta', round(delta, 4), [round(dlo, 4), round(dhi, 4)],
          'POSITIVE' if dlo > 0 else 'measured negative', flush=True)

out = {'prereg': 'docs/PREREG_JOINT_TRAINING_20260928.md', 'seed': SEED,
       'partition': 'StratifiedKFold(5, shuffle, seed) on cohort-x-event labels; one fixed partition per stratum',
       'bootstrap': 'B=2000 paired patient bootstrap on pooled OOF scores (no refit inside bootstrap); per-cohort descriptive B=500; one rng stream, HARM-DMFS then HARM-RFS',
       'integrity': 'committed LOCO loading path reused by source extraction; frozen per-cohort record reproduced to 4 decimals before any joint fit',
       'strata': out_strata,
       'limits': 'Joint training sees every cohort in training: a pooled-joint-method result, NOT new-cohort transport; does not repair the item-11 LOCO negative; one fixed partition per stratum.'}
json.dump(out, open('results/joint_training.json', 'w'), indent=1)
print('wrote results/joint_training.json', flush=True)
