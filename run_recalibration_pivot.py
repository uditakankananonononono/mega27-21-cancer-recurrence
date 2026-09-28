"""Pivot 6A: cohort-anchored slope recalibration, prereg 2026-09-28 07:23 IST.
Integrity gate: committed full-cohort slopes reproduce to 4 decimals before
any split is scored. Alpha fit on each cohort's own training half only;
held-out recalibrated slope with held-out-only patient bootstrap CI."""
import json, sys
import numpy as np, pandas as pd
from pathlib import Path
from lifelines import CoxPHFitter
from sklearn.model_selection import train_test_split
from run_transport_stability import load_cohort, COHORTS
sys.path.insert(0, 'src')
from recurscan.data.dataset import assemble

ROOT = Path(__file__).resolve().parent
SEED = 20260928
d = assemble(); g = list(d.gene_names)
f = pd.DataFrame({f'x{j}': d.X_expr[:, j].astype(float) for j in range(len(g))})
f['t'] = d.time.astype(float); f['e'] = d.event.astype(int)
beta = CoxPHFitter(penalizer=.05).fit(f, 't', 'e').params_.values

def slope(t, e, eta):
    return float(CoxPHFitter(penalizer=.05).fit(
        pd.DataFrame({'eta': eta, 't': t, 'e': e}), 't', 'e').params_['eta'])

# integrity gate: committed full-cohort slopes reproduce before any split
prior = json.loads((ROOT/'results/transport_calibration.json').read_text())['cohorts']
cohorts = {}
for name in COHORTS:
    c = load_cohort(name, g)
    eta = c['X'] @ beta
    s_full = slope(c['t'], c['e'], eta)
    assert abs(s_full - prior[name]['slope']) < 5e-5, (name, s_full, prior[name]['slope'])
    cohorts[name] = (c, eta)
print('integrity gate passed: committed full-cohort slopes reproduced', flush=True)

rng = np.random.default_rng(SEED)
rows = []
for name in COHORTS:
    c, eta = cohorts[name]
    t, e = c['t'], c['e']
    tr, te = train_test_split(np.arange(len(t)), test_size=.5,
                              random_state=SEED, stratify=e)
    if e[tr].sum() < 5:
        rows.append({'cohort': name, 'skipped': True, 'reason': 'training half <5 events'})
        continue
    alpha = slope(t[tr], e[tr], eta[tr])
    recal = slope(t[te], e[te], alpha * eta[te])
    frozen = slope(t[te], e[te], eta[te])
    dist, skipped = [], 0
    for b in range(1000):
        ix = rng.integers(0, len(te), len(te))
        if len(np.unique(e[te][ix])) < 2:
            skipped += 1; continue
        dist.append(slope(t[te][ix], e[te][ix], (alpha * eta[te])[ix]))
    rows.append({'cohort': name, 'skipped': False, 'n_train': int(len(tr)),
                 'n_heldout': int(len(te)), 'events_train': int(e[tr].sum()),
                 'events_heldout': int(e[te].sum()),
                 'alpha_train': alpha, 'heldout_slope_frozen': frozen,
                 'heldout_slope_recalibrated': recal,
                 'boot95': [float(x) for x in np.percentile(dist, [2.5, 97.5])],
                 'boot_skipped': skipped,
                 'pass_band_0.8_1.25': bool(0.8 <= recal <= 1.25)})
    print(name, 'alpha', round(alpha, 3), 'frozen', round(frozen, 3),
          'recal', round(recal, 3), 'pass', rows[-1]['pass_band_0.8_1.25'], flush=True)

valid = [r for r in rows if not r['skipped']]
out = {'protocol': 'docs/PREREG_RECALIBRATION_PIVOT_20260928.md',
       'seed': SEED, 'split': 'stratified 50/50 by event, one fixed split per cohort',
       'band': [0.8, 1.25], 'n_bootstrap': 1000, 'rows': rows,
       'n_pass': int(sum(r['pass_band_0.8_1.25'] for r in valid)),
       'n_valid': len(valid),
       'item_verdict_positive': bool(sum(r['pass_band_0.8_1.25'] for r in valid) >= 4 and len(valid) == 5),
       'limits': 'One fixed split per cohort; held-out halves ~23-55 events so CIs are wide and are reported, not gated; endpoint heterogeneity inherited from the committed record; pass claims only out-of-sample stability of local one-parameter recalibration, not clinical utility or discrimination transport.'}
(ROOT/'results/recalibration_pivot.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps({k: v for k, v in out.items() if k != 'rows'}, indent=1), flush=True)
