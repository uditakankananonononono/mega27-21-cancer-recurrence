"""Pre-registered nested-size DeepSurv vs penalized Cox sensitivity.
See docs/PREREG_N_THRESHOLD_20260928.md, committed before outcomes.
"""
import json, os, sys
from pathlib import Path
import numpy as np, pandas as pd, torch
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
sys.path.insert(0,'src')
from recurscan.data.dataset import assemble
from recurscan.benchmark import stratified_event_split, standardize
from recurscan.models.deepsurv import DeepSurv
from recurscan.train import set_seed, train_cox_fullbatch
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/n_threshold_checkpoint.json'
SIZES=[200,400,800,1200,1481]

def nested_train(tr,event,seed):
    rng=np.random.default_rng(seed)
    order=[]
    by={k:rng.permutation(tr[event[tr]==k]).tolist() for k in (0,1)}
    for n in range(1,len(tr)+1):
        target=int(round(n*sum(event[tr]==1)/len(tr)))
        if len(by[1]) and (sum(event[i] for i in order)<target or not by[0]):
            order.append(by[1].pop())
        else:order.append(by[0].pop())
    return np.array(order)

def run(ds,seed,n):
    tr,te=stratified_event_split(ds.event,.25,seed)
    assert len(tr)==1481 and len(te)==494
    sub=nested_train(tr,ds.event,seed)[:n]
    assert not set(sub)&set(te)
    # Match original DeepSurv input order: expression features, then clinical.
    Xe_tr,Xe_te=standardize(ds.X_expr[sub],ds.X_expr[te])
    Xc_tr,Xc_te=standardize(ds.X_clin[sub],ds.X_clin[te])
    Xtr=np.hstack([Xe_tr,Xc_tr]);Xte=np.hstack([Xe_te,Xc_te])
    # Cox is order-invariant, but use the same features and 0.1 penalizer.
    df=pd.DataFrame(Xtr,columns=[f'x{j}' for j in range(Xtr.shape[1])]);df['t']=ds.time[sub];df['e']=ds.event[sub]
    cox=CoxPHFitter(penalizer=.1).fit(df,'t','e')
    coxscore=cox.predict_log_partial_hazard(pd.DataFrame(Xte,columns=df.columns[:-2])).values.ravel()
    cc=float(concordance_index(ds.time[te],-coxscore,ds.event[te]))
    set_seed(seed)
    net=train_cox_fullbatch(DeepSurv(Xtr.shape[1]),Xtr,ds.time[sub],ds.event[sub],epochs=300,lr=1e-3,seed=seed)
    net.eval()
    with torch.no_grad():score=net(torch.tensor(Xte,dtype=torch.float32)).numpy()
    dc=float(concordance_index(ds.time[te],-score,ds.event[te]))
    return {'seed':seed,'n_train':n,'train_events':int(ds.event[sub].sum()),'test_n':len(te),'test_events':int(ds.event[te].sum()),'cox_cindex':cc,'deepsurv_cindex':dc,'delta_deep_minus_cox':dc-cc}

def main():
    ds=assemble();assert len(ds.time)==1975
    rows=json.loads(OUT.read_text())['rows'] if OUT.exists() else []
    done={(x['seed'],x['n_train']) for x in rows}
    for seed in range(5):
      for n in SIZES:
        if (seed,n) in done:continue
        row=run(ds,seed,n);rows.append(row)
        payload={'protocol':'docs/PREREG_N_THRESHOLD_20260928.md','partial':True,'rows':rows}
        tmp=OUT.with_suffix('.tmp');tmp.write_text(json.dumps(payload,indent=1)+'\n');os.replace(tmp,OUT)
        print(json.dumps(row),flush=True)
    summary={}
    for n in SIZES:
        d=[x['delta_deep_minus_cox'] for x in rows if x['n_train']==n]
        summary[str(n)]={'n_seeds':len(d),'mean_delta':float(np.mean(d)),'min_delta':float(np.min(d)),'max_delta':float(np.max(d)),
                        'descriptive_advantage':bool(np.mean(d)>0)}
    (ROOT/'results/n_threshold.json').write_text(json.dumps({'protocol':'docs/PREREG_N_THRESHOLD_20260928.md','rows':rows,'summary':summary},indent=1)+'\n')
    print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
