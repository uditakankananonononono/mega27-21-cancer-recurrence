"""Queue item 6: strictly all-four sign stable panel, prereg 2026-09-28."""
import json,sys
import numpy as np,pandas as pd
from pathlib import Path
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
sys.path.insert(0,'src')
from recurscan.data.dataset import assemble
from run_transport_stability import COHORTS,load_cohort

ROOT=Path(__file__).parent
old=json.loads((ROOT/'results/transport_stability.json').read_text())
ds=assemble();genes=list(ds.gene_names); assert set(old['gene_map'])==set(genes)
frame=pd.DataFrame({f'x{j}':ds.X_expr[:,j].astype(float) for j in range(len(genes))})
frame['t']=ds.time.astype(float);frame['e']=ds.event.astype(int)
beta=CoxPHFitter(penalizer=.05).fit(frame,'t','e').params_.values
assert all(abs(beta[j]-old['gene_map'][g]['metabric_coef'])<1e-8 for j,g in enumerate(genes)), 'original METABRIC coefficients differ'
cohorts={h:load_cohort(h,genes) for h in COHORTS}
ci=lambda c,s:float(concordance_index(c['t'],-s,c['e']))
full={h:ci(c,c['X']@beta) for h,c in cohorts.items()}
for row in old['folds']:
 h=row['name'];assert abs(full[h]-row['full_cindex'])<1e-8,(h,full[h],row['full_cindex'])
print('original full comparators integrity gate passed',flush=True)

folds=[]
for h,c in cohorts.items():
 mask=[]
 for j,g in enumerate(genes):
  signs=[old['gene_map'][g]['cohorts'][other]['univariate_coef'] for other in COHORTS if other!=h]
  mask.append(bool(beta[j]!=0 and np.isfinite(beta[j]) and all(x is not None and np.isfinite(x) and x!=0 and np.sign(x)==np.sign(beta[j]) for x in signs)))
 mask=np.asarray(mask,bool)
 d={'name':h,'n_selected':int(mask.sum()),'selected_genes':[genes[j] for j in np.flatnonzero(mask)],'full_cindex':full[h]}
 if d['n_selected']>=3:
  d['reduced_cindex']=ci(c,c['X']@(beta*mask));d['delta']=d['reduced_cindex']-full[h]
 else:d['reduced_cindex']=None;d['delta']=None
 folds.append(d)
 print(h,'selected',d['n_selected'],'delta',d['delta'],flush=True)

out={'protocol':'docs/PREREG_STABILITY_RERUN_STRICT_20260928.md','source_gene_map':'results/transport_stability.json',
     'rule':'all FOUR training cohorts must each agree with METABRIC coefficient sign; no missing or zero coefficient',
     'folds':folds,'n_bootstrap_requested':2000,'bootstrap_seed':1,
     'limits':'Same five cohorts and endpoints, retrospective masked METABRIC coefficients. No clinical, causal or independent validation claim.'}
if all(f['n_selected']>=3 for f in folds):
 out['mean_delta']=float(np.mean([f['delta'] for f in folds]))
 rng=np.random.default_rng(1);dist=[];skipped=0
 for b in range(2000):
  vals=[]
  for f in folds:
   c=cohorts[f['name']]; ix=rng.integers(0,len(c['t']),len(c['t']))
   if len(np.unique(c['e'][ix]))<2:break
   cc={'t':c['t'][ix],'e':c['e'][ix]};X=c['X'][ix]
   mask=np.isin(genes,f['selected_genes'])
   try:vals.append(ci(cc,X@(beta*mask))-ci(cc,X@beta))
   except (ZeroDivisionError,ValueError):break
  if len(vals)==5:dist.append(float(np.mean(vals)))
  else:skipped+=1
 assert dist,'all bootstrap draws degenerate'
 lo,hi=np.percentile(dist,[2.5,97.5]);out['paired_boot_ci95']=[float(lo),float(hi)]
 out['n_bootstrap_valid']=len(dist);out['n_bootstrap_skipped']=skipped
 out['two_sided_statistical_gate']=bool(lo>0 or hi<0)
 out['interpretation']='reduced panel loses concordance' if hi<0 else 'positive direction, descriptive only' if lo>0 else 'difference sign uncertain'
else:
 out['two_sided_statistical_gate']=False;out['interpretation']='at least one fold has fewer than three genes, delta undefined'
(ROOT/'results/transport_stability_rerun_strict.json').write_text(json.dumps(out,indent=1)+'\n')
print(json.dumps(out,indent=1),flush=True)
