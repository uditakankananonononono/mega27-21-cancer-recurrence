"""Pre-outcome round-1 transport-stability map.
Five leave-one-cohort-out evaluations prevent selecting genes on a test cohort.
This is a retrospective research analysis, not a clinical predictor.

Locked before result: gene direction stability = METABRIC coefficient sign
consistent with four training-cohort univariate Cox signs, each with
95% bootstrap interval excluding zero. Missing probe or insufficient events
makes gene unclassified. Stable-only panel is the METABRIC-trained full Cox
coefficients masked to stable genes (no external-cohort refit). Primary metric:
mean of five held-out C-index differences, stable minus full; one-sided
paired bootstrap lower bound >0 and alpha .01; 1000 random-gene-set
permutations matching each fold's stable-panel size, p<.01. If there are
no stable genes, the gate fails rather than redefining the threshold.
Provisional follow-up: no 5-year absolute calibration claim without
cohort-specific baseline survival; separate calibration audit required.
"""
import gzip, json, math, sys
from pathlib import Path
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
sys.path.insert(0,'src')
from recurscan.data.dataset import assemble

ROOT=Path(__file__).resolve().parent
GEO=ROOT/'data_cache/external/geo'
COHORTS=['GSE7390','GSE11121','GSE2990','GSE20685','GSE25066']
ENDPOINT={
 'GSE7390':('t.rfs','e.rfs',1./30.4375), 'GSE11121':('t.dmfs','e.dmfs',1.),
 'GSE2990':('time.rfs','event.rfs',12.),
 'GSE20685':('follow_up_duration (years)','event_metastasis',12.),
 'GSE25066':('drfs_even_time_years','drfs_1_event_0_censored',12.)}

def probes(platform, wanted):
    out={}; active=False; header=None
    with open(GEO/(platform+'.txt'),errors='replace') as f:
        for line in f:
            if line.startswith('!platform_table_begin'):active=True;continue
            if line.startswith('!platform_table_end'):break
            if not active:continue
            fields=line.rstrip('\n').split('\t')
            if header is None:header=fields;continue
            if len(fields)>header.index('Gene Symbol'):
                gene=fields[header.index('Gene Symbol')].split(' ///')[0].strip()
                if gene in wanted:out[fields[0]]=gene
    return out

def load_cohort(name, genes):
    platform='GPL570' if name=='GSE20685' else 'GPL96'
    lookup=probes(platform,set(genes))
    sup=None
    if name=='GSE2990':
        sup=pd.read_csv(GEO/'GSE2990_suppl_info.txt',sep='\t').set_index('geo_accn')
    rows=[]; expression=[]; cur=None; mode=None; header=None
    def flush():
        if cur is None:return
        c=cur['chars']; key_t,key_e,mult=ENDPOINT[name]
        if name=='GSE2990':
            if cur['id'] not in sup.index:return
            s=sup.loc[cur['id']]
            if isinstance(s,pd.DataFrame):s=s.iloc[0]
            val_t,val_e=s[key_t],s[key_e]
        else:val_t,val_e=c.get(key_t),c.get(key_e)
        try:
            t=float(val_t)*mult;e=int(float(val_e))
            if not(np.isfinite(t) and t>0 and e in (0,1)):return
        except (TypeError,ValueError):return
        vec=np.array([cur['genes'].get(g,[0.,0])[0]/cur['genes'][g][1]
                      if g in cur['genes'] and cur['genes'][g][1] else np.nan for g in genes],dtype=float)
        rows.append((t,e));expression.append(vec)
    with gzip.open(GEO/(name+'_family.soft.gz'),'rt',errors='replace') as f:
        for line in f:
            if line.startswith('^SAMPLE'):
                flush();cur={'id':line.split('=',1)[1].strip(),'chars':{},'genes':{}}
                mode=None;header=None
            elif cur is None:continue
            elif line.startswith('!Sample_characteristics_ch1'):
                v=line.split('=',1)[1].strip();k,sep,v=v.partition(':');
                if sep:cur['chars'][k.strip()]=v.strip()
            elif line.startswith('!sample_table_begin'):mode='table'
            elif line.startswith('!sample_table_end'):mode=None
            elif mode=='table':
                parts=line.rstrip('\n').split('\t')
                if header is None:header=parts;continue
                g=lookup.get(parts[0]);
                if not g or len(parts)<2:continue
                try:v=float(parts[1])
                except ValueError:continue
                a=cur['genes'].setdefault(g,[0.,0]);a[0]+=v;a[1]+=1
        flush()
    if not rows:raise RuntimeError(name+' yielded no endpoint rows')
    T,E=np.array(rows).T;X=np.array(expression)
    mapped=np.isfinite(X).any(0)
    # Match the original transport scripts: within-cohort z-score each
    # measured probe mean. Missing genes contribute zero, explicitly flagged.
    for j in range(len(genes)):
        if mapped[j]:
            v=X[:,j];v[~np.isfinite(v)]=np.nanmedian(v)
            X[:,j]=(v-v.mean())/(v.std()+1e-9)
        else:X[:,j]=0
    print(name,len(T),'events',int(E.sum()),'mapped',int(mapped.sum()),flush=True)
    return {'t':T,'e':E.astype(int),'X':X,'mapped':mapped}

def fit_univariate(cohort,j):
    if not cohort['mapped'][j] or cohort['e'].sum()<10:return np.nan
    d=pd.DataFrame({'x':cohort['X'][:,j],'t':cohort['t'],'e':cohort['e']})
    try:return float(CoxPHFitter(penalizer=.05).fit(d,'t','e').params_.iloc[0])
    except Exception:return np.nan

def ci(c,score):return float(concordance_index(c['t'],-score,c['e']))

def main():
    ds=assemble();genes=list(ds.gene_names)
    print('METABRIC',len(ds.time),'events',int(ds.event.sum()),flush=True)
    d={f'x{j}':ds.X_expr[:,j].astype(float) for j in range(len(genes))}
    d['t']=ds.time.astype(float);d['e']=ds.event.astype(int)
    met=CoxPHFitter(penalizer=.05).fit(pd.DataFrame(d),'t','e')
    beta=met.params_.values
    cohorts={name:load_cohort(name,genes) for name in COHORTS}
    coef={name:np.array([fit_univariate(cohorts[name],j) for j in range(len(genes))]) for name in COHORTS}
    # Bootstrap intervals of the per-gene direction estimates, fixed seed.
    rng=np.random.default_rng(0);boot={}
    for name,c in cohorts.items():
        B=np.full((200,len(genes)),np.nan)
        for b in range(200):
            idx=rng.integers(0,len(c['t']),len(c['t']))
            cc={k:(v[idx] if k in ('t','e','X') else v) for k,v in c.items()}
            for j in range(len(genes)):
                if c['mapped'][j]:B[b,j]=fit_univariate(cc,j)
        boot[name]=np.nanpercentile(B,[2.5,97.5],axis=0)
        print('bootstrap',name,flush=True)
    folds=[]
    for name,c in cohorts.items():
        train=[z for z in COHORTS if z!=name]
        stable=np.array([np.all([np.isfinite(coef[z][j]) and
             np.sign(coef[z][j])==np.sign(beta[j]) and
             boot[z][0,j]*boot[z][1,j]>0 and
             np.sign(boot[z][0,j])==np.sign(beta[j]) for z in train])
             for j in range(len(genes))],dtype=bool)
        full=c['X']@beta;reduced=c['X']@(beta*stable)
        folds.append({'name':name,'stable_genes':[genes[j] for j in np.flatnonzero(stable)],
                      'n_stable':int(stable.sum()),'full_cindex':ci(c,full),
                      'stable_cindex':ci(c,reduced) if stable.any() else None})
    valid=all(x['n_stable']>0 for x in folds)
    out={'protocol':'five cohort holdouts; METABRIC Cox coefficients masked; fixed 200 bootstrap intervals',
         'folds':folds,'gene_map':{genes[j]:{'metabric_coef':float(beta[j]),
           'cohorts':{z:{'univariate_coef':float(coef[z][j]) if np.isfinite(coef[z][j]) else None,
                         'boot95':[float(v) if np.isfinite(v) else None for v in boot[z][:,j]]}
                       for z in COHORTS}} for j in range(len(genes))}}
    if valid:
        delta=np.mean([x['stable_cindex']-x['full_cindex'] for x in folds])
        out['mean_delta']=float(delta)
        # Paired sample bootstrap: resample patients within each held-out
        # cohort; fixed masks and frozen coefficients.
        rng=np.random.default_rng(1);dist=[]
        for b in range(2000):
            vals=[]
            for f in folds:
                c=cohorts[f['name']];ix=rng.integers(0,len(c['t']),len(c['t']))
                sub={'t':c['t'][ix],'e':c['e'][ix]};X=c['X'][ix]
                mask=np.isin(genes,f['stable_genes'])
                if sum(sub['e'])<5:break
                vals.append(ci(sub,X@(beta*mask))-ci(sub,X@beta))
            if len(vals)==5:dist.append(np.mean(vals))
        out['paired_boot95']=[float(x) for x in np.percentile(dist,[2.5,97.5])]
        out['paired_boot_one_sided_p']=(1+sum(x<=0 for x in dist))/(1+len(dist))
        rng=np.random.default_rng(2);null=[]
        for b in range(1000):
            vals=[]
            for f in folds:
                c=cohorts[f['name']];sel=rng.choice(len(genes),f['n_stable'],replace=False)
                mask=np.zeros(len(genes),bool);mask[sel]=True
                vals.append(ci(c,c['X']@(beta*mask))-f['full_cindex'])
            null.append(float(np.mean(vals)))
        out['stability_label_permutation_p']=(1+sum(x>=delta for x in null))/(1+len(null))
        out['gate_pass']=bool(delta>0 and out['paired_boot_one_sided_p']<.01
                              and out['stability_label_permutation_p']<.01)
    else:out['gate_pass']=False;out['reason']='at least one held-out fold has zero stable genes'
    with open(ROOT/'results/transport_stability.json','w') as f:json.dump(out,f,indent=2)
    print('verdict',out['gate_pass'],out.get('mean_delta'),flush=True)

if __name__=='__main__':main()
