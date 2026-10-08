"""EXT-058 extension batch (EXTENDS lane 21). Protocol: docs/PREREG_EXT058_20261008.md (committed before this ran)."""
import sys,os,json,math,warnings,numpy as np,pandas as pd
warnings.filterwarnings('ignore')
os.environ.setdefault('RECURSCAN_CACHE','/tmp/mb_cache'); sys.path.insert(0,'src')
from recurscan.data.dataset import assemble
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from lifelines.statistics import proportional_hazard_test
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sksurv.ensemble import RandomSurvivalForest
from sksurv.util import Surv
d=assemble('/tmp/mb_cache'); t=d.time.astype(float); e=d.event.astype(int); n=len(t)
Xc=d.X_clin.astype(float); Xe=d.X_expr.astype(float); cn=d.clin_names; sid=d.sample_ids
clin=json.load(open('/tmp/mb_cache/clinical.json'))
raw5=np.array([[float(clin[s].get(a,'nan')) if str(clin[s].get(a,'')).replace('.','',1).replace('-','',1).isdigit() else np.nan for a in ['AGE_AT_DIAGNOSIS','TUMOR_SIZE','GRADE','LYMPH_NODES_EXAMINED_POSITIVE','NPI']] for s in sid])
OUT='results/ext058.json'; R={'n':n,'events':int(e.sum())}
def save(): json.dump(R,open(OUT,'w'),indent=1,default=float)
def std(tr,te):
    mu=tr.mean(0); sd=tr.std(0); sd[sd<1e-8]=1; return (tr-mu)/sd,(te-mu)/sd
def cox(Xtr,ttr,etr,Xte,pen=0.1,strata=None):
    A,B=std(Xtr,Xte); df=pd.DataFrame(A,columns=[f'f{i}' for i in range(A.shape[1])]); df['T']=ttr; df['E']=etr
    if strata is not None: df['S']=strata[0]
    m=CoxPHFitter(penalizer=pen).fit(df,'T','E',strata=['S'] if strata is not None else None)
    dt=pd.DataFrame(B,columns=df.columns[:A.shape[1]]); 
    if strata is not None: dt['S']=strata[1]
    return np.asarray(m.predict_partial_hazard(dt)).ravel()
def rsf(Xtr,ttr,etr,Xte,seed=0):
    A,B=std(Xtr,Xte); m=RandomSurvivalForest(n_estimators=200,min_samples_leaf=15,n_jobs=2,random_state=seed).fit(A,Surv.from_arrays(etr.astype(bool),ttr)); return m.predict(B)
def oof(X,tt,ee,fn,seed=0,k=5,**kw):
    o=np.zeros(len(tt))
    for tr,te in StratifiedKFold(k,shuffle=True,random_state=seed).split(X,ee): o[te]=fn(X[tr],tt[tr],ee[tr],X[te],**kw)
    return o
C=lambda tt,ee,r:float(concordance_index(tt,-r,ee))
XB=np.hstack([Xc,Xe]); 
o_base=oof(XB,t,e,cox); R['base_oof_C']=round(C(t,e,o_base),4); print('base',R['base_oof_C'],flush=True); save()
# X1
cs=[];rs=[]
for s in range(3):
    cs.append(C(t,e,oof(XB,t,e,cox,seed=s))); rs.append(C(t,e,oof(XB,t,e,rsf,seed=s))); print('X1',s,cs[-1],rs[-1],flush=True)
R['X1']=dict(cox=[round(x,4) for x in cs],rsf=[round(x,4) for x in rs],mean_cox=round(float(np.mean(cs)),4),mean_rsf=round(float(np.mean(rs)),4),gate_rsf_ge_cox_plus_0_01=bool(np.mean(rs)>=np.mean(cs)+0.01)); save()
# X2
age=Xc[:,0]; g2={'<=50':age<=50,'51-65':(age>50)&(age<=65),'>65':age>65}
R['X2']={k:dict(n=int(v.sum()),events=int(e[v].sum()),C=round(C(t[v],e[v],o_base[v]),4)) for k,v in g2.items()}
R['X2_gate_min_stratum_ge_overall_minus_0.05']=bool(min(x['C'] for x in R['X2'].values())>=R['base_oof_C']-0.05); save()
# X3
sub={nm.split('=')[1]:Xc[:,j]==1 for j,nm in enumerate(cn) if nm.startswith('CLAUDIN_SUBTYPE=')}
R['X3']={k:dict(n=int(v.sum()),events=int(e[v].sum()),C=(round(C(t[v],e[v],o_base[v]),4) if e[v].sum()>=2 and v.sum()>5 else None),gated=bool(e[v].sum()>=30)) for k,v in sub.items()}
R['X3_gate']=bool(all(x['C']>=R['base_oof_C']-0.05 for x in R['X3'].values() if x['gated'])); save()
# X4
lm=t>36; o4=oof(XB[lm],t[lm]-36,e[lm],cox); R['X4']=dict(n=int(lm.sum()),events=int(e[lm].sum()),C=round(C(t[lm],e[lm],o4),4)); R['X4_gate_C_ge_0.60']=bool(R['X4']['C']>=0.60); print('X4',R['X4'],flush=True); save()
# X5
r=np.random.default_rng(0); full=cox(XB,t,e,XB); app=C(t,e,full); opt=[]
for b in range(200):
    ix=r.integers(0,n,n); 
    try: Xb,tb,eb=XB[ix],t[ix],e[ix]; pb=cox(Xb,tb,eb,np.vstack([Xb,XB])); opt.append(C(tb,eb,pb[:n])-C(t,e,pb[n:]))
    except Exception as ex: pass
R['X5']=dict(apparent_C=round(app,4),mean_optimism=round(float(np.mean(opt)),4),corrected_C=round(app-float(np.mean(opt)),4),n_boot_ok=len(opt)); R['X5_gate_optimism_le_0.02']=bool(np.mean(opt)<=0.02); print('X5',R['X5'],flush=True); save()
# X6
grp={'age':[0],'size':[1],'grade':[2],'nodes':[3],'NPI':[4],'ER_PR_HER2':[5,6,7],'subtype':list(range(8,14)),'treatment':[14,15,16]}
R['X6']={}
for k,cols in grp.items():
    keep=[j for j in range(XB.shape[1]) if j not in cols]; c=C(t,e,oof(XB[:,keep],t,e,cox)); R['X6'][k]=dict(C_without=round(c,4),delta=round(c-R['base_oof_C'],4))
print('X6',R['X6'],flush=True); save()
# X7
mut=json.load(open('/tmp/mb_cache/mutations.json')); genes=json.load(open('/tmp/mb_cache/genes.json')); e2s={v:k for k,v in genes.items()}; pos={s:i for i,s in enumerate(sid)}
M=np.zeros((n,len(e2s)));gl=sorted(e2s)
for m in mut:
    if m['sampleId'] in pos and m['entrezGeneId'] in gl: M[pos[m['sampleId']],gl.index(m['entrezGeneId'])]=1
keep=np.where(M.mean(0)>=0.03)[0]; R['X7_genes']=[e2s[gl[j]] for j in keep]; XM=np.hstack([XB,M[:,keep]])
dd=[]
for s in range(3): dd.append(C(t,e,oof(XM,t,e,cox,seed=s))-C(t,e,oof(XB,t,e,cox,seed=s)))
R['X7']=dict(n_mut_features=len(keep),delta_C_per_repeat=[round(x,4) for x in dd],mean_delta=round(float(np.mean(dd)),4)); R['X7_gate_mean_delta_ge_0.005']=bool(np.mean(dd)>=0.005); print('X7',R['X7'],flush=True); save()
# X8
ok=(e==1)&(t<=60)|(t>60); y5=((e==1)&(t<=60)).astype(int); idx=np.where(ok)[0]; XS=XB[idx]; ys=y5[idx]; er=Xc[idx,cn.index('ER_STATUS=Positive')]==1; st=np.array([[nm.split('=')[1] for j,nm in enumerate(cn) if nm.startswith('CLAUDIN_SUBTYPE=') and Xc[i,j]==1][:1] or ['NA'] for i in idx]).ravel()
cov=[];size=[];sg={k:[0,0] for k in ['ER+','ER-']+sorted(set(st))}
for s in range(5):
    rg=np.random.default_rng(s); p=rg.permutation(len(idx)); a=len(p)//2; b=a+int(0.3*len(p)); tr,ca,te=p[:a],p[a:b],p[b:]
    A,B=std(XS[tr],XS[np.r_[ca,te]]); lrm=LogisticRegression(C=0.05,max_iter=500).fit(A,ys[tr]); pr=lrm.predict_proba(B); pc,pt=pr[:len(ca)],pr[len(ca):]
    sc=1-pc[np.arange(len(ca)),ys[ca]]; q=np.quantile(sc,min(1,math.ceil((len(ca)+1)*0.9)/len(ca))); sets=(1-pt)<=q
    hit=sets[np.arange(len(te)),ys[te]]; cov.append(hit.mean()); size.append(sets.sum(1).mean())
    for lab,mk in [('ER+',er[te]),('ER-',~er[te])]+[(k,st[te]==k) for k in set(st)]: sg[lab][0]+=int(hit[mk].sum()); sg[lab][1]+=int(mk.sum())
R['X8']=dict(n=len(idx),events5y=int(ys.sum()),marginal_coverage=round(float(np.mean(cov)),4),mean_set_size=round(float(np.mean(size)),3),subgroup_coverage={k:dict(n=v[1],cov=round(v[0]/v[1],3)) for k,v in sg.items() if v[1]>0})
R['X8']['max_subgroup_deviation']=round(max(abs(v['cov']-0.9) for v in R['X8']['subgroup_coverage'].values() if v['n']>=50),3); R['X8_gate_marginal_in_0.87_0.93']=bool(0.87<=np.mean(cov)<=0.93); print('X8',R['X8'],flush=True); save()
# X9
R['X9_missing_fraction']={a:round(float(np.isnan(raw5[:,j]).mean()),4) for j,a in enumerate(['age','size','grade','nodes','NPI'])}
cc=~np.isnan(raw5).any(1); Xcc=XB[cc]; 
def itfn(Xtr,ttr,etr,Xte):  # imputation inside the fold
    return None
oi=np.zeros(n)
for tr,te in StratifiedKFold(5,shuffle=True,random_state=0).split(XB,e):
    Z=np.hstack([raw5,Xc[:,5:]]); imp=IterativeImputer(max_iter=10,random_state=0).fit(Z[tr]); Zi=imp.transform(Z); Xi=np.hstack([Zi,Xe]); oi[te]=cox(Xi[tr],t[tr],e[tr],Xi[te])
occ=oof(Xcc,t[cc],e[cc],cox)
R['X9']=dict(C_median_all=R['base_oof_C'],C_iterative_all=round(C(t,e,oi),4),n_complete_case=int(cc.sum()),C_completecase_model=round(C(t[cc],e[cc],occ),4),C_median_model_on_cc_rows=round(C(t[cc],e[cc],o_base[cc]),4))
R['X9_gate_imputation_change_le_0.01']=bool(abs(R['X9']['C_iterative_all']-R['X9']['C_median_all'])<=0.01); print('X9',R['X9'],flush=True); save()
# X10
Xcl=(Xc-Xc.mean(0))/np.where(Xc.std(0)<1e-8,1,Xc.std(0)); df=pd.DataFrame(Xcl,columns=cn); df['T']=t; df['E']=e
mdl=CoxPHFitter(penalizer=0.1).fit(df,'T','E'); ph=proportional_hazard_test(mdl,df,time_transform='rank').summary
pv=ph['p'].to_dict(); R['X10_ph_p']={k:float(v) for k,v in pv.items()}; R['X10_violators_bonferroni_0.05_over_17']=[k for k,v in pv.items() if v<0.05/17]
worst=min(pv,key=pv.get); wj=cn.index(worst); col=Xc[:,wj]; 
bins=pd.qcut(col,4,labels=False,duplicates='drop') if len(set(col))>2 else col.astype(int)
keep=[j for j in range(XB.shape[1]) if j!=wj]; os_=np.zeros(n)
for tr,te in StratifiedKFold(5,shuffle=True,random_state=0).split(XB,e): os_[te]=cox(XB[tr][:,keep],t[tr],e[tr],XB[te][:,keep],strata=(np.asarray(bins)[tr],np.asarray(bins)[te]))
R['X10']=dict(worst_violator=worst,p=float(pv[worst]),C_stratified_oof=round(C(t,e,os_),4),note='partial hazards from a stratified model are not comparable across strata; C is indicative only',C_base=R['base_oof_C']); print('X10',R['X10'],flush=True); save()
print(json.dumps(R,indent=1,default=float)[:200])
