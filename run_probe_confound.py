"""Judge round 2 pre-outcome: technical-probe confounding in stability map.
For each gene, derive probe multiplicity on GPL96 and GPL570 (available
without outcomes), and within-cohort expression SD. After the stability-map
result, test whether technical variables alone predict stability labels;
leave-one-gene-out logistic AUC and Spearman associations. With 70 genes,
these are descriptive sensitivity analyses, not causal tests. Do not claim
biological mechanism if technical signals explain much of the map.
"""
import json, re
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from sklearn.model_selection import LeaveOneOut
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from run_transport_stability import GEO,COHORTS,load_cohort
from recurscan.data.dataset import assemble

def multiplicity(platform,genes):
    cnt={g:0 for g in genes};active=False;header=None
    with open(GEO/(platform+'.txt'),errors='replace') as f:
        for line in f:
            if line.startswith('!platform_table_begin'):active=True;continue
            if line.startswith('!platform_table_end'):break
            if not active:continue
            row=line.rstrip('\n').split('\t')
            if header is None:header=row;continue
            if len(row)<=header.index('Gene Symbol'):continue
            # Count only the first annotated symbol, matching the transport parser.
            g=row[header.index('Gene Symbol')].split(' ///')[0].strip()
            if g in cnt:cnt[g]+=1
    return np.array([cnt[g] for g in genes])

def main():
    d=assemble();genes=list(d.gene_names)
    j=json.load(open('results/transport_stability.json'))
    co={n:load_cohort(n,genes) for n in COHORTS}
    q96=multiplicity('GPL96',genes);q570=multiplicity('GPL570',genes)
    sd=np.array([np.std(co[n]['X'],axis=0) for n in COHORTS])
    # Z-normalized expression SD is exactly ~1 whenever mapped: it is not a
    # valid intensity/noise measure. Report availability and probe count only.
    mapped=np.array([co[n]['mapped'] for n in COHORTS])
    tech=np.column_stack([np.log1p(q96),np.log1p(q570),mapped.sum(0)])
    labels=np.array([sum(g in f['stable_genes'] for f in j['folds']) for g in genes])
    binary=(labels>=3).astype(int)
    out={'n_genes':len(genes),'technical_predictors':['log1p GPL96 probe count','log1p GPL570 probe count','mapped cohort count'],
         'warning':'cohort z-scoring destroys raw intensity/variance; this audit cannot separate every batch effect or biology',
         'genes':{g:{'GPL96_probe_count':int(q96[i]),'GPL570_probe_count':int(q570[i]),
                     'mapped_cohorts':int(mapped[:,i].sum()),'stable_folds':int(labels[i])}
                  for i,g in enumerate(genes)}}
    out['correlations']={k:{'rho':float(spearmanr(tech[:,i],labels).statistic),
                            'p':float(spearmanr(tech[:,i],labels).pvalue)}
                         for i,k in enumerate(out['technical_predictors'])}
    if len(set(binary))>1:
        pred=np.zeros(len(genes))
        for tr,te in LeaveOneOut().split(tech):
            if len(set(binary[tr]))<2:pred[te]=binary[tr].mean();continue
            model=LogisticRegression(max_iter=300).fit(tech[tr],binary[tr])
            pred[te]=model.predict_proba(tech[te])[:,1]
        out['technical_LOO_AUROC']=float(roc_auc_score(binary,pred))
    else:out['technical_LOO_AUROC']=None
    with open('results/probe_confound.json','w') as f:json.dump(out,f,indent=2)
    print('technical AUC',out['technical_LOO_AUROC'],flush=True)
if __name__=='__main__':main()
