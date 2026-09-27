"""Render evidence tables from committed results, without re-estimating outcomes."""
import json
from pathlib import Path
R=Path(__file__).resolve().parent.parent/'results'
O=Path(__file__).resolve().parent/'audit_tables.tex'
load=lambda name:json.loads((R/name).read_text())
bench=load('metabric_rfs_benchmark.json')
trans=load('transport_stability.json')
cal=load('transport_calibration.json')
probe=load('probe_confound.json')
lines=[]
def add(*s):lines.extend(s)
def esc(s):return str(s).replace('_',r'\_')
add(r'\subsection{Reproducibility audit of the internal benchmark}',
    'The stratified 25\\% holdout contains 494 patients and 200 events in each of the five seeds. '+
    'Table~\\ref{tab:seed-audit} reports every test-set concordance as stored in '+
    r'\texttt{results/metabric\_rfs\_benchmark.json}. Repeated splits of the same cohort are '+
    'not five external replications; their spread measures sensitivity to this split protocol, not population transfer. '+
    'The fold-specific percentile intervals below are the original bootstrap outputs, not confidence bounds on a five-seed average.',
    r'\begin{table}[htbp]\centering\small',
    r'\begin{tabular}{rccccc}\hline',
    r'Seed & Clinical Cox & Clinical+expression Cox & DeepSurv & CoxCNN & CoxGCN \\ \hline')
keys=['coxph_clin','coxph_full','deepsurv','coxcnn','coxgcn']
for i in range(5):
    row=[f'{bench["models"][k]["per_seed"][i]["c_index"]:.4f}' for k in keys]
    add(str(i)+' & '+' & '.join(row)+r' \\')
add(r'\hline\end{tabular}',
    r'\caption{All five held-out seed results for the original METABRIC experiment. These are repeated partitions of one cohort.}\label{tab:seed-audit}',r'\end{table}',
    'The full Cox model exceeds the clinical-only model in all five partitions, '+
    'but the mean difference is small relative to the uncertainty of any individual test set. '+
    'No deep architecture exceeds the full Cox mean. The cohort-level result cannot be used as a '+
    'replacement for independent verification of a new signature.',
    r'\begin{table}[htbp]\centering\small',
    r'\begin{tabular}{rccc}\hline',
    r'Seed & Clinical+expression Cox & 95\% bootstrap CI & Test events \\ \hline')
for z in bench['models']['coxph_full']['per_seed']:
    lo,hi=z['c_index_ci95']
    add(f'{z["seed"]} & {z["c_index"]:.4f} & ({lo:.4f}, {hi:.4f}) & {z["events"]} '+r'\\')
add(r'\hline\end{tabular}',r'\caption{Individual holdout bootstrap intervals for the full Cox benchmark, retained without pooling them into a spurious external-replication interval.}\label{tab:seed-ci}',r'\end{table}')
add(r'\subsection{External cohort accounting and comparator limits}',
    'Table~\\ref{tab:transport-audit} unpacks the five external comparisons. '+
    'The frozen expression coefficients were transported without fitting to the evaluation outcomes, '+
    'but within-cohort expression standardization and probe mapping remain consequential preprocessing. '+
    'Point-estimate ordering is reported without a claim that overlapping intervals establish superiority. '+
    'The cohorts have different event definitions, treatment settings and mapped-gene counts; '+
    'pooling their concordances as if they were one common endpoint would conceal those differences.',
    r'\begin{table}[htbp]\centering\footnotesize',
    r'\begin{tabular}{lrrllp{1.25in}}\hline',
    r'Cohort & Patients & Events & Frozen panel & Comparator & Interpretation \\ \hline')
rows=[('GSE7390',196,89,'0.6009','GGI 0.5366; Veridex-76 0.5757','RFS point-estimate beat'),
      ('GSE25066',508,111,'0.6400','GGI class 0.6025','DRFS; chemotherapy cohort'),
      ('GSE11121',200,46,'0.6979','Grade 0.6294','DMFS; node-negative'),
      ('GSE2990',187,67,'0.6559','GGI continuous 0.6651','RFS; comparator ahead'),
      ('GSE20685',327,83,'0.6264','Nodal stage 0.6996','Metastasis; comparator ahead')]
for name,n,e,score,other,note in rows:
    add(f'{name} & {n} & {e} & {score} & {other} & {note} '+r'\\')
add(r'\hline\end{tabular}',
    r'\caption{Cohort-specific endpoint and comparator audit from the committed transport records. Sample numbers for the calibration audit can differ by a filtering step; the analyses are not interchangeable.}\label{tab:transport-audit}',r'\end{table}',
    '\\par In GSE25066 the DLDA-30 metadata direction failed an internal check against the observed outcome '+
    'and was excluded; the chemotherapy-sensitivity comparator (0.6422) is effectively tied with the panel. '+
    'The GSE2990 GGI score is continuous rather than a two-level published category. '+
    'GSE20685 nodal stage wins by a substantial point-estimate margin, so the expression score is not '+
    'a general substitute for measured clinical risk.')
add(r'\subsection{Absolute-risk transfer audit}',
    'Discrimination does not establish calibration. The 60-month METABRIC baseline survival was '+
    'applied without cohort-specific refitting, and Table~\\ref{tab:calibration-audit} gives the '+
    'stored slope estimates. Only GSE7390 and GSE2990 have a comparable recurrence-free endpoint; '+
    'the other rows are descriptive stress tests rather than absolute-risk validations. '+
    'Cohort-wise standardization further limits interpretation of transported absolute risk.',
    r'\begin{table}[htbp]\centering\small',
    r'\begin{tabular}{lrlccc}\hline',
    r'Cohort & $n$ & RFS comparable? & Slope (95\% CI) & Predicted 5-y event & KM observed \\ \hline')
for name,c in cal['cohorts'].items():
    lo,hi=c['slope_ci95']
    add(f'{name} & {c["n"]} & {"yes" if c["endpoint_RFS_comparable"] else "no"} & '+
        f'{c["slope"]:.3f} ({lo:.3f}, {hi:.3f}) & {c["mean_predicted"]:.3f} & {c["observed_KM_5yr_event"]:.3f} '+r'\\')
add(r'\hline\end{tabular}',
    r'\caption{Descriptive five-year calibration audit from \texttt{results/transport\_calibration.json}; non-RFS endpoints cannot validate RFS absolute probability.}\label{tab:calibration-audit}',r'\end{table}',
    '\\par GSE7390 has a slope of 0.464 (95\\% CI 0.077--0.852); its overall mean predicted '+
    'five-year event probability is 0.264 versus Kaplan--Meier estimate 0.286. '+
    'Mean agreement conceals a shallow slope and non-monotone grouped observations. '+
    'GSE2990 has slope 1.030 with a wide interval (0.408--1.652). '+
    'Neither outcome justifies individual-patient risk prediction.')
add(r'\subsection{Grouped calibration diagnostics for comparable endpoints}',
    'The two recurrence-free cohorts allow descriptive quartile comparisons at 60 months. '+
    'Table~\\ref{tab:calibration-quartiles} retains all grouped predictions and Kaplan--Meier '+
    'estimates; risk-set counts expose censoring at the target horizon. The groups are '+
    'within-cohort quantiles, not common clinical cutoffs. They cannot validate the score '+
    'for individual decisions or be mixed with the three differently defined endpoints.',
    r'\begin{table}[htbp]\centering\small',
    r'\begin{tabular}{llrrrr}\hline',
    r'Cohort & Quantile & $n$ & At risk at 5 y & Mean predicted & KM observed \\ \hline')
for name in ['GSE7390','GSE2990']:
    assert cal['cohorts'][name]['endpoint_RFS_comparable']
    for i,g in enumerate(cal['cohorts'][name]['calibration_groups'],1):
        add(f'{name} & Q{i} & {g["n"]} & {g["at_risk_5yr"]} & '+
            f'{g["mean_predicted_5yr_event"]:.3f} & {g["km_observed_5yr_event"]:.3f} '+r'\\')
add(r'\hline\end{tabular}',
    r'\caption{All recorded calibration quartiles for the two recurrence-free external cohorts. Five-year risk-set counts can be small; compare groups cautiously.}\label{tab:calibration-quartiles}',
    r'\end{table}',
    'In GSE7390, the highest-risk quartile has lower observed recurrence at five years '+
    'than the third quartile even though its predicted risk is higher, consistent with '+
    'the slope warning. In GSE2990 the top quartile separates more clearly, but '+
    'the middle quartiles are close. These observed curves include censoring uncertainty '+
    'not displayed as confidence bands here; this descriptive table cannot substitute '+
    'for prespecified calibration tests on an independently matched treatment cohort.')
add(r'\subsection{Leave-one-cohort-out stability falsification}',
    'The stability experiment trains its panel-selection rule on METABRIC and four external cohorts '+
    'while withholding the fifth external cohort from gene selection. A gene must match the METABRIC '+
    'coefficient direction with a bootstrap interval excluding zero in all four training cohorts. '+
    'The reduced score masks unstable METABRIC coefficients; it is not re-trained on the holdout. '+
    'Table~\\ref{tab:stability-audit} is the complete five-fold outcome, including missing-by-design '+
    'reduced scores where no gene passed the rule.',
    r'\begin{table}[htbp]\centering\small',
    r'\begin{tabular}{lrrrl}\hline',
    r'Held-out cohort & Stable genes & Full $C$ & Reduced $C$ & Selected genes \\ \hline')
for f in trans['folds']:
    add(f'{f["name"]} & {f["n_stable"]} & {f["full_cindex"]:.4f} & '+
        (f'{f["stable_cindex"]:.4f}' if f['stable_cindex'] is not None else 'undefined')+
        ' & '+(', '.join(esc(g) for g in f['stable_genes']) or 'none')+r' \\')
add(r'\hline\end{tabular}',r'\caption{Locked stability gate. A zero-gene fold has no reduced score; this is gate failure, not a zero-valued concordance.}\label{tab:stability-audit}',r'\end{table}',
    '\\par Three zero-gene folds make the predeclared five-fold mean delta and the planned bootstrap and '+
    'random-size-matched gene-label permutation comparisons undefined. In the remaining folds, '+
    'the reduced score also performs worse. Computing a mean over only those two folds would change '+
    'the target of inference after seeing outcomes. The result rejects this stable-subpanel route '+
    'to a biological discovery; it does not erase the separate full-panel comparator audit.',
    r'\subsection{Biology-versus-measurement boundary}',
    'The technical audit counted annotated GPL96 and GPL570 probes for the panel genes and compared '+
    'those counts with fold-stability frequency. Spearman correlations were not significant: '+
    'GPL96 $\\rho=-0.127$, $p=0.296$; GPL570 $\\rho=-0.181$, $p=0.134$; '+
    'mapped-cohort count $\\rho=0.092$, $p=0.447$. Yet no gene was stable in three folds, '+
    'so binary technical leave-one-gene-out AUROC is undefined. '+
    'A non-significant correlation cannot rule out batch or probe bias, especially after '+
    'within-cohort z-scoring has discarded raw intensity and variance information. '+
    'The separately locked Hallmark coherence test likewise cannot run past its zero-gene-fold gate. '+
    'Neither test licenses pathway enrichment or a causal mechanism claim.')
assert not probe['technical_LOO_AUROC'] and not trans['gate_pass']
O.write_text('\n'.join(lines)+'\n')
print(f'wrote {O}, {len(lines)} lines')


# --- verdict item: gene-level readout of the fitted model (bounded) ---
risk=json.loads((Path(__file__).resolve().parent.parent/'results'/'risk_model.json').read_text())
feats=risk['features']; coef=risk['coef']; n_feat=len(feats)
rows=sorted(((coef['x%d'%i],feats[i]) for i in range(n_feat)), key=lambda t:-t[0])
def esc(s): return s.replace('_', '\\_')
lines2=['\\subsection{Fitted-model coefficient audit (bounded gene-level readout)}',
 'The deployed risk model is the penalized Cox fit on 1,580 METABRIC patients with 17 clinical '
 'covariates and the 70-gene expression panel (87 features; held-out $n=395$, $C=0.6891$, '
 '\\texttt{results/risk\\_model.json}). Table~\\ref{tab:coef-audit} lists the twelve largest '
 'risk-increasing and twelve largest risk-decreasing standardized coefficients. Three caveats '
 'bound every reading: the fit is penalized, so correlated features split and shrink weight; '
 'expression inputs are z-scored, so a coefficient compares one-SD expression shifts, not '
 'presence or absence of a gene; and coefficients are associational. A counterintuitive sign '
 'on any single gene is not evidence about mechanism and is not used as one anywhere in this paper.',
 '\\begin{table}[htbp]\\centering\\footnotesize',
 '\\begin{tabular}{lr|lr}\\hline',
 '\\multicolumn{2}{c|}{Risk-increasing} & \\multicolumn{2}{c}{Risk-decreasing} \\\\ ',
 'Feature & $\\hat\\beta$ & Feature & $\\hat\\beta$ \\\\ \\hline']
for i in range(12):
    pv,pn=rows[i]; nv,nn=rows[-(12-i)]
    lines2.append('%s & %+.4f & %s & %+.4f \\\\' % (esc(pn), pv, esc(nn), nv))
lines2+=['\\hline\\end{tabular}',
 '\\caption{Largest standardized coefficients of the committed risk model. Signs and magnitudes are exactly as stored; no feature was re-estimated for this table.}\\label{tab:coef-audit}',
 '\\end{table}',
 'The pattern is consistent with the rest of the paper rather than a new finding: clinical '
 'burden variables (positive lymph-node count, NPI, tumor size) dominate the risk-increasing '
 'side, joined by mitotic-panel genes (UBE2C, CCNB1, MELK), while hormone-receptor-pathway '
 'expression (ESR1, PGR, GATA3, MAPT) and recorded hormone and radiation treatment sit on the '
 'risk-decreasing side. The treatment indicators encode assignment, not benefit, and their '
 'signs reflect confounding by indication; they are reported, not interpreted as effects.']
with open(O,'a') as fh: fh.write('\n'.join(lines2)+'\n')
print('appended coefficient audit,', len(lines2), 'lines')

# ---- five-year Brier + decision-curve audit (results/brier_dca.json)
bd=load('brier_dca.json')
n0=len(lines)
add(r'\subsection{Five-year Brier score and decision-curve audit}',
    'Discrimination and calibration slopes do not say whether acting on the score helps. '+
    'Using the frozen panel and the transported 60-month baseline '+
    r'(\texttt{run\_brier\_dca.py}, \texttt{results/brier\_dca.json}; design locked in the same commit as first run), '+
    'Table~\\ref{tab:brier-audit} reports the IPCW Brier score at 60 months against a null model '+
    'that assigns every patient the cohort Kaplan--Meier five-year event probability. '+
    'Inverse-probability-of-censoring weights follow the standard time-point form; patients censored '+
    'before 60 months contribute through the censoring model only. The METABRIC row is apparent '+
    'performance on the training cohort and is not validation.',
    r'\begin{table}[htbp]\centering\small',
    r'\begin{tabular}{llrrrr}\hline',
    r'Cohort & Role & $n$ & Brier (IPCW) & Null Brier & Skill \\ \hline')
for name,c in bd['cohorts'].items():
    role='apparent' if name.startswith('METABRIC') else 'external'
    add(f"{esc(name)} & {role} & {c['n']} & {c['brier_ipcw_60m']:.4f} & {c['brier_null_60m']:.4f} & {c['brier_skill_vs_null']:.3f} "+r'\\')
add(r'\hline\end{tabular}',
    r'\caption{IPCW Brier scores at 60 months. Skill is $1 - B/B_{\mathrm{null}}$; values near zero mean the score barely beats assigning everyone the cohort average.}\label{tab:brier-audit}',
    r'\end{table}',
    'On the two RFS-comparable external cohorts the skill over the null is small: 0.026 in GSE7390 '+
    'and 0.072 in GSE2990. The score adds a little absolute-risk information beyond the cohort average, '+
    'not enough to support individual treatment decisions. Decision-curve analysis over threshold '+
    'probabilities 0.05--0.60 (same IPCW weights; Table~\\ref{tab:dca-audit}) shows net benefit above '+
    'treat-all mainly in the 0.15--0.25 threshold band; at low thresholds the model and treat-all are '+
    'nearly identical, and at 0.30 and above treat-all goes negative while the model stays positive but small.',
    r'\begin{table}[htbp]\centering\small',
    r'\begin{tabular}{lrrrr}\hline',
    r'Cohort & NB @ 0.15 & NB @ 0.20 & NB @ 0.25 & NB @ 0.30 \\ \hline')
for name,c in bd['cohorts'].items():
    nb={r['pt']:r['net_benefit'] for r in c['decision_curve']}
    add(f"{esc(name)} & {nb[0.15]:.3f} & {nb[0.20]:.3f} & {nb[0.25]:.3f} & {nb[0.30]:.3f} "+r'\\')
add(r'\hline\end{tabular}',
    r'\caption{Decision-curve net benefit at four thresholds (treat-none is zero by construction). Descriptive audit only; thresholds for care are a clinical, not statistical, decision.}\label{tab:dca-audit}',
    r'\end{table}',
    'These curves are descriptive: they use a transported baseline without cohort refit, assume '+
    'independent censoring given time, and cannot recommend a clinical threshold. They bound, rather '+
    'than promote, the score\'s decision value.')
with open(O,'a') as fh: fh.write('\n'.join(lines[n0:])+'\n')
print('appended brier/dca audit,', len(lines)-n0, 'lines')
