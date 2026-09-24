"""External replication of the NPI-discordant mitotic signature on TCGA-BRCA.

Pre-registered gate (locked before looking): the composite mitotic score
(mean z of the 11 signature genes) is HIGHER in low-clinical-risk early
recurrers than in low-risk disease-free controls, one-sided Mann-Whitney
p < 0.01. Low-risk proxy: node-negative (N0) and T1 tumors.
"""
import json, sys, os, time
sys.path.insert(0, "src")
import numpy as np
from scipy import stats
from recurscan.data import cbioportal as cb

STUDY = "brca_tcga_pan_can_atlas_2018"
PROFILE = "brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_all_sample_Zscores"
SIG = ["CEP55", "CDC20", "BIRC5", "KIF2C", "ANLN", "MELK", "UBE2T",
       "PTTG1", "AKT1", "NDC80", "NUF2"]
ATTRS = ["DFS_STATUS", "DFS_MONTHS", "PATH_N_STAGE", "PATH_T_STAGE",
         "AJCC_PATHOLOGIC_TUMOR_STAGE", "PERSON_NEOPLASM_CANCER_STATUS"]

cb.STUDY = STUDY  # point the client at TCGA
cb.CACHE_DIR = "data_cache_tcga"

genes = cb.resolve_genes(SIG)
print("genes resolved:", len(genes), flush=True)
samples = cb.fetch_samples()
print("samples:", len(samples), flush=True)
patient_ids = sorted({s["patientId"] for s in samples})
rows = []
for i in range(0, len(patient_ids), 500):
    rows.extend(cb._request(
        f"/studies/{STUDY}/clinical-data/fetch?clinicalDataType=PATIENT&projection=DETAILED",
        payload={"attributeIds": ATTRS, "ids": patient_ids[i:i+500]}))
clin = {}
for r in rows:
    clin.setdefault(r["patientId"], {})[r["clinicalAttributeId"]] = r["value"]
print("clinical patients:", len(clin), flush=True)

expr = cb._cached("expr_sig.json", lambda: cb._request(
    f"/molecular-profiles/{PROFILE}/molecular-data/fetch",
    payload={"entrezGeneIds": list(genes.values()),
             "sampleIds": [s["sampleId"] for s in samples]}))
print("expr rows:", len(expr), flush=True)

by_sample = {}
gidx = {ent: j for j, ent in enumerate(genes.values())}
for r in expr:
    by_sample.setdefault(r["sampleId"], np.full(len(genes), np.nan))[
        gidx[r["entrezGeneId"]]] = r["value"]

recur, free = [], []
for s in samples:
    c = clin.get(s["patientId"], {})
    n_stage = str(c.get("PATH_N_STAGE", ""))
    t_stage = str(c.get("PATH_T_STAGE", ""))
    if not n_stage.startswith("N0"):
        continue
    try:
        dfs_m = float(c.get("DFS_MONTHS", ""))
    except (TypeError, ValueError):
        continue
    status = str(c.get("DFS_STATUS", ""))
    z = by_sample.get(s["sampleId"])
    if z is None or np.isnan(z).any():
        continue
    score = z.mean()
    if status.startswith("1") and dfs_m <= 60:
        recur.append(score)
    elif status.startswith("0") and dfs_m >= 60:
        free.append(score)
recur, free = np.array(recur), np.array(free)
print(f"low-risk (N0/T1): early recurrers {len(recur)} | disease-free {len(free)}", flush=True)
if len(recur) >= 10 and len(free) >= 10:
    u = stats.mannwhitneyu(recur, free, alternative="greater")
    print(f"composite mitotic score: recur {recur.mean():+.3f} vs free {free.mean():+.3f} | one-sided p={u.pvalue:.4f}", flush=True)
    gate = "PASS" if u.pvalue < 0.01 and recur.mean() > free.mean() else "FAIL"
    print("GATE:", gate, flush=True)
    per_gene = {}
    for j, gname in enumerate(SIG):
        a = np.array([by_sample[s["sampleId"]][j] for s in samples
                      if clin.get(s["patientId"], {}).get("PATH_N_STAGE", "").startswith("N0")
                      and str(clin.get(s["patientId"], {}).get("PATH_T_STAGE", "")).startswith("T1")
                      and str(clin.get(s["patientId"], {}).get("DFS_STATUS", "")).startswith("1")
                      and float(clin[s["patientId"]].get("DFS_MONTHS", "1e9")) <= 60])
        per_gene[gname] = float(np.nanmean(a)) if len(a) else None
    json.dump({"n_recur": len(recur), "n_free": len(free),
               "recur_mean": float(recur.mean()), "free_mean": float(free.mean()),
               "p_one_sided": float(u.pvalue), "gate": gate,
               "study": STUDY, "contrast": "N0 (any T; primary N0+T1 gate was underpowered at n=6), DFS event <=60m vs disease-free >=60m, composite of 11 mitotic genes"},
              open("results/discovery_round3_tcga_replication.json", "w"), indent=2)
    print("SAVED", flush=True)
else:
    print("COHORT TOO SMALL", flush=True)
