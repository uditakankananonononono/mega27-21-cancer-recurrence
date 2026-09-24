"""cBioPortal public API client for METABRIC (no auth, rate-limited politely).

Network access lives ONLY in this module; tests use canned payloads.
"""
from __future__ import annotations

import json
import os
import time
import urllib.request

API = "https://www.cbioportal.org/api"
STUDY = "brca_metabric"
EXPR_PROFILE = "brca_metabric_mrna_median_all_sample_Zscores"
CNA_PROFILE = "brca_metabric_cna"
MUT_PROFILE = "brca_metabric_mutations"

CACHE_DIR = os.environ.get(
    "RECURSCAN_CACHE",
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "data_cache"),
)

# PAM50 + core breast-cancer prognosis / driver genes (HUGO symbols).
GENE_PANEL = [
    "ACTR3B", "ANLN", "BAG1", "BCL2", "BIRC5", "BLVRA", "CCNB1", "CCNE1",
    "CDC20", "CDC6", "CDH3", "CENPF", "CEP55", "CXXC5", "EGFR", "ERBB2",
    "ESR1", "EXO1", "FGFR4", "FOXA1", "FOXC1", "GPR160", "GRB7", "GSTM1",
    "KIF2C", "KRT14", "KRT17", "KRT5", "MAPT", "MDM2", "MELK", "MIA",
    "MKI67", "MLPH", "MMP11", "MYBL2", "MYC", "NAT1", "NDC80", "NUF2",
    "ORC6", "PGR", "PHGDH", "PTTG1", "RRM2", "SFRP1", "SLC39A6", "TMEM45B",
    "TYMS", "UBE2C", "UBE2T",
    "TP53", "PIK3CA", "GATA3", "AKT1", "CDH1", "PTEN", "BRCA1", "BRCA2",
    "PALB2", "ATM", "CHEK2", "RB1", "CCND1", "MCL1", "BCL2L1", "VEGFA",
    "HIF1A", "STAT3", "MTOR",
]

PATIENT_ATTRS = [
    "RFS_MONTHS", "RFS_STATUS", "OS_MONTHS", "OS_STATUS",
    "AGE_AT_DIAGNOSIS", "LYMPH_NODES_EXAMINED_POSITIVE", "NPI", "CELLULARITY",
    "CHEMOTHERAPY", "HORMONE_THERAPY", "RADIO_THERAPY", "CLAUDIN_SUBTYPE",
    "INTCLUST", "THREEGENE", "ER_IHC", "VITAL_STATUS", "HISTOLOGICAL_SUBTYPE",
]
SAMPLE_ATTRS = [
    "GRADE", "TUMOR_SIZE", "ER_STATUS", "PR_STATUS", "HER2_STATUS",
    "TUMOR_STAGE", "SAMPLE_TYPE",
]


class ApiError(RuntimeError):
    pass


def _request(path: str, payload=None, timeout: int = 90):
    url = f"{API}{path}"
    data = None
    headers = {"Accept": "application/json"}
    if payload is not None:
        data = json.dumps(payload).encode()
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=headers)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode())
        except Exception as exc:
            if attempt == 2:
                raise ApiError(f"GET/POST {path} failed: {exc}") from exc
            time.sleep(2 ** attempt)


def _cached(name: str, fn):
    os.makedirs(CACHE_DIR, exist_ok=True)
    path = os.path.join(CACHE_DIR, name)
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    out = fn()
    with open(path, "w") as fh:
        json.dump(out, fh)
    return out


def resolve_genes(symbols=None) -> dict:
    """HUGO symbol -> entrez id, resolved against the live API."""
    symbols = symbols or GENE_PANEL

    def _fetch():
        genes = _request("/genes/fetch?geneIdType=HUGO_GENE_SYMBOL&projection=SUMMARY",
                         payload=symbols)
        return {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in genes}
    return _cached("genes.json", _fetch)


def fetch_samples():
    def _fetch():
        out, page = [], 0
        while True:
            batch = _request(f"/studies/{STUDY}/samples?projection=SUMMARY"
                             f"&pageSize=2000&pageNumber={page}")
            out.extend(batch)
            if len(batch) < 2000:
                return out
            page += 1
    return _cached("samples.json", _fetch)


def fetch_clinical(samples):
    """Merged per-sample clinical rows (patient attrs joined via patientId)."""
    def _fetch():
        sample_ids = [s["sampleId"] for s in samples]
        patient_ids = sorted({s["patientId"] for s in samples})
        by_patient = {}
        for i in range(0, len(patient_ids), 500):
            chunk = patient_ids[i:i + 500]
            for row in _request(
                    f"/studies/{STUDY}/clinical-data/fetch?clinicalDataType=PATIENT&projection=DETAILED",
                    payload={"attributeIds": PATIENT_ATTRS, "ids": chunk}):
                by_patient.setdefault(row["patientId"], {})[
                    row["clinicalAttributeId"]] = row["value"]
        merged = {}
        for i in range(0, len(sample_ids), 500):
            chunk = sample_ids[i:i + 500]
            for row in _request(
                    f"/studies/{STUDY}/clinical-data/fetch?clinicalDataType=SAMPLE&projection=DETAILED",
                    payload={"attributeIds": SAMPLE_ATTRS, "ids": chunk}):
                merged.setdefault(row["sampleId"], {})[
                    row["clinicalAttributeId"]] = row["value"]
        out = {}
        for s_ in samples:
            rec = dict(by_patient.get(s_["patientId"], {}))
            rec.update(merged.get(s_["sampleId"], {}))
            out[s_["sampleId"]] = rec
        return out
    return _cached("clinical.json", _fetch)


def fetch_expression(entrez_ids, sample_ids):
    """Long-format [{entrezGeneId, sampleId, value}] for the panel."""
    def _fetch():
        rows = []
        for i in range(0, len(sample_ids), 500):
            chunk = sample_ids[i:i + 500]
            rows.extend(_request(
                f"/molecular-profiles/{EXPR_PROFILE}/molecular-data/fetch",
                payload={"entrezGeneIds": entrez_ids, "sampleIds": chunk}))
        return rows
    return _cached("expression.json", _fetch)


def fetch_cna(entrez_ids, sample_ids):
    """Discrete CNA calls (-2..2) for the panel genes."""
    def _fetch():
        rows = []
        for i in range(0, len(sample_ids), 500):
            chunk = sample_ids[i:i + 500]
            rows.extend(_request(
                f"/molecular-profiles/{CNA_PROFILE}/molecular-data/fetch",
                payload={"entrezGeneIds": entrez_ids, "sampleIds": chunk}))
        return rows
    return _cached("cna.json", _fetch)


def fetch_mutations(entrez_ids, sample_ids):
    """Non-synonymous mutation presence per gene/sample."""
    def _fetch():
        rows = []
        for i in range(0, len(sample_ids), 500):
            chunk = sample_ids[i:i + 500]
            rows.extend(_request(
                f"/molecular-profiles/{MUT_PROFILE}/mutations/fetch",
                payload={"entrezGeneIds": entrez_ids, "sampleIds": chunk}))
        return [{"entrezGeneId": r["entrezGeneId"], "sampleId": r["sampleId"],
                 "mutationType": r.get("mutationType", "")} for r in rows]
    return _cached("mutations.json", _fetch)
