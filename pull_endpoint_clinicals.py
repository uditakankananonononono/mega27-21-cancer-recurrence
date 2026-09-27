"""Pull raw per-patient clinical data for the endpoint-harmonization study
(docs/PREREG_ENDPOINT_HARMONIZATION_20260928.md). Writes one committed raw
clinical CSV per cohort under data_cache/external/geo/ plus a manifest with
sha256 + source. Harmonized endpoints are derived later, in the analysis
script, from these committed raw files. Resumable: skips cohorts whose CSV
already exists. No GEOparse dependency: series-matrix headers are streamed
from NCBI over HTTPS and parsing stops before the expression table."""
import gzip, hashlib, json, re, urllib.request
import pandas as pd

GEO = "data_cache/external/geo"
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
COHORTS = ["GSE2990", "GSE7390", "GSE11121", "GSE20685", "GSE25066"]
SUPPL = {"GSE2990": f"{BASE}/GSE2nnn/GSE2990/suppl/GSE2990_suppl_info.txt"}

def prefix(gse):
    return re.sub(r"\d{3}$", "nnn", gse)

def matrix_files(gse):
    idx = urllib.request.urlopen(f"{BASE}/{prefix(gse)}/{gse}/matrix/", timeout=60).read().decode()
    return sorted(set(re.findall(r'(GSE\d+[^"]*series_matrix\.txt\.gz)', idx)))

def stream_clinical(url):
    """Parse series-matrix header lines; stop at the expression table."""
    resp = urllib.request.urlopen(url, timeout=120)
    fh = gzip.open(resp, "rt", encoding="utf-8", errors="replace")
    acc, chars = [], {}
    for line in fh:
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!Sample_geo_accession"):
            acc = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            k0 = vals[0].partition(":")[0].strip() if vals and ":" in vals[0] else None
            key = k0 or f"char_{len(chars)}"
            if key in chars:  # repeated key on multi-platform files
                key = f"{key}__{len(chars)}"
            chars[key] = vals
    resp.close()
    rows = []
    for i, gsm in enumerate(acc):
        rec = {"geo_accn": gsm}
        for k, vals in chars.items():
            rec[k] = vals[i] if i < len(vals) else ""
        rows.append(rec)
    return rows

man_path = f"{GEO}/clinical_pull_manifest.json"
try:
    man = json.load(open(man_path))
except Exception:
    man = {"files": {}}

for gse in COHORTS:
    out_csv = f"{GEO}/clinical_{gse}.csv"
    try:
        pd.read_csv(out_csv)
        print(gse, "already pulled, skipping", flush=True)
        continue
    except Exception:
        pass
    if gse in SUPPL:
        urllib.request.urlretrieve(SUPPL[gse], f"{GEO}/{gse}_suppl_info.txt")
        df = pd.read_csv(f"{GEO}/{gse}_suppl_info.txt", sep="\t")
        src = SUPPL[gse]
    else:
        rows = []
        files = matrix_files(gse)
        for mf in files:
            rows.extend(stream_clinical(f"{BASE}/{prefix(gse)}/{gse}/matrix/{mf}"))
        df = pd.DataFrame(rows)
        src = f"{BASE}/{prefix(gse)}/{gse}/matrix/ {files} (series-matrix Sample_characteristics_ch1, verbatim)"
    df.to_csv(out_csv, index=False)
    h = hashlib.sha256(open(out_csv, "rb").read()).hexdigest()
    man["files"][gse] = {"csv": out_csv, "sha256": h, "source": src, "n_rows": int(len(df))}
    json.dump(man, open(man_path, "w"), indent=1)
    print(gse, "rows:", len(df), "cols:", len(df.columns), flush=True)
print("manifest updated:", man_path)
