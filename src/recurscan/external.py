"""External research-tool clients for recurrence-signature verification.
Live calls only outside CI (hermetic tests mock _get/_post)."""
import json, os, time, urllib.parse, urllib.request

CACHE = os.path.join(os.path.dirname(__file__), "..", "..", "data_cache", "external")


def _req(url, data=None, timeout=25):
    body = json.dumps(data).encode() if data is not None else None
    req = urllib.request.Request(url, data=body,
        headers={"User-Agent": "recurscan/0.1", "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode())


def _get(url, timeout=25):
    return _req(url, None, timeout)


def _cached(name, fetcher):
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, name + ".json")
    if os.path.exists(path):
        return json.load(open(path))
    d = fetcher()
    json.dump(d, open(path, "w"))
    time.sleep(0.2)
    return d


def gprofiler_enrich(genes, sources=("KEGG", "REAC", "GO:BP")):
    """g:Profiler API: pathway enrichment for a gene list."""
    d = _req("https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
             {"organism": "hsapiens", "query": genes, "sources": list(sources)})
    return [{"source": r["source"], "term_name": r["name"], "p_value": r["p_value"]}
            for r in d.get("result", [])]


def ensembl_lookup(symbol):
    """Ensembl REST: gene id + biotype for a human symbol."""
    url = (f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}"
           "?content-type=application/json")
    return _get(url)


def string_enrichment(genes, species=9606):
    """STRING API: functional enrichment among genes."""
    ids = "%0d".join(genes)
    url = (f"https://string-db.org/api/json/enrichment?identifiers={ids}"
           f"&species={species}")
    return _get(url)


def europepmc_search(query, page_size=5):
    """Europe PMC REST: literature corroboration."""
    q = urllib.parse.quote(query)
    url = (f"https://www.ebi.ac.uk/europepmc/webservices/rest/search"
           f"?query={q}&format=json&pageSize={page_size}")
    return _get(url).get("resultList", {}).get("result", [])


def geo_series_matrix(gse, destdir=None):
    """GEOparse: download + parse a GEO series matrix (external validation cohorts)."""
    import GEOparse
    g = GEOparse.get_GEO(geo=gse, destdir=destdir or os.path.join(CACHE, "geo"),
                         annotate_gpl=False, silent=True)
    return g


def kegg_pathways(symbol):
    """KEGG REST: pathways for a human gene symbol (exact match)."""
    import urllib.request as u
    def _txt(url):
        with u.urlopen(u.Request(url, headers={"User-Agent": "recurscan/0.1"}), timeout=20) as r:
            return r.read().decode()
    hits = [l for l in _txt(f"https://rest.kegg.jp/find/genes/{symbol}").strip().split("\n")
            if l.split("\t")[0].startswith("hsa:")]
    exact = [l for l in hits if l.split("\t")[1].split(";")[0].split(",")[0].strip() == symbol]
    pick = (exact or hits or [None])[0]
    if pick is None:
        return {"kegg_id": None, "pathways": []}
    kegg_id = pick.split("\t")[0]
    txt = _txt(f"https://rest.kegg.jp/link/pathway/{kegg_id}")
    return {"kegg_id": kegg_id,
            "pathways": [l.split("\t")[1] for l in txt.strip().split("\n") if "\t" in l]}


def hgnc_symbol(symbol):
    """HGNC REST: approved symbol/name."""
    import urllib.request as u
    url = f"https://rest.genenames.org/fetch/symbol/{symbol}"
    req = u.Request(url, headers={"Accept": "application/json", "User-Agent": "recurscan/0.1"})
    with u.urlopen(req, timeout=20) as r:
        docs = json.loads(r.read().decode()).get("response", {}).get("docs", [])
    return docs[0] if docs else {}


def clinicaltrials_search(query, page_size=5):
    """ClinicalTrials.gov API v2."""
    import urllib.parse
    q = urllib.parse.quote(query)
    url = (f"https://clinicaltrials.gov/api/v2/studies?query.term={q}"
           f"&pageSize={page_size}&fields=NCTId,BriefTitle,OverallStatus")
    return _get(url).get("studies", [])


def ucsc_xena_datasets(host="tcga.xenahubs.net", hub_dataset="TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"):
    """UCSC Xena: field list for a hosted dataset (proves hub access)."""
    url = f"https://{host}/data/"
    import urllib.request as u
    req = u.Request(url + hub_dataset, headers={"User-Agent": "recurscan/0.1"})
    with u.urlopen(req, timeout=25) as r:
        return r.read().decode()[:2000]


def reactome_analyze(genes, page_size=10):
    """Reactome Analysis Service: over-representation for a gene list."""
    import urllib.request as u
    url = (f"https://reactome.org/AnalysisService/identifiers/"
           f"?pageSize={page_size}&page=1&sortBy=ENTITIES_FDR&order=ASC")
    req = u.Request(url, data="\n".join(genes).encode(),
                    headers={"Content-Type": "text/plain", "User-Agent": "recurscan/0.1"})
    with u.urlopen(req, timeout=30) as r:
        d = json.loads(r.read().decode())
    return [{"stId": p["stId"], "name": p["name"],
             "fdr": p["entities"]["fdr"]} for p in d.get("pathways", [])]


def opentargets_target(symbol, size=5):
    """Open Targets GraphQL: top disease associations for a gene symbol."""
    import urllib.request as u
    q = {"query": '{ search(queryString: "%s", entityNames: ["target"]) { hits { id name } } }' % symbol}
    req = u.Request("https://api.platform.opentargets.org/api/v4/graphql",
                    data=json.dumps(q).encode(),
                    headers={"Content-Type": "application/json", "User-Agent": "recurscan/0.1"})
    with u.urlopen(req, timeout=30) as r:
        hits = json.loads(r.read().decode())["data"]["search"]["hits"]
    ensembl_id = hits[0]["id"]
    q2 = {"query": '{ target(ensemblId: "%s") { approvedSymbol associatedDiseases(page: {index: 0, size: %d}) { rows { disease { name } score } } } }' % (ensembl_id, size)}
    req2 = u.Request("https://api.platform.opentargets.org/api/v4/graphql",
                     data=json.dumps(q2).encode(),
                     headers={"Content-Type": "application/json", "User-Agent": "recurscan/0.1"})
    with u.urlopen(req2, timeout=30) as r:
        d = json.loads(r.read().decode())["data"]["target"]
    return {"ensembl_id": ensembl_id, "symbol": d["approvedSymbol"],
            "top_diseases": [{"name": r2["disease"]["name"], "score": r2["score"]}
                             for r2 in d["associatedDiseases"]["rows"]]}


def intact_interactions(query, max_lines=1000):
    """IntAct via PSICQUIC REST: binary interaction partners (tab25)."""
    import urllib.request as u
    from urllib.parse import quote
    url = ("https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices"
           f"/current/search/query/{quote(query)}?format=tab25")
    with u.urlopen(u.Request(url, headers={"User-Agent": "recurscan/0.1"}), timeout=30) as r:
        lines = r.read().decode().strip().split("\n")[:max_lines]
    pairs = set()
    for ln in lines:
        f = ln.split("\t")
        if len(f) >= 2:
            pairs.add((f[0].replace("uniprotkb:", ""), f[1].replace("uniprotkb:", "")))
    return {"query": query, "n_rows": len(lines), "unique_pairs": len(pairs),
            "partners": sorted({b for a, b in pairs} | {a for a, b in pairs})}


def hpa_gene_prognostics(ensembl_id):
    """Human Protein Atlas per-gene TSV: TCGA cancer prognostics columns."""
    import urllib.request as u, csv, io
    url = f"https://www.proteinatlas.org/{ensembl_id}.tsv"
    with u.urlopen(u.Request(url, headers={"User-Agent": "recurscan/0.1"}), timeout=30) as r:
        text = r.read().decode()
    rows = list(csv.reader(io.StringIO(text), delimiter="\t"))
    hdr, rec = rows[0], rows[1]
    return {"gene": rec[hdr.index("Gene")],
            "prognostics": {h: rec[i] for i, h in enumerate(hdr)
                            if "prognostics" in h.lower()}}


def alphafold_prediction(uniprot_acc):
    """AlphaFold DB API: predicted-structure metadata for a UniProt accession."""
    import urllib.request as u
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_acc}"
    with u.urlopen(u.Request(url, headers={"User-Agent": "recurscan/0.1"}), timeout=30) as r:
        d = json.loads(r.read().decode())
    if not d:
        raise KeyError(uniprot_acc)
    m = d[0]
    return {"accession": uniprot_acc, "model": m["modelEntityId"],
            "gene": m.get("gene"), "mean_plddt": m["globalMetricValue"],
            "frac_very_high": m.get("fractionPlddtVeryHigh")}


def quickgo_annotations(uniprot_acc, aspect="biological_process", limit=50):
    """QuickGO (EBI GOA): GO annotations for a UniProt accession."""
    import urllib.request as u
    url = (f"https://www.ebi.ac.uk/QuickGO/services/annotation/search"
           f"?geneProductId=UniProtKB:{uniprot_acc}&goAspect={aspect}&limit={limit}")
    with u.urlopen(u.Request(url, headers={"User-Agent": "recurscan/0.1"}), timeout=30) as r:
        d = json.loads(r.read().decode())
    return {"accession": uniprot_acc, "n_hits": d["numberOfHits"],
            "go_ids": sorted({x["goId"] for x in d["results"]})}


def mobidb_entry(uniprot_acc):
    """MobiDB API: disorder/evidence summary for a UniProt accession."""
    import urllib.request as u
    url = f"https://mobidb.org/api/download?acc={uniprot_acc}"
    with u.urlopen(u.Request(url, headers={"User-Agent": "recurscan/0.1"}), timeout=30) as r:
        d = json.loads(r.read().decode())
    if not d:
        raise KeyError(uniprot_acc)
    m = d[0]
    dis = m.get("prediction-disorder-priority", {})
    pfam = m.get("homology-domain-pfam", {})
    return {"accession": uniprot_acc, "gene": m.get("gene"), "length": m.get("length"),
            "disorder_fraction": dis.get("content_fraction"),
            "pfam_domains": pfam.get("regions_names", [])}


def rcsb_search(text, rows=10):
    """RCSB search API v2: full-text structure search."""
    import urllib.request as u
    from urllib.parse import quote
    q = json.dumps({"query": {"type": "terminal", "service": "full_text",
                              "parameters": {"value": text}},
                    "return_type": "entry",
                    "request_options": {"paginate": {"start": 0, "rows": rows}}})
    url = f"https://search.rcsb.org/rcsbsearch/v2/query?json={quote(q)}"
    with u.urlopen(u.Request(url, headers={"User-Agent": "recurscan/0.1"}), timeout=30) as r:
        d = json.loads(r.read().decode())
    return {"query": text, "total": d["total_count"],
            "top_ids": [x["identifier"] for x in d.get("result_set", [])]}


def interpro_domains(uniprot_acc, page_size=25):
    """InterPro API: integrated domain/family entries for a UniProt protein."""
    import urllib.request as u
    url = (f"https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/"
           f"{uniprot_acc}/?page_size={page_size}")
    with u.urlopen(u.Request(url, headers={"User-Agent": "recurscan/0.1"}), timeout=30) as r:
        d = json.loads(r.read().decode())
    return {"accession": uniprot_acc, "count": d["count"],
            "domains": [{"accession": x["metadata"]["accession"],
                         "name": x["metadata"]["name"],
                         "type": x["metadata"]["type"]} for x in d["results"]]}


def monarch_gene_diseases(symbol):
    """Monarch Initiative v3 API: search + causal disease associations."""
    import urllib.request as u
    from urllib.parse import quote
    def _get(url):
        with u.urlopen(u.Request(url, headers={"User-Agent": "recurscan/0.1"}), timeout=30) as r:
            return json.loads(r.read().decode())
    s = _get(f"https://api.monarchinitiative.org/v3/api/search?q={quote(symbol)}&limit=3")
    hit = next((i for i in s["items"] if i["name"] == symbol and i["category"] == "biolink:Gene"), s["items"][0])
    a = _get("https://api.monarchinitiative.org/v3/api/association"
             f"?category=biolink:CausalGeneToDiseaseAssociation&entity={hit['id']}&limit=10")
    return {"symbol": symbol, "monarch_id": hit["id"], "xrefs": hit.get("xref", []),
            "causal_diseases": [x.get("object_label") or x.get("object") for x in a["items"]]}


def ncbi_gene(symbol, organism="human"):
    """NCBI E-utilities: GeneID + summary for a gene symbol."""
    import urllib.request as u
    from urllib.parse import quote
    def _get(url):
        with u.urlopen(u.Request(url, headers={"User-Agent": "recurscan/0.1"}), timeout=30) as r:
            return json.loads(r.read().decode())
    q = quote(f"{symbol}[sym] AND {organism}[orgn]")
    s = _get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={q}&retmode=json")
    ids = s["esearchresult"]["idlist"]
    if not ids:
        raise KeyError(symbol)
    summ = _get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={ids[0]}&retmode=json")
    doc = summ["result"][ids[0]]
    return {"symbol": symbol, "gene_id": ids[0], "name": doc.get("name"),
            "description": doc.get("description"), "chromosome": doc.get("chromosome")}
