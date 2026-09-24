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
