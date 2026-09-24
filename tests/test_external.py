"""Hermetic tests for recurscan.external (HTTP mocked)."""
from unittest import mock
from recurscan import external


def test_gprofiler_parses(monkeypatch):
    monkeypatch.setattr(external, "_req",
        lambda url, data=None, timeout=25: {"result": [
            {"source": "KEGG", "name": "Cell cycle", "p_value": 1e-8}]})
    out = external.gprofiler_enrich(["CEP55", "CDC20"])
    assert out[0]["term_name"] == "Cell cycle" and out[0]["p_value"] == 1e-8


def test_ensembl_lookup(monkeypatch):
    monkeypatch.setattr(external, "_get",
        lambda url, timeout=25: {"id": "ENSG00000138107", "biotype": "protein_coding"})
    assert external.ensembl_lookup("CEP55")["id"].startswith("ENSG")


def test_string_enrichment(monkeypatch):
    monkeypatch.setattr(external, "_get",
        lambda url, timeout=25: [{"term": "GO:0000278", "p_value": 1e-6}])
    assert external.string_enrichment(["CDC20", "BIRC5"])[0]["term"].startswith("GO:")


def test_europepmc(monkeypatch):
    monkeypatch.setattr(external, "_get",
        lambda url, timeout=25: {"resultList": {"result": [{"id": "PM1"}]}})
    assert external.europepmc_search("CEP55 breast cancer")[0]["id"] == "PM1"


def test_cache(tmp_path, monkeypatch):
    monkeypatch.setattr(external, "CACHE", str(tmp_path))
    n = {"i": 0}
    def f():
        n["i"] += 1
        return [1]
    external._cached("a", f); external._cached("a", f)
    assert n["i"] == 1


def test_cli_risk_scoring(tmp_path, capsys):
    import json
    from recurscan.cli import main
    model = {"features": ["NPI", "GRADE"], "coef": {"x0": 0.5, "x1": 0.1},
             "means": [3.0, 2.0], "stds": [1.0, 1.0], "tercile_cut_log": [-0.3, 0.3],
             "heldout_cindex": 0.68, "trained_on": "synthetic"}
    p = tmp_path / "m.json"
    p.write_text(json.dumps(model))
    main(["risk", "--model", str(p), "--features", '{"NPI": 5.0}'])
    out = json.loads(capsys.readouterr().out)
    # (5-3)/1*0.5 + (2-2)/1*0.1 = 1.0 -> high band
    assert out["log_hazard"] == 1.0 and out["band"] == "high"


def test_kegg_exact_match(monkeypatch):
    import recurscan.external as ex
    class R:
        def __init__(self, b): self.b = b
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self): return self.b
    def fake(req, timeout=20):
        url = req.full_url
        if "find/genes" in url:
            return R(b"hsa:136647\tMPLKIP, C7orf11; PLK1-interacting\nhsa:5347\tPLK1, PLK; kinase PLK1\n")
        return R(b"hsa:5347\tpath:hsa04110\n")
    monkeypatch.setattr("urllib.request.urlopen", fake)
    out = ex.kegg_pathways("PLK1")
    assert out["kegg_id"] == "hsa:5347"


def test_hgnc_and_trials(monkeypatch):
    import recurscan.external as ex
    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self):
            return b'{"response": {"docs": [{"symbol": "CEP55", "name": "centrosomal protein 55"}]}}'
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=20: R())
    assert ex.hgnc_symbol("CEP55")["symbol"] == "CEP55"
    monkeypatch.setattr(ex, "_get", lambda url, timeout=25: {"studies": [{"nctId": "NCT1"}]})
    assert len(ex.clinicaltrials_search("CEP55")) == 1


def test_reactome_analyze(monkeypatch):
    import recurscan.external as ex
    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self):
            return b'{"pathways": [{"stId": "R-HSA-1", "name": "Mitosis", "entities": {"fdr": 1e-9}}]}'
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R())
    out = ex.reactome_analyze(["CDC20"])
    assert out[0]["stId"] == "R-HSA-1" and out[0]["fdr"] == 1e-9


def test_opentargets_target(monkeypatch):
    import recurscan.external as ex, json as j
    responses = iter([
        {"data": {"search": {"hits": [{"id": "ENSG1", "name": "CDC20"}]}}},
        {"data": {"target": {"approvedSymbol": "CDC20", "associatedDiseases": {"rows": [{"disease": {"name": "breast carcinoma"}, "score": 0.4}]}}}},
    ])
    class R:
        def __init__(self, payload): self.payload = payload
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self): return j.dumps(self.payload).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R(next(responses)))
    out = ex.opentargets_target("CDC20")
    assert out["ensembl_id"] == "ENSG1" and out["top_diseases"][0]["name"] == "breast carcinoma"


def test_intact_interactions(monkeypatch):
    import recurscan.external as ex
    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self):
            return b"uniprotkb:Q12834\tuniprotkb:Q9NS23-2\tintact:EBI-1\nuniprotkb:Q12834\tuniprotkb:O76009\tintact:EBI-2"
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R())
    out = ex.intact_interactions("CDC20 AND taxidA:9606")
    assert out["n_rows"] == 2 and "Q12834" in out["partners"] and "O76009" in out["partners"]


def test_hpa_gene_prognostics(monkeypatch):
    import recurscan.external as ex
    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self):
            return b'Gene\tCancer prognostics - Breast Invasive Carcinoma (TCGA)\nCDC20\tunprognostic (9.83e-3)'
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R())
    out = ex.hpa_gene_prognostics("ENSG00000117399")
    assert out["gene"] == "CDC20"
    assert any("Breast" in k for k in out["prognostics"])


def test_alphafold_prediction(monkeypatch):
    import recurscan.external as ex, json as j
    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self):
            return j.dumps([{"modelEntityId": "AF-Q12834-F1", "gene": "CDC20",
                             "globalMetricValue": 84.12, "fractionPlddtVeryHigh": 0.615}]).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R())
    out = ex.alphafold_prediction("Q12834")
    assert out["gene"] == "CDC20" and out["mean_plddt"] == 84.12


def test_quickgo_annotations(monkeypatch):
    import recurscan.external as ex, json as j
    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self):
            return j.dumps({"numberOfHits": 2, "results": [{"goId": "GO:0000278"}, {"goId": "GO:0000278"}]}).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R())
    out = ex.quickgo_annotations("Q12834")
    assert out["go_ids"] == ["GO:0000278"] and out["n_hits"] == 2


def test_mobidb_entry(monkeypatch):
    import recurscan.external as ex, json as j
    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self):
            return j.dumps([{"acc": "Q12834", "gene": "CDC20", "length": 499,
                             "prediction-disorder-priority": {"content_fraction": 0.21},
                             "homology-domain-pfam": {"regions_names": ["WD40"]}}]).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R())
    out = ex.mobidb_entry("Q12834")
    assert out["disorder_fraction"] == 0.21 and out["pfam_domains"] == ["WD40"]


def test_rcsb_search(monkeypatch):
    import recurscan.external as ex, json as j
    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self):
            return j.dumps({"total_count": 2, "result_set": [{"identifier": "4GGA"}, {"identifier": "4GGC"}]}).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R())
    out = ex.rcsb_search("CDC20 anaphase-promoting complex")
    assert out["total"] == 2 and out["top_ids"] == ["4GGA", "4GGC"]


def test_interpro_domains(monkeypatch):
    import recurscan.external as ex, json as j
    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self):
            return j.dumps({"count": 1, "results": [{"metadata": {
                "accession": "IPR001680", "name": "WD40 repeat", "type": "repeat"}}]}).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R())
    out = ex.interpro_domains("Q12834")
    assert out["domains"][0]["name"] == "WD40 repeat"


def test_monarch_gene_diseases(monkeypatch):
    import recurscan.external as ex, json as j
    calls = iter([
        {"items": [{"id": "HGNC:1723", "name": "CDC20", "category": "biolink:Gene", "xref": ["OMIM:603618"]}]},
        {"items": [{"object_label": "oocyte maturation defect 14"}]},
    ])
    class R:
        def __init__(self, payload): self.payload = payload
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self): return j.dumps(self.payload).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R(next(calls)))
    out = ex.monarch_gene_diseases("CDC20")
    assert out["monarch_id"] == "HGNC:1723" and out["causal_diseases"] == ["oocyte maturation defect 14"]


def test_ncbi_gene(monkeypatch):
    import recurscan.external as ex, json as j
    calls = iter([
        {"esearchresult": {"idlist": ["991"]}},
        {"result": {"991": {"name": "CDC20", "description": "cell division cycle 20", "chromosome": "1"}}},
    ])
    class R:
        def __init__(self, payload): self.payload = payload
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self): return j.dumps(self.payload).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R(next(calls)))
    out = ex.ncbi_gene("CDC20")
    assert out["gene_id"] == "991" and out["chromosome"] == "1"


def test_ols_search(monkeypatch):
    import recurscan.external as ex, json as j
    class R:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self):
            return j.dumps({"response": {"docs": [{"obo_id": "MONDO:0004989", "label": "breast carcinoma"}]}}).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R())
    out = ex.ols_search("breast carcinoma")
    assert out["terms"][0]["id"] == "MONDO:0004989"


def test_gtex_median_expression_parses(monkeypatch):
    from recurscan.external import gtex_median_expression
    import io, json as J

    class R(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    payloads = iter([
        J.dumps({"data": [{"gencodeId": "ENSG00000117399.13",
                           "geneSymbolUpper": "CDC20"}]}).encode(),
        J.dumps({"data": [{"median": 1.635, "unit": "TPM",
                           "tissueSiteDetailId": "Breast_Mammary_Tissue"}]}).encode(),
    ])
    monkeypatch.setattr("urllib.request.urlopen",
                        lambda req, timeout=30: R(next(payloads)))
    out = gtex_median_expression("CDC20")
    assert out["gencode_id"] == "ENSG00000117399.13"
    assert out["median_tpm"] == 1.635
    assert out["unit"] == "TPM"


def test_gwas_catalog_snps_parses(monkeypatch):
    from recurscan.external import gwas_catalog_snps
    import io, json as J

    class R(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    payload = J.dumps({"_embedded": {"singleNucleotidePolymorphisms": [
        {"rsId": "rs1", "functionalClass": "missense_variant"},
        {"rsId": "rs2", "functionalClass": "intron_variant"},
        {"rsId": "rs3", "functionalClass": "missense_variant"}]}}).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R(payload))
    out = gwas_catalog_snps("CDC20")
    assert out["n_snps"] == 3
    assert out["functional_classes"] == {"missense_variant": 2, "intron_variant": 1}
    assert out["sample_rsids"] == ["rs1", "rs2", "rs3"]


def test_ucsc_canonical_parses(monkeypatch):
    from recurscan.external import ucsc_canonical_transcript
    import io, json as J

    class R(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    payload = J.dumps({"knownCanonical": [
        {"transcript": "ENST00000310955.11", "protein": "ENSG00000117399",
         "chromStart": 43358980, "chromEnd": 43363203}]}).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R(payload))
    out = ucsc_canonical_transcript("chr1", 43358954, 43363203)
    assert out["transcripts"][0]["transcript"] == "ENST00000310955.11"


def test_gtex_gene_coords_parses(monkeypatch):
    from recurscan.external import gtex_gene_coords
    import io, json as J

    class R(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    payload = J.dumps({"data": [{"chromosome": "chr1", "start": 43358955,
                                 "end": 43363203, "gencodeId": "ENSG00000117399.13",
                                 "geneSymbolUpper": "CDC20"}]}).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R(payload))
    out = gtex_gene_coords("CDC20")
    assert out["gencode_id"] == "ENSG00000117399.13"
    assert out["chrom"] == "chr1"


def test_mygene_query_parses(monkeypatch):
    from recurscan.external import mygene_query
    import io, json as J

    class R(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    payload = J.dumps({"total": 1, "hits": [
        {"symbol": "CDC20", "entrezgene": "991",
         "ensembl": {"gene": "ENSG00000117399"}, "alias": ["p55CDC"]}]}).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R(payload))
    out = mygene_query("CDC20")
    assert out["hits"][0]["ensembl"] == "ENSG00000117399"
    assert out["hits"][0]["n_alias"] == 1


def test_openalex_works_parses(monkeypatch):
    from recurscan.external import openalex_works
    import io, json as J

    class R(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    payload = J.dumps({"meta": {"count": 2008}, "results": [
        {"id": "https://openalex.org/W1", "display_name": "CDC20 in breast cancer",
         "publication_year": 2020, "cited_by_count": 42}]}).encode()
    monkeypatch.setattr("urllib.request.urlopen", lambda req, timeout=30: R(payload))
    out = openalex_works("CDC20 breast cancer")
    assert out["count"] == 2008
    assert out["top"][0]["cited_by"] == 42
