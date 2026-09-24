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
