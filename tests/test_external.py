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
