"""Live external pulls: enrichment of the R3 mitotic signature + lit check."""
import json, sys
sys.path.insert(0, "src")
from recurscan import external

SIG = ["CEP55","CDC20","BIRC5","KIF2C","ANLN","MELK","UBE2T","PTTG1","AKT1","NDC80","NUF2"]
out = {}
out["gprofiler"] = external._cached("gp_r3sig", lambda: external.gprofiler_enrich(SIG))
out["string_enrich"] = external._cached("se_r3sig", lambda: external.string_enrichment(SIG))
out["ensembl"] = {g: external._cached("ens_" + g, lambda g=g: external.ensembl_lookup(g)) for g in SIG}
out["epmc"] = external._cached("epmc_r3", lambda: external.europepmc_search(
    '(CEP55 OR CDC20 OR KIF2C) AND "breast cancer" AND (recurrence OR relapse)'))
json.dump(out, open("results/external_pull.json", "w"), default=str)
print("gprofiler terms:", len(out["gprofiler"]), "| string terms:", len(out["string_enrich"]))
print("top:", [(t["source"], t["term_name"], t["p_value"]) for t in out["gprofiler"][:4]])
print("ensembl ok:", sum(1 for v in out["ensembl"].values() if v.get("id")), "/", len(SIG))
