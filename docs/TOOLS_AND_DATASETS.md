# External tools and datasets - item 21 (honest registry, updated as used)

Counting rule: accession-level datasets actually used; one study's condition
matrix = one dataset. Tiers: USED IN RESULTS vs STAGED.

## External research tools (web resources + databases)
Used in verified results:
1. cBioPortal API - METABRIC + TCGA-BRCA PanCanAtlas pulls (benchmark, R3 replication attempt)
2. lifelines 0.30.0 - Cox PLL cross-verification (<1e-3), logrank

Staged (client live + tested, cached):
3. g:Profiler gost API - R3 signature enrichment: 100 terms, top = mitotic
   cell cycle process (biological coherence confirmed)
4. STRING API enrichment - 173 functional terms for R3 signature
5. Europe PMC REST - recurrence literature corroboration
6. GEOparse 2.x (installed) - GEO cohort validation route (GSE2034, GSE2990,
   GSE7390, GSE11121, GSE96058 queued for external validation #2+)
7. Ensembl REST - BLOCKED live (timeouts 9:13 PM; client + tests ready, retry queued)

## Packages
8. PyTorch  9. NumPy  10. pandas  11. scikit-learn  12. SciPy
13. networkx  14. matplotlib  15. pytest (11 hermetic tests)

Honest tool count: 15 (7 resources, 8 packages)
Path to 40: MSigDB, COSMIC, DepMap, GDSC, Reactome, KEGG, WikiPathways,
Enrichr, UCSC Xena, GDC API, GEO (per-cohort), Oncomine successors, HGNC,
BioGRID, IntAct, OncoKB, CIViC, DGIdb, PharmGKB, ClinicalTrials.gov,
gnomAD, dbSNP, Ensembl VEP, Protein Atlas (HPA), KM-plotter, cBioPortal
studies beyond BRCA.

## Datasets (accession-level, actually pulled)
Used in results: METABRIC (1), TCGA-BRCA PanCanAtlas (1) = 2
Staged: 100 g:Profiler term records + 173 STRING enrichment records +
5 Europe PMC articles + GEOparse GEO route (0 cohorts pulled yet) = 278 records
but by rule these count as 3 dataset pulls (g:Profiler result set, STRING
result set, EuropePMC result set).
Honest dataset count: 5 (2 used in results, 3 staged)
Path to 120+: GEO breast cohorts (30+ accessions queued), cBioPortal
non-BRCA recurrence cohorts (~10), DepMap/CCLE files, METABRIC data layers
counted once per rule.
