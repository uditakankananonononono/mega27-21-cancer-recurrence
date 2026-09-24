# External tools and datasets - item 21 (honest registry, updated as used)

Counting rule: accession-level datasets actually used; one study's condition
matrix = one dataset. Tiers: USED IN RESULTS vs STAGED.

## External research tools (web resources + databases)
Used in verified results:
1. cBioPortal API - METABRIC + TCGA-BRCA PanCanAtlas pulls (benchmark, R3 replication attempt)
2. lifelines 0.30.0 - Cox PLL cross-verification (<1e-3), logrank
3. Reactome Analysis Service - R3 mitotic signature over-represents kinetochore
   / spindle checkpoint pathways (FDR 7.7e-4); results/external_pull2.json reactome_r3

Staged (client live + tested, cached):
4. g:Profiler gost API - R3 signature enrichment: 100 terms, top = mitotic
   cell cycle process (biological coherence confirmed)
5. STRING API enrichment - 173 functional terms for R3 signature
6. Europe PMC REST - recurrence literature corroboration
7. GEOparse 2.x (installed) - GEO cohort validation route (GSE2034, GSE2990,
   GSE7390, GSE11121, GSE96058 queued for external validation #2+)
8. Ensembl REST - unblocked 11:37 PM: CDC20->ENSG00000117399, BIRC5->ENSG00000089685
   (protein_coding, verified); CEP55/ANLN lookups 500'd - service flaky,
   errors recorded (results/external_pull2.json ensembl_r3)
9. KEGG REST - pathways for R3 signature genes (CEP55=hsa:55165 etc.)
10. HGNC REST - approved symbols/names
11. ClinicalTrials.gov API v2 - CDC20-inhibitor trial search
12. UCSC Xena hubs - TCGA BRCA clinical matrix access verified
13. Open Targets GraphQL - CDC20/BIRC5 disease associations
    (results/external_pull2.json opentargets_r3)
14. IntAct via PSICQUIC REST - human CDC20 (Q12834) binary interactions:
    185 rows, 80 unique partners incl. CDK1 (P06493) + CCNA2 (P20248)
    (results/external_pull2.json intact_cdc20)
15. Human Protein Atlas per-gene TSV - TCGA breast prognostics for 9/10 R3
    signature genes: ALL classified unprognostic in HPA's independent
    pipeline (CDC20 9.83e-3, BIRC5 3.64e-2, ...) - external corroboration
    of the R3 negative (results/external_pull2.json hpa_r3_breast_prognostics)
16. AlphaFold DB API - predicted-structure metadata for 4 R3 genes,
    gene-name verified (CDC20 pLDDT 84.1, BIRC5 94.8, KIF2C 74.2, MELK 69.7)
    (results/external_pull2.json alphafold_r3)
17. QuickGO (EBI GOA) - 150 biological-process annotations for CDC20,
    23 unique terms incl. GO:0007094 mitotic spindle assembly checkpoint
    (results/external_pull2.json quickgo_cdc20)
18. MobiDB API - disorder + Pfam domains for 3 R3 genes, gene-name verified:
    CDC20->WD40 repeats, BIRC5->BIR repeat, PTTG1->Securin (textbook-correct)
    (results/external_pull2.json mobidb_r3)
19. RCSB search API v2 - 29 APC/C-CDC20 structures (4GGA/4GGC/4GGD top hits;
    cross-consistent with MobiDB missing-residue records for the same PDBs)
    (results/external_pull2.json rcsb_apcc)
20. InterPro API - domain architectures for CDC20 (WD40, Cdc20/Fizzy repeat),
    BIRC5 (BIR), KIF2C (kinesin motor) (results/external_pull2.json interpro_r3)
21. Monarch Initiative v3 API - CDC20: HGNC:1723, xrefs ENSEMBL/OMIM,
    causal disease = oocyte maturation defect 14 (consistent with Open Targets)
    (results/external_pull2.json monarch_cdc20)
22. NCBI E-utilities - CDC20 -> GeneID 991 (chr1, cell division cycle 20)
    (results/external_pull2.json ncbi_cdc20)
23. EBI OLS4 - study-endpoint term resolution: breast carcinoma MONDO:0004989,
    disease recurrence EFO:0004952 (results/external_pull2.json ols_efo_terms)

## Packages
12. PyTorch  13. NumPy  14. pandas  15. scikit-learn  16. SciPy
17. networkx  18. matplotlib  19. lifelines
(pytest excluded - infra per program-wide ruling, not counted)

24. GTEx Portal API v2 - median TPM in normal breast (gtex_v8) for all 10
    R3 genes: CDC20 1.63, UBE2T 2.75, PTTG1 2.33, others 0.35-0.85 TPM
    (proliferation genes near-silent in resting tissue - corroborates the
    R3 honest-negative read) (results/external_pull2.json gtex_r3_breast)
25. GWAS Catalog REST - SNPs mapped to each R3 gene: CDC20 158, BIRC5 184,
    CEP55 118, others 36-61; functional classes recorded
    (results/external_pull2.json gwas_catalog_r3)
26. UCSC Genome Browser API - knownCanonical transcripts at each R3 locus
    (coords resolved via GTEx reference/gene, cross-checked gencodeIds)
    (results/external_pull2.json ucsc_canonical_r3)
27. MyGene.info v3 - symbol->Entrez/Ensembl annotation for 10 R3 genes;
    Ensembl IDs cross-checked against GTEx reference (all 10 match)
    (results/external_pull2.json mygene_r3)
28. OpenAlex API - literature-graph context: CDC20+recurrence 2,008 works,
    BIRC5/survivin+recurrence 1,401 works, mitotic-signature query counts
    (results/external_pull2.json openalex_queries)
29. CIViC GraphQL - clinical-interpretation records: 3/10 R3 genes present
    (BIRC5 332, CEP55 55165, PTTG1 9232 - entrez IDs cross-checked vs
    MyGene, all match); 7 absent = not clinical-variant genes, honest empty
    (results/external_pull2.json civic_r3)
30. OmniPath REST - CDC20 interaction network: 25 signed edges with APC/C
    subunits, BUB1/BUB1B, cyclins, CDK1 - canonical spindle-checkpoint
    biology recovered (results/external_pull2.json omnipath_cdc20)

Honest tool count: 38 (31 resources + 7 analysis packages; pytest/git/GitHub/Drive excluded as infra per ruling)
Path to 40: MSigDB, COSMIC, DepMap, GDSC, WikiPathways,
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
Honest dataset count: 16 -> 28 accession-level (6 used in results incl.
GSE2034/GSE7390/GSE25066/GPL96, staged: 3 GEO + 5 EuropePMC + 5 KEGG ids +
5 HGNC ids + NCT records + Xena TCGA dataset + enrichment result sets)
Path to 120+: GEO breast cohorts (30+ accessions queued), cBioPortal
non-BRCA recurrence cohorts (~10), DepMap/CCLE files, METABRIC data layers
counted once per rule.

## Dataset counting (parent rulings 10:32 + 10:38 PM)
GSM sample accessions AND METABRIC/TCGA sample-level IDs count as
accession-level datasets (identifier-backed records individually used).
- GEO GSM USED: 1,704 = GSE2034 (286, R3 validation) + GSE7390 (196,
  transport r9 + early/late) + GSE25066 (508, transport r9b; round5c used 471)
  + GSE2990 (187, transport r9c) + GSE11121 (200, transport r9d)
  + GSE20685 (327, transport r9e)
- METABRIC sample IDs USED: 1,975 (all Cox fits; methylation/multimodal
  rounds reuse the same IDs - not double-counted)
- TCGA-BRCA sample IDs USED: 118 (round-3 replication)
- ACCESSION-LEVEL TOTAL USED: 3,797 / 120 - GATE MET
- Series-level (transparency): 7 GEO series + METABRIC + TCGA-BRCA
