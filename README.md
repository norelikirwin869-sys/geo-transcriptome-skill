# GEO Transcriptome Full-Stack Analysis Skill

A reusable ChatGPT Skill for end-to-end transcriptome analysis based on GEO datasets. It is designed for disease-oriented bioinformatics projects such as diabetic sarcopenia, diabetes-related skeletal muscle dysfunction, cancer, neurodegenerative disease, inflammation, and other public-transcriptome mining scenarios.

## What this skill does

This skill helps an assistant or analyst build a reproducible GEO transcriptomics workflow covering:

- Disease and control definition
- GEO dataset search strategy
- Dataset inclusion and exclusion criteria
- Microarray and RNA-seq preprocessing
- Batch-effect correction
- Differential expression analysis
- GO, KEGG, Reactome, Hallmark and GSEA analysis
- WGCNA co-expression network analysis
- PPI network and hub-gene discovery
- Machine-learning biomarker screening
- ROC validation
- Immune infiltration analysis
- Metabolism, mitochondrial and cell-death pathway analysis
- miRNA, TF, lncRNA and ceRNA regulatory-network construction
- Drug prediction and molecular docking planning
- External validation and manuscript-ready method writing

## Recommended repository structure

```text
geo-transcriptome-skill/
├── SKILL.md
├── README.md
├── LICENSE
├── docs/
│   └── workflow.md
├── examples/
│   └── diabetic_sarcopenia_prompt.md
├── scripts/
│   ├── 01_geo_download_limma.R
│   ├── 02_deseq2_rnaseq.R
│   ├── 03_enrichment_gsea.R
│   ├── 04_wgcna_template.R
│   └── 05_ml_biomarker_template.R
└── templates/
    ├── project_config.yaml
    ├── dataset_screening_table.csv
    └── manuscript_methods_template.md
```

## Quick start

1. Choose a disease or phenotype.
2. Search GEO with disease, tissue, species and platform keywords.
3. Fill in `templates/project_config.yaml`.
4. Use the R templates in `scripts/` to run analysis.
5. Summarize methods using `templates/manuscript_methods_template.md`.

## Example research topic

> Integrated bioinformatics and machine-learning analysis identifies immune- and mitochondrial-related biomarkers in diabetic sarcopenia based on GEO transcriptome datasets.

## Disclaimer

This repository provides research workflow templates only. Results from public transcriptome mining require independent validation through external datasets and wet-lab experiments before biological or clinical interpretation.
