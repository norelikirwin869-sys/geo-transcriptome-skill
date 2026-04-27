# GEO Transcriptome Full-Stack Analysis Skill

## Purpose

Use this skill when the user wants to perform, design, explain, document, or reproduce an end-to-end transcriptome analysis based on GEO or similar public expression datasets.

Typical use cases include:

- Disease mechanism exploration
- Differential-expression analysis
- GEO-based biomarker discovery
- Immune-infiltration analysis
- WGCNA module discovery
- Machine-learning diagnostic-gene screening
- Drug repositioning from expression signatures
- Manuscript methods and results drafting

The skill is disease-agnostic. Examples include diabetic sarcopenia, diabetes, cancer, cardiovascular disease, neurodegeneration, inflammatory disease, aging, metabolic disease, and tissue injury.

## Required inputs

Ask for or infer the following:

1. Disease or phenotype
2. Species, preferably human or mouse
3. Tissue, cell type, or model
4. Data type: microarray, RNA-seq counts, single-cell, or mixed
5. Desired analysis depth
6. Whether the user wants code, methods text, figures, or a full project plan
7. Whether external validation or wet-lab validation should be included

If the user provides only a disease name, generate a complete default workflow and clearly mark assumptions.

## Recommended workflow

### 1. Define the scientific question

Create a concise research aim. For example:

> To identify key genes, biological pathways, immune characteristics, and diagnostic biomarkers associated with diabetic sarcopenia using integrated GEO transcriptome analysis.

### 2. GEO search strategy

Generate English search terms combining:

- Disease terms
- Tissue or cell terms
- Species terms
- Platform terms
- Phenotype-related terms

Example for diabetic sarcopenia:

```text
diabetes skeletal muscle expression profiling
type 2 diabetes skeletal muscle RNA-seq
sarcopenia skeletal muscle transcriptome
muscle atrophy diabetes GEO
aging skeletal muscle expression profiling
```

### 3. Dataset screening

Use clear inclusion criteria:

- Disease and control groups are available
- Sample grouping is explicit
- Raw counts or normalized expression matrix is available
- Sample size is sufficient for exploratory analysis
- Platform annotation is available
- Tissue or model is biologically relevant

Use clear exclusion criteria:

- Missing phenotype labels
- Very small sample size
- Unclear tissue source
- Poor platform annotation
- Duplicated or non-independent samples
- Severe confounding that cannot be handled

### 4. Data preprocessing

For microarray:

- Download GEO matrix or CEL files
- Map probes to gene symbols
- Remove probes without gene symbols
- Collapse multiple probes per gene
- Perform log2 transformation if needed
- Normalize expression matrix

For RNA-seq:

- Use raw count matrix when available
- Filter low-expression genes
- Normalize with DESeq2, edgeR, or limma-voom
- Convert identifiers where needed

### 5. Batch-effect correction

If combining datasets:

- Use ComBat or removeBatchEffect
- Evaluate with PCA, UMAP, boxplot, and sample clustering
- Avoid merging platforms when biological or technical differences are too large

### 6. Differential-expression analysis

Microarray:

- Use limma
- Typical threshold: adjusted P value < 0.05 and absolute log2FC > 0.5 or 1

RNA-seq:

- Use DESeq2, edgeR, or limma-voom
- Report log2FC, P value, adjusted P value, and regulation direction

### 7. Functional enrichment

Perform:

- GO Biological Process
- GO Cellular Component
- GO Molecular Function
- KEGG
- Reactome
- Hallmark gene sets
- GSEA based on ranked genes

### 8. WGCNA

Use WGCNA when sample size is adequate. Identify disease-associated modules and intersect module genes with DEGs.

### 9. PPI and hub genes

Use STRING and Cytoscape concepts:

- Build PPI network
- Use MCODE for functional clusters
- Use cytoHubba algorithms such as MCC, Degree, MNC, EPC, Closeness and Betweenness
- Prefer consensus hub genes across multiple methods

### 10. Machine-learning biomarker discovery

Recommended algorithms:

- LASSO logistic regression
- Random forest
- SVM-RFE
- XGBoost, when sample size supports it

Use intersection or consensus ranking to select biomarkers. Validate with ROC curves and independent GEO datasets.

### 11. Immune infiltration

For bulk tissue data, consider:

- ssGSEA
- xCell
- MCP-counter
- CIBERSORT, when input requirements are met

Analyze:

- Group differences in immune-cell scores
- Correlation between hub genes and immune cells
- Immune-related pathway activation

### 12. Mechanistic extensions

Depending on disease biology, evaluate:

- Mitochondrial dysfunction
- Oxidative phosphorylation
- Fatty-acid metabolism
- Glycolysis
- Insulin signaling
- AMPK/mTOR/FoxO signaling
- Autophagy
- Apoptosis
- Ferroptosis
- Pyroptosis
- Senescence

### 13. Regulatory network

Build candidate networks:

- TF-mRNA
- miRNA-mRNA
- lncRNA-miRNA-mRNA ceRNA

Mention databases such as TRRUST, JASPAR, ChEA3, TargetScan, miRDB, miRWalk, starBase, ENCORI and LncBase.

### 14. Drug prediction

Use disease up- and down-regulated genes for drug repositioning.

Possible resources:

- CMap
- L1000FWD
- DGIdb
- DSigDB
- DrugBank
- CTD
- HERB
- SymMap
- TCMSP

### 15. Validation

Recommend at least one of:

- Independent GEO validation cohort
- qRT-PCR
- Western blot
- Immunofluorescence
- Immunohistochemistry
- Animal or cell model validation
- Protein expression validation using HPA when relevant

## Output style

When responding to the user, provide:

1. A project title
2. Scientific hypothesis
3. Dataset search terms
4. Complete workflow
5. Recommended tools and R packages
6. Code templates when requested
7. Expected figures and tables
8. Manuscript-ready methods text if requested
9. Caveats and validation recommendations

## Quality-control rules

- Do not overclaim causality from observational public data.
- Distinguish discovery, validation, and mechanistic confirmation.
- State assumptions about disease group, control group, species and tissue.
- Recommend external validation whenever biomarkers are proposed.
- For clinical conclusions, emphasize that experimental and clinical validation are required.
