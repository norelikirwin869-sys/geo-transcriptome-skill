# Full-Stack GEO Transcriptome Workflow

## 1. Research design

Define:

- Disease or phenotype
- Species
- Tissue or cell type
- Case and control groups
- Primary outcome
- Secondary analyses

Example:

> Disease: diabetic sarcopenia  
> Tissue: skeletal muscle  
> Data source: GEO bulk transcriptome datasets  
> Goal: identify diagnostic biomarkers and immune-metabolic mechanisms.

## 2. Dataset discovery

Search GEO using combined keywords:

```text
disease + tissue + expression profiling
disease + tissue + RNA-seq
disease + phenotype + transcriptome
```

For diabetic sarcopenia:

```text
type 2 diabetes skeletal muscle RNA-seq
sarcopenia skeletal muscle transcriptome
diabetic muscle atrophy GEO
aging skeletal muscle expression profiling
```

## 3. Dataset screening table

Create a table with:

- GEO accession
- Platform
- Species
- Tissue
- Case number
- Control number
- Data type
- Group definition
- Inclusion status
- Reason for exclusion if applicable

## 4. Preprocessing

### Microarray

- Load GEO matrix or raw CEL files
- Normalize expression values
- Map probes to gene symbols
- Remove unannotated probes
- Collapse duplicate gene symbols

### RNA-seq

- Use raw counts when available
- Filter low-count genes
- Normalize with DESeq2, edgeR, or limma-voom
- Transform data for visualization

## 5. Quality control

Recommended figures:

- Boxplot before and after normalization
- PCA or UMAP plot
- Sample clustering heatmap
- Density plot

## 6. Differential expression

Recommended methods:

- limma for microarray
- DESeq2 or edgeR for RNA-seq counts
- limma-voom for RNA-seq with linear modeling needs

Common thresholds:

```text
adjusted P value < 0.05
absolute log2FC > 0.5 or 1
```

## 7. Enrichment analysis

Perform:

- GO BP, CC, MF
- KEGG
- Reactome
- Hallmark
- GSEA

Report:

- Term name
- Gene count
- Adjusted P value
- Enrichment score if applicable
- Representative genes

## 8. WGCNA

Use WGCNA when sample size is adequate. Recommended steps:

1. Filter low-variance genes
2. Cluster samples
3. Pick soft-threshold power
4. Build adjacency and TOM matrices
5. Detect modules
6. Correlate modules with phenotype
7. Intersect disease-related module genes with DEGs

## 9. PPI and hub genes

Use STRING to build protein-protein interaction networks and Cytoscape for visualization.

Hub-gene prioritization strategies:

- Degree
- MCC
- MNC
- EPC
- Closeness
- Betweenness
- Intersection with WGCNA module genes
- Validation in independent datasets

## 10. Machine learning

Recommended pipeline:

```text
Candidate genes
→ train/test split or cross-validation
→ LASSO
→ Random forest
→ SVM-RFE
→ consensus biomarkers
→ ROC validation
```

Avoid overfitting when sample size is small.

## 11. Immune infiltration

Bulk transcriptome options:

- ssGSEA
- xCell
- MCP-counter
- CIBERSORT

Analyze:

- Immune-cell differences between groups
- Correlation between hub genes and immune cells
- Association with inflammatory pathways

## 12. Mechanistic modules

Common disease-relevant modules:

- Mitochondrial function
- Oxidative phosphorylation
- Glycolysis
- Fatty-acid metabolism
- Insulin signaling
- AMPK signaling
- PI3K-Akt signaling
- FoxO signaling
- Autophagy
- Apoptosis
- Ferroptosis
- Pyroptosis
- Senescence

## 13. Regulatory network

Build:

- TF-mRNA network
- miRNA-mRNA network
- lncRNA-miRNA-mRNA ceRNA network

## 14. Drug prediction

Use upregulated and downregulated DEGs to query drug-reversal databases such as CMap or L1000FWD. Follow with target validation and molecular docking when appropriate.

## 15. Validation

Recommended validation layers:

- Independent GEO dataset
- qRT-PCR
- Western blot
- Immunohistochemistry
- Immunofluorescence
- Animal model
- Cell model

## 16. Final outputs

Deliver:

- Dataset screening table
- DEG table
- Volcano plot
- Heatmap
- PCA plot
- GO/KEGG/GSEA figures
- WGCNA module figures
- PPI network
- Hub-gene table
- ML biomarker table
- ROC curves
- Immune infiltration figures
- Regulatory networks
- Drug-prediction results
- Manuscript-ready methods
