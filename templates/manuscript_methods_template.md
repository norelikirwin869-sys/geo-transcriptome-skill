# Manuscript Methods Template

## GEO dataset acquisition

Gene-expression datasets were retrieved from the Gene Expression Omnibus database. Search terms included disease-related keywords, tissue-specific keywords and transcriptome-related keywords. Datasets were screened according to predefined inclusion and exclusion criteria. Studies were included when they contained clearly defined disease and control samples, relevant tissue or cell sources, accessible expression matrices or raw count data, and adequate platform annotation.

## Data preprocessing

For microarray datasets, expression matrices were downloaded and normalized. Probe identifiers were mapped to official gene symbols using the corresponding platform annotation. Probes without valid gene symbols were removed, and multiple probes corresponding to the same gene were collapsed. For RNA-seq datasets, raw count matrices were used when available. Low-expression genes were filtered before normalization and differential-expression analysis.

## Differential-expression analysis

Differentially expressed genes between disease and control groups were identified using limma for microarray data or DESeq2 for RNA-seq count data. P values were adjusted using the Benjamini-Hochberg method. Genes with adjusted P value < 0.05 and absolute log2 fold change > 0.5 were considered differentially expressed, unless otherwise specified.

## Functional enrichment analysis

Gene Ontology and pathway enrichment analyses were performed using clusterProfiler. Enriched biological processes, cellular components, molecular functions and KEGG pathways were identified. Gene set enrichment analysis was conducted using ranked gene lists and curated gene sets such as Hallmark, KEGG or Reactome.

## Weighted gene co-expression network analysis

Weighted gene co-expression network analysis was performed to identify disease-associated gene modules. After filtering low-variance genes, a scale-free co-expression network was constructed. Module eigengenes were correlated with disease traits, and genes from disease-associated modules were selected for further analysis.

## PPI network and hub-gene identification

Protein-protein interaction networks were constructed using public interaction resources such as STRING. Hub genes were prioritized according to topological algorithms including Degree, MCC, MNC, EPC, Closeness and Betweenness. Candidate hub genes were further integrated with differential-expression and WGCNA results.

## Machine-learning biomarker screening

Candidate biomarkers were screened using machine-learning algorithms including LASSO logistic regression, random forest and SVM-RFE. Consensus genes identified by multiple methods were selected as diagnostic biomarker candidates. Diagnostic performance was evaluated using receiver operating characteristic curves and area under the curve values.

## Immune infiltration analysis

Immune-cell infiltration patterns were estimated using bulk transcriptome deconvolution or single-sample enrichment methods such as ssGSEA, xCell or MCP-counter. Differences in immune-cell scores between disease and control groups were assessed, and correlations between biomarker genes and immune-cell infiltration were calculated.

## Validation

Candidate biomarkers and pathways were validated using independent datasets when available. Additional experimental validation, such as qRT-PCR, western blotting, immunohistochemistry, immunofluorescence, animal models or cell models, was recommended to confirm biological significance.
