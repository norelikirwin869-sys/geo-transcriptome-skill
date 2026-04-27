# 02_deseq2_rnaseq.R
# RNA-seq count matrix differential-expression template using DESeq2

suppressPackageStartupMessages({
  library(DESeq2)
})

# -----------------------------
# Input files
# -----------------------------
count_file <- "data/count_matrix.csv"      # genes x samples
metadata_file <- "data/sample_metadata.csv" # columns: sample_id, group
out_dir <- "results/deseq2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load data
# -----------------------------
counts <- read.csv(count_file, row.names = 1, check.names = FALSE)
meta <- read.csv(metadata_file)
rownames(meta) <- meta$sample_id

counts <- counts[, rownames(meta)]
meta$group <- factor(meta$group, levels = c("Control", "Disease"))

# -----------------------------
# Filter low-expression genes
# -----------------------------
keep <- rowSums(counts >= 10) >= 2
counts <- counts[keep, ]

# -----------------------------
# DESeq2
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData = meta,
  design = ~ group
)

dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "Disease", "Control"))
res <- as.data.frame(res[order(res$padj), ])

write.csv(res, file.path(out_dir, "DESeq2_DEG.csv"))

# Variance-stabilized matrix for visualization and downstream analysis
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)
write.csv(vsd_mat, file.path(out_dir, "DESeq2_vst_expression.csv"))
