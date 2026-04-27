# 03_enrichment_gsea.R
# GO, KEGG and Hallmark GSEA template

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(enrichplot)
})

# -----------------------------
# Input
# -----------------------------
deg_file <- "results/limma/GSEXXXX_limma_DEG.csv"
out_dir <- "results/enrichment"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# DEG table should contain gene symbols in row names or a SYMBOL column and logFC/log2FoldChange.
deg <- read.csv(deg_file, row.names = 1, check.names = FALSE)

if (!"SYMBOL" %in% colnames(deg)) {
  deg$SYMBOL <- rownames(deg)
}

fc_col <- if ("logFC" %in% colnames(deg)) "logFC" else "log2FoldChange"
p_col <- if ("adj.P.Val" %in% colnames(deg)) "adj.P.Val" else "padj"

sig <- deg[!is.na(deg[[p_col]]) & deg[[p_col]] < 0.05 & abs(deg[[fc_col]]) > 0.5, , drop = FALSE]

# -----------------------------
# ID conversion
# -----------------------------
gene_df <- bitr(sig$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# -----------------------------
# GO enrichment
# -----------------------------
ego_bp <- enrichGO(
  gene = gene_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

write.csv(as.data.frame(ego_bp), file.path(out_dir, "GO_BP.csv"))

# -----------------------------
# KEGG enrichment
# -----------------------------
ekk <- enrichKEGG(
  gene = gene_df$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

write.csv(as.data.frame(ekk), file.path(out_dir, "KEGG.csv"))

# -----------------------------
# Hallmark GSEA
# -----------------------------
ranks <- deg[[fc_col]]
names(ranks) <- deg$SYMBOL
ranks <- ranks[!is.na(ranks)]
ranks <- sort(ranks, decreasing = TRUE)

hallmark <- msigdbr(species = "Homo sapiens", category = "H")
term2gene <- hallmark[, c("gs_name", "gene_symbol")]

gsea_h <- GSEA(
  geneList = ranks,
  TERM2GENE = term2gene,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

write.csv(as.data.frame(gsea_h), file.path(out_dir, "Hallmark_GSEA.csv"))
