# 04_wgcna_template.R
# WGCNA co-expression network template

suppressPackageStartupMessages({
  library(WGCNA)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

expr_file <- "results/deseq2/DESeq2_vst_expression.csv" # genes x samples
trait_file <- "data/sample_metadata.csv"                # columns: sample_id, group
out_dir <- "results/wgcna"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
traits <- read.csv(trait_file)
rownames(traits) <- traits$sample_id

# WGCNA expects samples x genes
datExpr <- as.data.frame(t(expr))
datExpr <- datExpr[rownames(traits), ]

# Filter low-variance genes
vars <- apply(datExpr, 2, var, na.rm = TRUE)
datExpr <- datExpr[, vars > quantile(vars, 0.25, na.rm = TRUE)]

gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

trait_df <- data.frame(
  Disease = ifelse(traits[rownames(datExpr), "group"] == "Disease", 1, 0)
)
rownames(trait_df) <- rownames(datExpr)

powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

softPower <- 6
net <- blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
moduleTraitCor <- cor(MEs, trait_df, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

write.csv(moduleTraitCor, file.path(out_dir, "module_trait_cor.csv"))
write.csv(moduleTraitPvalue, file.path(out_dir, "module_trait_pvalue.csv"))
write.csv(data.frame(gene = colnames(datExpr), module = moduleColors), file.path(out_dir, "gene_modules.csv"), row.names = FALSE)
