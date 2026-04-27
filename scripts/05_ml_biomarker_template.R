# 05_ml_biomarker_template.R
# Machine-learning biomarker screening template

suppressPackageStartupMessages({
  library(glmnet)
  library(randomForest)
  library(e1071)
  library(caret)
  library(pROC)
})

expr_file <- "results/deseq2/DESeq2_vst_expression.csv" # genes x samples
candidate_file <- "data/candidate_genes.txt"             # one gene symbol per line
metadata_file <- "data/sample_metadata.csv"              # sample_id, group
out_dir <- "results/ml"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
meta <- read.csv(metadata_file)
rownames(meta) <- meta$sample_id
candidates <- scan(candidate_file, what = character(), quiet = TRUE)

common_genes <- intersect(rownames(expr), candidates)
x <- t(expr[common_genes, rownames(meta), drop = FALSE])
y <- factor(meta$group, levels = c("Control", "Disease"))

set.seed(123)

# -----------------------------
# LASSO logistic regression
# -----------------------------
cvfit <- cv.glmnet(as.matrix(x), y, family = "binomial", alpha = 1)
coef_mat <- coef(cvfit, s = "lambda.min")
lasso_genes <- rownames(coef_mat)[as.numeric(coef_mat) != 0]
lasso_genes <- setdiff(lasso_genes, "(Intercept)")
writeLines(lasso_genes, file.path(out_dir, "lasso_genes.txt"))

# -----------------------------
# Random forest
# -----------------------------
rf_dat <- data.frame(group = y, x, check.names = FALSE)
rf <- randomForest(group ~ ., data = rf_dat, importance = TRUE, ntree = 1000)
rf_imp <- importance(rf)
rf_imp <- data.frame(gene = rownames(rf_imp), rf_imp, check.names = FALSE)
rf_imp <- rf_imp[order(rf_imp$MeanDecreaseGini, decreasing = TRUE), ]
write.csv(rf_imp, file.path(out_dir, "random_forest_importance.csv"), row.names = FALSE)

# -----------------------------
# ROC for LASSO genes
# -----------------------------
roc_summary <- data.frame()
for (gene in lasso_genes) {
  roc_obj <- roc(y, as.numeric(x[, gene]), levels = c("Control", "Disease"), quiet = TRUE)
  roc_summary <- rbind(roc_summary, data.frame(gene = gene, AUC = as.numeric(auc(roc_obj))))
}
write.csv(roc_summary, file.path(out_dir, "lasso_gene_auc.csv"), row.names = FALSE)

# -----------------------------
# Consensus candidates
# -----------------------------
top_rf <- head(rf_imp$gene, 30)
consensus <- intersect(lasso_genes, top_rf)
writeLines(consensus, file.path(out_dir, "consensus_biomarkers.txt"))
