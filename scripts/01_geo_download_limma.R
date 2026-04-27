# 01_geo_download_limma.R
# GEO microarray download and limma differential-expression template

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
})

# -----------------------------
# User settings
# -----------------------------
gse_id <- "GSEXXXX"
case_pattern <- "Disease"
control_pattern <- "Control"
out_dir <- "results/limma"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Download GEO matrix
# -----------------------------
gse <- getGEO(gse_id, GSEMatrix = TRUE)
eset <- gse[[1]]
expr <- exprs(eset)
pheno <- pData(eset)

# Optional log2 transform
qx <- as.numeric(quantile(expr, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
log_needed <- qx[5] > 100 || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (log_needed) {
  expr[expr <= 0] <- NA
  expr <- log2(expr)
}

expr <- normalizeBetweenArrays(expr)

# -----------------------------
# Define groups manually after inspecting pheno
# -----------------------------
print(colnames(pheno))
print(head(pheno[, seq_len(min(6, ncol(pheno)))]))

# Replace this with project-specific grouping.
# Example:
# group <- ifelse(grepl(case_pattern, pheno$characteristics_ch1, ignore.case = TRUE), "Disease", "Control")
# group <- factor(group, levels = c("Control", "Disease"))

stop("Edit group assignment before running differential-expression analysis.")

# -----------------------------
# Differential expression
# -----------------------------
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
contrast <- makeContrasts(Disease - Control, levels = design)

fit <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

deg <- topTable(fit2, adjust.method = "BH", number = Inf)
write.csv(deg, file.path(out_dir, paste0(gse_id, "_limma_DEG.csv")))

# Save normalized expression
write.csv(expr, file.path(out_dir, paste0(gse_id, "_normalized_expression.csv")))
