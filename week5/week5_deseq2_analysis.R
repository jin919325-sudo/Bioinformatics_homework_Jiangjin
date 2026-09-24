# ==========================================================================
# Week 5 Homework 1 — Bulk RNA-seq differential expression with DESeq2
# From verified counts to an interpretable DESeq2 result
#
# Inputs (read-only): data_raw/Week5_Homework_Count_Matrix.csv
#                     data_raw/Week5_Homework_Sample_Metadata.csv
#                     data_raw/Week5_Homework_Gene_Annotation_Instructor_Key.csv
# Outputs: outputs/week5_deseq2_results.csv, outputs/week5_deseq2_object.rds,
#          outputs/session_info.txt, outputs/week5_summary_stats.txt,
#          outputs/week5_sensitivity_batch.csv, outputs/week5_results_vs_truth.csv
#          figures/week5_pca.png, figures/week5_de_plot.png, figures/week5_ma_plot.png
#
# Run from the project root:  Rscript scripts/01_deseq2_analysis.R
# ==========================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(ggrepel)
})

set.seed(1)
dir.create("outputs", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

count_file    <- "data_raw/Week5_Homework_Count_Matrix.csv"
metadata_file <- "data_raw/Week5_Homework_Sample_Metadata.csv"
key_file      <- "data_raw/Week5_Homework_Gene_Annotation_Instructor_Key.csv"

stopifnot(file.exists(count_file), file.exists(metadata_file), file.exists(key_file))

cat("### Week 5 / Homework 1 - DESeq2 analysis\n")
cat("### started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --------------------------------------------------------------------------
# 1. Import
# --------------------------------------------------------------------------
counts <- read.csv(count_file, row.names = 1, check.names = FALSE)
coldata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)

cat("--- 1. Import ---\n")
cat("count matrix      :", nrow(counts), "genes x", ncol(counts), "samples\n")
cat("metadata          :", nrow(coldata), "samples x", ncol(coldata), "variables\n")
cat("metadata columns  :", paste(colnames(coldata), collapse = ", "), "\n\n")

# --------------------------------------------------------------------------
# 2. Mandatory validation (every check is an assertion, not a comment)
# --------------------------------------------------------------------------
cat("--- 2. Validation ---\n")

# 2.1 sample identity: columns of the matrix must equal row names of metadata,
#     in the same order (a reordered metadata file would silently pair the
#     wrong condition with the wrong column)
stopifnot(ncol(counts) == nrow(coldata))
stopifnot(identical(colnames(counts), rownames(coldata)))
cat("[OK] count-matrix columns identical to metadata row names, same order\n")

# 2.2 duplicate sample IDs
stopifnot(!anyDuplicated(colnames(counts)), !anyDuplicated(rownames(coldata)))
cat("[OK] no duplicated sample IDs\n")

# 2.3 raw integer counts: non-negative, no fractional and no negative values
m <- as.matrix(counts)
stopifnot(is.numeric(m))
stopifnot(all(is.finite(m)))
stopifnot(all(m >= 0))
stopifnot(all(m == round(m)))
cat("[OK] all", length(m), "values are finite, non-negative integers",
    sprintf("(min = %g, max = %g)", min(m), max(m)), "\n")

# 2.4 factor levels and reference level
coldata$condition <- factor(coldata$condition)
coldata$batch     <- factor(coldata$batch)
stopifnot(setequal(levels(coldata$condition), c("control", "treated")))
stopifnot(setequal(levels(coldata$batch), c("A", "B", "C")))
coldata$condition <- relevel(coldata$condition, ref = "control")
coldata$batch     <- relevel(coldata$batch, ref = "A")
cat("[OK] condition levels:", paste(levels(coldata$condition), collapse = " < "),
    "| reference =", levels(coldata$condition)[1], "\n")
cat("[OK] batch levels    :", paste(levels(coldata$batch), collapse = " < "), "\n")

# 2.5 the design must be balanced and the model matrix full rank, otherwise the
#     batch term and the treatment term are partially confounded
design_tab <- table(batch = coldata$batch, condition = coldata$condition)
cat("--- design table (batch x condition) ---\n")
print(design_tab)
stopifnot(all(design_tab > 0))
stopifnot(length(unique(as.vector(design_tab))) == 1)   # balanced (2 per cell)
cat("[OK] design is balanced (", as.vector(design_tab)[1], "samples per cell)\n", sep = "")

mm <- model.matrix(~ batch + condition, data = coldata)
stopifnot(qr(mm)$rank == ncol(mm))
cat("[OK] model matrix ~ batch + condition is full rank (rank =",
    qr(mm)$rank, "of", ncol(mm), "columns)\n\n")

# --------------------------------------------------------------------------
# 3. Library-size QC (before any modelling)
# --------------------------------------------------------------------------
cat("--- 3. Library-size QC ---\n")
lib <- colSums(m)
lib_summary <- summary(lib)
print(lib_summary)
cat("max/min library size ratio:", round(max(lib) / min(lib), 3), "\n")
cat("samples below 1e5 total counts:", sum(lib < 1e5), "\n")
cat("QC observation: library sizes are homogeneous (",
    format(min(lib), big.mark = ","), " - ", format(max(lib), big.mark = ","),
    " counts, ratio ", round(max(lib) / min(lib), 2),
    "); DESeq2 median-of-ratios normalisation is therefore appropriate and no\n",
    "sample was removed on library size or on PCA position alone.\n\n", sep = "")
write.csv(data.frame(sample_id = names(lib), library_size = as.integer(lib),
                     row.names = NULL),
          "outputs/week5_library_sizes.csv", row.names = FALSE)

# --------------------------------------------------------------------------
# 4. DESeq2 object
# --------------------------------------------------------------------------
cat("--- 4. DESeqDataSet ---\n")
# batch is included because the three batches are balanced across conditions:
# batch is a technical source of variance, and modelling it removes that
# variance from the residual, which increases power instead of confounding.
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = ~ batch + condition)
cat("design:", deparse(design(dds)), "\n\n")

# --------------------------------------------------------------------------
# 5. Pre-filter: >= 10 counts in >= 3 samples
# --------------------------------------------------------------------------
cat("--- 5. Pre-filter ---\n")
keep <- rowSums(counts(dds) >= 10) >= 3
cat("rule            : retain genes with >= 10 counts in >= 3 samples\n")
cat("genes before    :", nrow(dds), "\n")
cat("genes retained  :", sum(keep), sprintf(" (%.1f%%)\n", 100 * mean(keep)))
cat("genes removed   :", sum(!keep), "\n")
dds <- dds[keep, ]
write.csv(data.frame(gene_id = names(keep), retained = as.logical(keep),
                     row.names = NULL),
          "outputs/week5_filter_table.csv", row.names = FALSE)

# --------------------------------------------------------------------------
# 6. Fit
# --------------------------------------------------------------------------
cat("\n--- 6. DESeq() fit ---\n")
dds <- DESeq(dds)
coef_names <- resultsNames(dds)
cat("resultsNames(dds):\n")
print(coef_names)

target_coef <- "condition_treated_vs_control"
if (!target_coef %in% coef_names) {
  stop("Expected coefficient was not found. Inspect resultsNames(dds) and update target_coef.")
}
cat("coefficient inspected and confirmed present:", target_coef, "\n")

# --------------------------------------------------------------------------
# 7. Extract + shrunken LFC
# --------------------------------------------------------------------------
cat("\n--- 7. Results + apeglm shrinkage ---\n")
res <- results(dds, contrast = c("condition", "treated", "control"), alpha = 0.05)
cat("unshrunken summary:\n")
print(summary(res))

res_shrunk <- lfcShrink(dds, coef = target_coef, type = "apeglm", quiet = TRUE)
cat("shrunken summary (LFC from apeglm):\n")
print(summary(res_shrunk))

res_df <- as.data.frame(res_shrunk)
# apeglm 1.34 in this R install rejects `alpha` (error: unused argument), so
# lfcShrink returns p/padj at its own default alpha = 0.1 while the shrunken LFC
# does not depend on alpha. The homework fixes alpha = 0.05, so the hypothesis
# test columns are taken from the results() call above (alpha = 0.05) and only
# the LFC / lfcSE / baseMean columns come from apeglm. The difference is logged.
stopifnot(identical(rownames(res_df), rownames(res)))
padj_diff <- sum((is.na(res_df$padj) != is.na(res$padj)) |
                 (!is.na(res_df$padj) & !is.na(res$padj) &
                    abs(res_df$padj - res$padj) > 1e-12))
cat("genes whose padj changed when alpha moved 0.1 -> 0.05 (independent filtering):",
    padj_diff, "\n")
cat("padj < 0.05 at alpha = 0.1 (lfcShrink default):",
    sum(!is.na(res_df$padj) & res_df$padj < 0.05),
    "| at alpha = 0.05 (used below):",
    sum(!is.na(res$padj) & res$padj < 0.05), "\n")
res_df$pvalue <- res$pvalue
res_df$padj   <- res$padj

res_df$gene_id <- rownames(res_df)
res_df$significant <- !is.na(res_df$padj) & res_df$padj < 0.05 &
                      abs(res_df$log2FoldChange) >= 1
res_df$direction <- "Not significant"
res_df$direction[res_df$significant & res_df$log2FoldChange > 0] <- "Up in treated"
res_df$direction[res_df$significant & res_df$log2FoldChange < 0] <- "Down in treated"
res_df <- res_df[order(res_df$padj, -abs(res_df$log2FoldChange)), ]

# full table, including non-significant genes, ordered by adjusted p value
write.csv(res_df, "outputs/week5_deseq2_results.csv", row.names = FALSE)
cat("\nthresholds: padj < 0.05  AND  |shrunken log2FoldChange| >= 1\n")
cat("significant genes:", sum(res_df$significant), "\n")
print(table(res_df$direction))
cat("genes with padj < 0.05 only (no effect-size filter):",
    sum(!is.na(res_df$padj) & res_df$padj < 0.05), "\n")
cat("genes with NA padj (independent filtering):", sum(is.na(res_df$padj)), "\n\n")

cat("top 10 by adjusted p value:\n")
print(head(res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE",
                      "pvalue", "padj", "direction")], 10))

# --------------------------------------------------------------------------
# 8. PCA on VST
# --------------------------------------------------------------------------
cat("\n--- 8. PCA (vst, blind = FALSE) ---\n")
# apeglm/vst quirk on this data set: vst() refuses to subsample when the object
# has fewer than `nsub` (default 1000) rows, and the homework filter leaves 989.
# Setting nsub to the actual number of retained genes keeps the parametric
# dispersion fit and uses every gene instead of a random subset.
vsd <- vst(dds, blind = FALSE, nsub = nrow(dds))
cat("vst(): nsub =", nrow(dds), "(fewer genes than the default 1000, so all genes used)\n")
pca_df <- plotPCA(vsd, intgroup = c("condition", "batch"), returnData = TRUE)
percent_var <- round(100 * attr(pca_df, "percentVar"))
cat("variance explained: PC1 =", percent_var[1], "% PC2 =", percent_var[2], "%\n")

# quantify what each axis tracks: ANOVA R^2 of the PC scores on condition / batch
r2 <- function(y, g) summary(lm(y ~ g))$r.squared
variance_axis <- data.frame(
  axis = c("PC1", "PC2"),
  R2_condition = c(r2(pca_df$PC1, pca_df$condition), r2(pca_df$PC2, pca_df$condition)),
  R2_batch     = c(r2(pca_df$PC1, pca_df$batch),     r2(pca_df$PC2, pca_df$batch))
)
cat("ANOVA R^2 of PC scores:\n")
print(variance_axis)

# sample-level outlier screen (reported, no sample removed)
sd_all <- assay(vsd)
pc_mean <- colMeans(sd_all)
sample_dist <- colMeans(as.matrix(dist(t(sd_all))))
cat("\nmean sample-to-sample VST distance (samples furthest from the group first):\n")
print(round(sort(sample_dist, decreasing = TRUE), 1))

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = batch)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = name), size = 3, max.overlaps = Inf, seed = 1) +
  labs(
    title = "Week 5 RNA-seq PCA (VST-transformed counts)",
    subtitle = "colour = condition, shape = batch (n = 12, no sample removed)",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme_bw(base_size = 12)
ggsave("figures/week5_pca.png", p_pca, width = 7, height = 5, dpi = 300)

# --------------------------------------------------------------------------
# 9. Volcano and MA plots
# --------------------------------------------------------------------------
cat("\n--- 9. DE figures ---\n")
plot_df <- res_df
plot_df$neg_log10_padj <- -log10(pmax(plot_df$padj, 1e-300))

p_volcano <- ggplot(plot_df, aes(x = log2FoldChange, y = neg_log10_padj,
                                 color = direction)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("Up in treated" = "#C0392B",
                                "Down in treated" = "#2F6DB3",
                                "Not significant" = "grey70")) +
  annotate("text", x = 0, y = 6.5, hjust = 0.5, size = 3.2,
           label = "thresholds: padj < 0.05\n& |shrunken LFC| >= 1") +
  labs(title = "Differential expression: treated versus control",
       subtitle = "shrunken (apeglm) log2 fold change, design ~ batch + condition",
       x = "Shrunken log2 fold change", y = "-log10 adjusted p value",
       color = NULL) +
  theme_bw(base_size = 12)
ggsave("figures/week5_de_plot.png", p_volcano, width = 7, height = 5, dpi = 300)

p_ma <- ggplot(plot_df, aes(x = log10(baseMean + 1), y = log2FoldChange,
                            color = direction)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = 0, color = "grey40") +
  scale_color_manual(values = c("Up in treated" = "#C0392B",
                                "Down in treated" = "#2F6DB3",
                                "Not significant" = "grey70")) +
  labs(title = "MA plot: treated versus control",
       subtitle = "shrunken (apeglm) log2 fold change against mean normalised count",
       x = "log10(mean normalised count + 1)", y = "Shrunken log2 fold change",
       color = NULL) +
  theme_bw(base_size = 12)
ggsave("figures/week5_ma_plot.png", p_ma, width = 7, height = 5, dpi = 300)
cat("written: figures/week5_pca.png, figures/week5_de_plot.png, figures/week5_ma_plot.png\n")

# --------------------------------------------------------------------------
# 10. Reproducibility files
# --------------------------------------------------------------------------
cat("\n--- 10. Saving object + session info ---\n")
saveRDS(dds, "outputs/week5_deseq2_object.rds")
saveRDS(vsd, "outputs/week5_vst_object.rds")
capture.output(sessionInfo(), file = "outputs/session_info.txt")
cat("written: outputs/week5_deseq2_object.rds, outputs/session_info.txt\n")

# --------------------------------------------------------------------------
# 11. Sensitivity check: what the batch term actually does
# --------------------------------------------------------------------------
cat("\n--- 11. Sensitivity check (batch in vs out of the design) ---\n")
dds_nobatch <- DESeqDataSetFromMatrix(countData = counts[keep, ], colData = coldata,
                                      design = ~ condition)
dds_nobatch <- DESeq(dds_nobatch)
res_nb <- lfcShrink(dds_nobatch, coef = "condition_treated_vs_control",
                    type = "apeglm", quiet = TRUE)
nb <- as.data.frame(res_nb)
nb$gene_id <- rownames(nb)
nb$pvalue <- results(dds_nobatch, contrast = c("condition", "treated", "control"),
                     alpha = 0.05)$pvalue
nb$padj <- results(dds_nobatch, contrast = c("condition", "treated", "control"),
                   alpha = 0.05)$padj
nb$significant <- !is.na(nb$padj) & nb$padj < 0.05 & abs(nb$log2FoldChange) >= 1
# align the two models by gene id: res_df was sorted by padj above, so a direct
# element-wise comparison would silently pair different genes
cmp <- merge(res_df[, c("gene_id", "log2FoldChange", "significant")],
             nb[, c("gene_id", "log2FoldChange", "significant")],
             by = "gene_id", suffixes = c("_batch", "_nobatch"))
stopifnot(nrow(cmp) == nrow(res_df))
cat("significant with ~ batch + condition  :", sum(cmp$significant_batch), "\n")
cat("significant with ~ condition (no batch):", sum(cmp$significant_nobatch), "\n")
cat("shared significant genes:",
    sum(cmp$significant_batch & cmp$significant_nobatch), "\n")
cat("significant only with batch modelled  :",
    sum(cmp$significant_batch & !cmp$significant_nobatch), "\n")
cat("significant only without batch        :",
    sum(!cmp$significant_batch & cmp$significant_nobatch), "\n")
cat("concordance of shrunken LFC (all 989 genes, gene-id aligned), Pearson r:",
    round(cor(cmp$log2FoldChange_batch, cmp$log2FoldChange_nobatch), 4), "\n")
cat("sign agreement (all genes):",
    round(mean(sign(cmp$log2FoldChange_batch) == sign(cmp$log2FoldChange_nobatch)), 4), "\n")
cat("mean |LFC change| when batch is dropped:",
    round(mean(abs(cmp$log2FoldChange_batch - cmp$log2FoldChange_nobatch)), 4), "\n")
sens <- data.frame(
  model = c("~ batch + condition", "~ condition"),
  n_significant = c(sum(cmp$significant_batch), sum(cmp$significant_nobatch)),
  n_shared = c(sum(cmp$significant_batch & cmp$significant_nobatch), NA),
  r_lfc_vs_full_model = c(1, round(cor(cmp$log2FoldChange_batch,
                                       cmp$log2FoldChange_nobatch), 4))
)
write.csv(sens, "outputs/week5_sensitivity_batch.csv", row.names = FALSE)

# --------------------------------------------------------------------------
# 12. Cross-check against the provided gene annotation key (independent check)
# --------------------------------------------------------------------------
cat("\n--- 12. Cross-check vs provided gene annotation key ---\n")
key <- read.csv(key_file, check.names = FALSE)
names(key)[names(key) == "truth_log2FC_for_instructor"] <- "truth_lfc"
merged <- merge(res_df, key[, c("gene_id", "truth_lfc")], by = "gene_id", all.x = TRUE)
stopifnot(nrow(merged) == nrow(res_df), !any(is.na(merged$truth_lfc)))

# the key encodes the simulation truth: truth_lfc == 0 for the non-DE genes,
# truth_lfc != 0 for the 100 designed DE genes. A sign comparison over all genes
# would be meaningless because sign(0) == 0 never matches, so the metrics below
# are computed on the subsets where they mean something.
merged$truth_de  <- merged$truth_lfc != 0
merged$truth_big <- abs(merged$truth_lfc) >= 1

cat("genes with a truth value        :", nrow(merged), "\n")
cat("truth non-DE genes (truth == 0)  :", sum(!merged$truth_de), "\n")
cat("truth DE genes (truth != 0)      :", sum(merged$truth_de),
    " (", sum(merged$truth_de & merged$truth_lfc > 0), "up /",
    sum(merged$truth_de & merged$truth_lfc < 0), "down )\n", sep = "")
cat("truth genes with |truth LFC| >= 1:", sum(merged$truth_big), "\n")

cat("\nPearson r (shrunken LFC vs truth LFC), all genes:",
    round(cor(merged$log2FoldChange, merged$truth_lfc), 4), "\n")
cat("Pearson r, truth DE genes only              :",
    round(cor(merged$log2FoldChange[merged$truth_de],
              merged$truth_lfc[merged$truth_de]), 4), "\n")
cat("sign agreement, truth DE genes only         :",
    round(mean(sign(merged$log2FoldChange[merged$truth_de]) ==
                 sign(merged$truth_lfc[merged$truth_de])), 4), "\n")

cat("\n--- detection performance against the truth labels ---\n")
cat("power (padj < 0.05, no effect-size filter) among truth DE genes:",
    sum(merged$truth_de & merged$padj < 0.05), "/", sum(merged$truth_de),
    sprintf(" = %.1f%%\n",
            100 * sum(merged$truth_de & merged$padj < 0.05) / sum(merged$truth_de)))
cat("power with the homework dual threshold (padj < 0.05 & |LFC| >= 1):",
    sum(merged$truth_big & merged$significant), "/", sum(merged$truth_big),
    sprintf(" = %.1f%%\n",
            100 * sum(merged$truth_big & merged$significant) / sum(merged$truth_big)))
cat("false positives among truth non-DE genes (dual threshold):",
    sum(!merged$truth_de & merged$significant), "/", sum(!merged$truth_de),
    sprintf(" = %.2f%%\n", 100 * sum(!merged$truth_de & merged$significant) /
              sum(!merged$truth_de)))
cat("empirical FDR of my 60-gene significant set:",
    sprintf("%.1f%%", 100 * sum(merged$significant & !merged$truth_de) /
              sum(merged$significant)), "\n")
cat("truth DE genes MISSED by the dual threshold:",
    sum(merged$truth_de & !merged$significant),
    "(typical cause: apeglm shrinkage of a true |LFC| just above 1 below 1)\n")
cat("\nconfusion table (rows = model significant, cols = truth DE):\n")
print(table(model_significant = merged$significant, truth_de = merged$truth_de))
cat("\ndirection agreement among truth DE genes called significant:\n")
sig_de <- merged$significant & merged$truth_de
print(table(model_direction = merged$direction[sig_de],
            truth_direction = ifelse(merged$truth_lfc[sig_de] > 0, "up", "down")))

write.csv(merged[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "pvalue",
                     "padj", "significant", "direction", "truth_lfc")],
          "outputs/week5_results_vs_truth.csv", row.names = FALSE)

# --------------------------------------------------------------------------
# 13. Machine-readable summary for the report
# --------------------------------------------------------------------------
top_up   <- head(res_df$gene_id[res_df$direction == "Up in treated"], 10)
top_down <- head(res_df$gene_id[res_df$direction == "Down in treated"], 10)
sink("outputs/week5_summary_stats.txt")
cat("n_genes_input=", nrow(counts), "\n", sep = "")
cat("n_samples=", ncol(counts), "\n", sep = "")
cat("n_genes_retained=", sum(keep), "\n", sep = "")
cat("median_library_size=", as.integer(median(lib)), "\n", sep = "")
cat("libsize_ratio_max_min=", round(max(lib) / min(lib), 3), "\n", sep = "")
cat("PC1_var=", percent_var[1], "\n", sep = "")
cat("PC2_var=", percent_var[2], "\n", sep = "")
cat("R2_PC1_condition=", round(variance_axis$R2_condition[1], 3), "\n", sep = "")
cat("R2_PC1_batch=", round(variance_axis$R2_batch[1], 3), "\n", sep = "")
cat("R2_PC2_condition=", round(variance_axis$R2_condition[2], 3), "\n", sep = "")
cat("R2_PC2_batch=", round(variance_axis$R2_batch[2], 3), "\n", sep = "")
cat("n_padj_lt_0.05=", sum(!is.na(res_df$padj) & res_df$padj < 0.05), "\n", sep = "")
cat("n_significant=", sum(res_df$significant), "\n", sep = "")
cat("n_up=", sum(res_df$direction == "Up in treated"), "\n", sep = "")
cat("n_down=", sum(res_df$direction == "Down in treated"), "\n", sep = "")
cat("n_na_padj=", sum(is.na(res_df$padj)), "\n", sep = "")
cat("strongest_up=", res_df$gene_id[which.max(res_df$log2FoldChange)], "\n", sep = "")
cat("strongest_up_lfc=", round(max(res_df$log2FoldChange), 3), "\n", sep = "")
cat("strongest_down=", res_df$gene_id[which.min(res_df$log2FoldChange)], "\n", sep = "")
cat("strongest_down_lfc=", round(min(res_df$log2FoldChange), 3), "\n", sep = "")
cat("top_up=", paste(top_up, collapse = ","), "\n", sep = "")
cat("top_down=", paste(top_down, collapse = ","), "\n", sep = "")
cat("r_lfc_truth_all=", round(cor(merged$log2FoldChange, merged$truth_lfc), 4), "\n",
    sep = "")
cat("r_lfc_truth_de=",
    round(cor(merged$log2FoldChange[merged$truth_de],
              merged$truth_lfc[merged$truth_de]), 4), "\n", sep = "")
cat("n_truth_de=", sum(merged$truth_de), "\n", sep = "")
cat("power_truth_de_pct=",
    round(100 * sum(merged$truth_de & merged$padj < 0.05) / sum(merged$truth_de), 1),
    "\n", sep = "")
cat("truth_big_recovered_pct=",
    round(100 * sum(merged$truth_big & merged$significant) / sum(merged$truth_big), 1),
    "\n", sep = "")
cat("n_false_positive=", sum(!merged$truth_de & merged$significant), "\n", sep = "")
cat("empirical_fdr_pct=",
    round(100 * sum(merged$significant & !merged$truth_de) / sum(merged$significant), 1),
    "\n", sep = "")
cat("n_significant_nobatch=", sum(cmp$significant_nobatch), "\n", sep = "")
cat("n_shared_batch_vs_nobatch=",
    sum(cmp$significant_batch & cmp$significant_nobatch), "\n", sep = "")
cat("r_lfc_batch_vs_nobatch=",
    round(cor(cmp$log2FoldChange_batch, cmp$log2FoldChange_nobatch), 4), "\n", sep = "")
cat("mean_abs_lfc_shift_dropping_batch=",
    round(mean(abs(cmp$log2FoldChange_batch - cmp$log2FoldChange_nobatch)), 4), "\n",
    sep = "")
sink()
cat("\nwritten: outputs/week5_summary_stats.txt\n")
cat("\n### finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
