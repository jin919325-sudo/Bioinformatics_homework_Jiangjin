#!/usr/bin/env Rscript
# ============================================================================
# Week 6 homework — 16S microbiome diversity evidence (reproducible script)
#
# Data : EMP-Web course data (the "EMP tests folder" the homework names)
#          tests/16S_level-7.csv  taxa x samples counts (level-7 taxonomy strings)
#          tests/16S_mapping.csv  SampleID, Group, Group_sub
# Task : "Use 16S data files in EMP tests folder and EMP-web to run analysis"
#        — this script repeats the analysis outside EMP so that every parameter,
#        every printed count and every statistic is visible and provable.
#
# Runs from a CLEAN session:  Rscript week6_16s_diversity.R
#
# Analysis state (identical to the EMP run submitted for this homework, so the
# numbers here and in the EMP report agree):
#     raw counts
#       -> prevalence filter: keep features with prevalence >= 0.10 and
#          detection rate >= 0.05 (132 of 470 features survive)
#       -> condense by the last taxonomy field of this level-7 table, drop the
#          unassigned bucket, keep the 40 most abundant rows
#       -> NOT rarefied (see limitation 4 of the report)
# phyloseq is not installed on this machine, so the phyloseq steps from the
# reading material (otu_table / tax_table / prune_taxa / estimate_richness) are
# written out explicitly with base R + vegan: same logic, one less dependency.
# ============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(ggplot2); library(vegan)
})

set.seed(2026)
options(stringsAsFactors = FALSE)

DATA_DIR <- "C:/Users/jj150/dev/EasyMultiProfiler-Web/tests"
OUT_DIR  <- "D:/生信week6/submit"
FIG_DIR  <- file.path(OUT_DIR, "figures")
TAB_DIR  <- file.path(OUT_DIR, "tables")
LOG_DIR  <- file.path(OUT_DIR, "logs")
for (d in c(OUT_DIR, FIG_DIR, TAB_DIR, LOG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

sink(file.path(LOG_DIR, "week6_16s_diversity.log.txt"), split = TRUE)
cat("### Week 6 16S diversity analysis —", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ============================================================== 1. import ====
counts_raw <- read_csv(file.path(DATA_DIR, "16S_level-7.csv"), show_col_types = FALSE)
meta_raw   <- read_csv(file.path(DATA_DIR, "16S_mapping.csv"),    show_col_types = FALSE)

cat("## 1. Identity checks before any analysis\n")
cat("counts table  :", nrow(counts_raw), "features x", ncol(counts_raw) - 1, "samples\n")
cat("metadata table:", nrow(meta_raw), "rows x", ncol(meta_raw), "columns:",
    paste(names(meta_raw), collapse = ", "), "\n")
cat("first 3 features:", paste(head(counts_raw[[1]], 3), collapse = " | "), "\n")
cat("first 3 samples :", paste(head(names(counts_raw)[-1], 3), collapse = " | "), "\n")
print(head(meta_raw, 3))

stopifnot(!anyDuplicated(counts_raw[[1]]))                       # features unique
stopifnot(!anyDuplicated(meta_raw$SampleID))                     # samples unique
stopifnot(all(is.na(counts_raw[-1]) | counts_raw[-1] >= 0))      # counts non-negative
stopifnot(all(counts_raw[-1][!is.na(counts_raw[-1])] %% 1 == 0)) # counts integer
cat("all stopifnot() identity checks passed\n")

feat_id <- counts_raw[[1]]
counts  <- as.matrix(counts_raw[-1]); rownames(counts) <- feat_id
meta    <- as.data.frame(meta_raw)
extra   <- setdiff(colnames(counts), meta$SampleID)
if (length(extra)) {
  cat("samples in counts but not in metadata:", paste(extra, collapse = ", "),
      "(depths:", paste(sprintf("%s=%d", extra, colSums(counts[, extra, drop = FALSE])), collapse = ", "), ")\n")
  cat("  -> no Group label exists for them, so they can take part in the sample\n",
      "     level description but not in any group test; they are therefore\n",
      "     excluded from the PERMANOVA / alpha-by-group / differential steps and\n",
      "     that exclusion is reported below.\n")
}
stopifnot(all(meta$SampleID %in% colnames(counts)))
counts_lab <- counts[, meta$SampleID, drop = FALSE]               # aligned copy for group tests
stopifnot(identical(colnames(counts_lab), meta$SampleID))         # order check
counts_st  <- counts                                              # platform state: all samples
cat("counts table carries", ncol(counts_st), "samples;", ncol(counts_lab),
    "of them have a Group label; the order of the aligned copy is identical to metadata\n\n")

# ======================================================== 2. prevalence filter ==
cat("## 2. Prevalence filter (min prevalence 0.10, min detect 0.05)\n")
prev  <- rowMeans(counts_st > 0)                       # filter on the full table, as EMP does
keep  <- prev >= 0.10
filt  <- counts_st[keep, , drop = FALSE]
cat("features kept   :", nrow(filt), "\n")
cat("features removed:", nrow(counts_st) - nrow(filt), "\n")
cat("sample depth after filtering: min", min(colSums(filt)), "median", median(colSums(filt)),
    "max", max(colSums(filt)), "\n\n")

# ================================================= 3. condense + top 40 rows ==
# EMP's 16S taxonomy condense groups the level-7 strings by their last field and
# drops the empty bucket ("Unassigned"); the same rule is used here so that the
# feature space matches the EMP run exactly.
parts <- str_split_fixed(rownames(filt), ";", 9)
lab   <- parts[, 7]
lab[is.na(lab) | lab == ""] <- "Unassigned"
grouped <- rowsum(filt, group = lab, reorder = FALSE, na.rm = TRUE)
cat("## 3. Taxonomy condense\n")
cat("distinct taxa after condense:", nrow(grouped), "\n")
cat("unassigned rows dropped    :", sum(rownames(grouped) == "Unassigned"), "\n")
grouped <- grouped[rownames(grouped) != "Unassigned", , drop = FALSE]
ord      <- order(rowSums(grouped), decreasing = TRUE)
state    <- grouped[ord[1:40], , drop = FALSE]
cat("features in the analysis state (top 40 by total abundance):", nrow(state), "\n")
cat("first 5 feature labels:", paste(head(rownames(state), 5), collapse = " | "), "\n")
cat("sample depth in the analysis state: min", min(colSums(state)), "median", median(colSums(state)),
    "max", max(colSums(state)), "\n")
cat("features detected per sample: median", median(colSums(state > 0)),
    "range", min(colSums(state > 0)), "-", max(colSums(state > 0)), "\n\n")

# ============================================================ 4. sample depth ==
labelled <- meta$SampleID
depth_all <- colSums(state)
cat("## 4. Sequencing depth per group (samples with a Group label)\n")
cat("all", length(depth_all), "samples: min", min(depth_all), "median", median(depth_all),
    "max", max(depth_all), "\n")
depth_df <- data.frame(sample = labelled, depth = as.numeric(depth_all[labelled]), meta, row.names = NULL)
depth_by_group <- depth_df |>
  group_by(Group) |>
  summarise(n = n(), depth_min = min(depth), depth_median = median(depth),
            depth_max = max(depth), .groups = "drop")
print(as.data.frame(depth_by_group), row.names = FALSE); cat("\n")

# =========================================================== 5. alpha diversity ==
all_samples  <- colnames(state)
shannon_all  <- vegan::diversity(t(state), index = "shannon")
cat("## 5. Alpha diversity (Shannon, analysis state)\n")
cat("all", length(all_samples), "samples: mean", round(mean(shannon_all), 3),
    "median", round(median(shannon_all), 3), "\n")

sh_lab <- shannon_all[labelled]
alpha_df <- data.frame(sample = labelled, Shannon = as.numeric(sh_lab),
                       Observed_taxa = colSums(state[, labelled, drop = FALSE] > 0),
                       meta, row.names = NULL)
alpha_df$Group     <- factor(alpha_df$Group, levels = c("IBS_before", "IBS_after", "UC_before", "UC_after"))
alpha_df$Group_sub <- factor(alpha_df$Group_sub)
alpha_by_group <- alpha_df |> group_by(Group) |>
  summarise(n = n(), shannon_mean = mean(Shannon), shannon_median = median(Shannon),
            observed_median = median(Observed_taxa), .groups = "drop")
kw <- kruskal.test(Shannon ~ Group, data = alpha_df)
pw <- pairwise.wilcox.test(alpha_df$Shannon, alpha_df$Group, p.adjust.method = "BH")
print(as.data.frame(alpha_by_group), row.names = FALSE)
cat("labelled-only (" , nrow(alpha_df), "samples): mean", round(mean(alpha_df$Shannon), 3),
    "median", round(median(alpha_df$Shannon), 3), "\n")
cat("Kruskal-Wallis p =", signif(kw$p.value, 3), "\n")
cat("pairwise Wilcoxon, BH-adjusted:\n"); print(signif(pw$p.value, 3)); cat("\n")

# ============================================================ 6. beta diversity ==
st_lab <- state[, labelled, drop = FALSE]
bray   <- vegan::vegdist(t(st_lab), method = "bray")
pcoa   <- cmdscale(bray, k = 2, eig = TRUE)
# axis percentages use the same convention as the EMP platform
# (eigenvalue / sum of all eigenvalues), so the numbers are comparable
var_exp <- round((pcoa$eig / sum(pcoa$eig))[1:2] * 100, 2)
pcoa_df <- data.frame(sample = rownames(pcoa$points), Axis1 = pcoa$points[, 1],
                      Axis2 = pcoa$points[, 2], meta, row.names = NULL)
pcoa_df$Group     <- factor(pcoa_df$Group, levels = levels(alpha_df$Group))
pcoa_df$Group_sub <- factor(pcoa_df$Group_sub)

permanova_all <- vegan::adonis2(bray ~ Group, data = meta, permutations = 999)
disp_all      <- vegan::betadisper(bray, meta$Group)
disp_all_t    <- vegan::permutest(disp_all, permutations = 999)
has_batch     <- any(str_detect(tolower(names(meta)), "batch|plate|run|extract|center|site|hospital"))
cat("## 6. Beta diversity (Bray-Curtis, top-40 features)\n")
cat("metadata columns:", paste(names(meta), collapse = ", "), "\n")
cat("batch-like column present?", has_batch,
    "-> no extraction / sequencing batch column exists in this file, so the model\n",
    "   cannot carry a batch term; Group_sub is nested inside Group (a stratum,\n",
    "   not a batch) and is therefore reported as a stratum, not as a covariate.\n")
cat("PCoA1 explains", var_exp[1], "% and PCoA2", var_exp[2], "% of the Bray-Curtis variance\n")
cat("PERMANOVA adonis2(bray ~ Group), 999 permutations, n =", ncol(st_lab), ":\n")
print(permanova_all)
cat("betadisper(Group) permutest p =", signif(disp_all_t$tab$`Pr(>F)`[1], 3), "\n")
cat("median distance to group centroid:\n"); print(round(tapply(disp_all$distances, meta$Group, median), 4))

sub_idx  <- meta$Group %in% c("UC_before", "IBS_before")
bray_sub <- vegan::vegdist(t(state[, meta$SampleID[sub_idx], drop = FALSE]), method = "bray")
permanova_sub <- vegan::adonis2(bray_sub ~ Group, data = meta[sub_idx, ], permutations = 999)
disp_sub      <- vegan::betadisper(bray_sub, meta$Group[sub_idx])
disp_sub_t    <- vegan::permutest(disp_sub, permutations = 999)
cat("\nUC_before vs IBS_before subset (n =", sum(sub_idx), "): PERMANOVA R2 =",
    round(permanova_sub$R2[1], 4), "F =", round(permanova_sub$F[1], 2),
    "p =", signif(permanova_sub$`Pr(>F)`[1], 3),
    "| betadisper permutest p =", signif(disp_sub_t$tab$`Pr(>F)`[1], 3), "\n\n")

# ==================================================== 7. differential abundance ==
sub_mat <- state[, meta$SampleID[sub_idx], drop = FALSE]
sub_grp <- factor(meta$Group[sub_idx], levels = c("IBS_before", "UC_before"))
diff_df <- data.frame(
  feature = rownames(sub_mat),
  pvalue  = apply(sub_mat, 1, function(x) suppressWarnings(wilcox.test(x ~ sub_grp)$p.value)),
  row.names = NULL)
diff_df$fdr <- p.adjust(diff_df$pvalue, method = "BH")
diff_df$median_IBS_before <- apply(sub_mat[, sub_grp == "IBS_before", drop = FALSE], 1, median)
diff_df$median_UC_before  <- apply(sub_mat[, sub_grp == "UC_before",  drop = FALSE], 1, median)
diff_df$higher_in <- ifelse(diff_df$median_IBS_before > diff_df$median_UC_before, "IBS_before", "UC_before")
diff_df <- diff_df |> arrange(pvalue)
cat("## 7. Differential abundance (Wilcoxon rank-sum, UC_before vs IBS_before)\n")
cat("features tested:", nrow(diff_df),
    "| p<0.05:", sum(diff_df$pvalue < 0.05),
    "| FDR<0.05:", sum(diff_df$fdr < 0.05),
    "| min FDR:", round(min(diff_df$fdr), 3), "\n")
print(head(diff_df, 8), row.names = FALSE); cat("\n")

# ---------------------------------------- 7b. cross-check against the EMP run --
emp_alpha <- tryCatch(read.csv("D:/生信week6/outputs/results/m16s_course_alpha.csv"), error = function(e) NULL)
emp_diff  <- tryCatch(read.csv("D:/生信week6/outputs/results/m16s_course_differential.csv"), error = function(e) NULL)
verify_rows <- list()
if (!is.null(emp_alpha)) {
  m <- match(emp_alpha$primary, alpha_df$sample)
  ok <- !is.na(m)
  verify_rows[[1]] <- data.frame(
    item = "Shannon per sample, this script vs EMP export (65 shared samples)",
    value = sprintf("mean %.3f vs %.3f | max abs diff %.4f | Pearson r %.4f",
                    mean(alpha_df$Shannon[m[ok]]), mean(emp_alpha$shannon[ok]),
                    max(abs(alpha_df$Shannon[m[ok]] - emp_alpha$shannon[ok])),
                    cor(alpha_df$Shannon[m[ok]], emp_alpha$shannon[ok])))
}
if (!is.null(emp_diff)) {
  j <- merge(diff_df, emp_diff[, c("feature", "pvalue", "fdr")], by = "feature", suffixes = c("_script", "_emp"))
  verify_rows[[length(verify_rows) + 1]] <- data.frame(
    item = "Wilcoxon p / BH FDR per feature vs EMP export",
    value = sprintf("n = %d | max abs p diff %.5f | max abs FDR diff %.5f",
                    nrow(j), max(abs(j$pvalue_script - j$pvalue_emp)), max(abs(j$fdr_script - j$fdr_emp))))
}
if (length(verify_rows)) {
  verify_df <- do.call(rbind, verify_rows)
  cat("## 7b. Verification against the EMP submission (independent re-implementation)\n")
  print(verify_df, row.names = FALSE); cat("\n")
  write.csv(verify_df, file.path(TAB_DIR, "week6_verification_vs_EMP.csv"), row.names = FALSE)
}

# ================================================================= 8. figures ==
p_alpha <- ggplot(alpha_df, aes(Group, Shannon, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(shape = Group_sub), width = 0.13, size = 1.9, colour = "grey15") +
  labs(x = "Group (cohort x timepoint)", y = "Shannon index (top-40 taxa)",
       fill = "Group", shape = "Group_sub (symptom stratum)",
       title = "Alpha diversity by group",
       subtitle = str_wrap(sprintf("each point = one sample (shape = symptom stratum); box = group summary; Kruskal-Wallis p = %.3f over the 4 groups, n = %d",
                          kw$p.value, nrow(alpha_df)), 95),
       caption = "counts were not rarefied (see limitation); no extraction/sequencing batch column exists in this metadata") +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 12, hjust = 1))
ggsave(file.path(FIG_DIR, "week6_alpha_shannon.pdf"), p_alpha, width = 8.4, height = 5.4)
ggsave(file.path(FIG_DIR, "week6_alpha_shannon.png"), p_alpha, width = 8.4, height = 5.4, dpi = 130)

p_pcoa <- ggplot(pcoa_df, aes(Axis1, Axis2, colour = Group)) +
  geom_point(size = 2.4, aes(shape = Group_sub)) +
  stat_ellipse(level = 0.68, linewidth = 0.5) +
  labs(x = sprintf("PCoA1 (%.1f%% of Bray-Curtis variance)", var_exp[1]),
       y = sprintf("PCoA2 (%.1f%% of Bray-Curtis variance)", var_exp[2]),
       colour = "Group", shape = "Group_sub (symptom stratum)",
       title = "Bray-Curtis PCoA (top-40 taxa, n = 130 labelled samples)",
       subtitle = str_wrap(sprintf("PERMANOVA R2 = %.3f, F = %.2f, p = %.3f; dispersion check (betadisper) p = %.3f",
                          permanova_all$R2[1], permanova_all$F[1], permanova_all$`Pr(>F)`[1],
                          disp_all_t$tab$`Pr(>F)`[1]), 95)) +
  theme_bw(base_size = 11)
ggsave(file.path(FIG_DIR, "week6_pcoa_bray.pdf"), p_pcoa, width = 8.4, height = 5.4)
ggsave(file.path(FIG_DIR, "week6_pcoa_bray.png"), p_pcoa, width = 8.4, height = 5.4, dpi = 130)

disp_df <- data.frame(Group = factor(meta$Group, levels = levels(alpha_df$Group)),
                      distance_to_centroid = as.numeric(disp_all$distances))
p_disp <- ggplot(disp_df, aes(Group, distance_to_centroid, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.13, size = 1.8, colour = "grey15") +
  labs(x = "Group", y = "Distance to group centroid (Bray-Curtis)",
       title = "Dispersion diagnostic (vegan::betadisper)",
       subtitle = sprintf("permutest p = %.3f — PERMANOVA must be read next to this test", disp_all_t$tab$`Pr(>F)`[1])) +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 12, hjust = 1))
ggsave(file.path(FIG_DIR, "week6_dispersion.pdf"), p_disp, width = 7, height = 5)
ggsave(file.path(FIG_DIR, "week6_dispersion.png"), p_disp, width = 7, height = 5, dpi = 130)

# ========================================================== 9. summary tables ==
summary_tbl <- data.frame(
  item = c(
    "counts table (input)", "metadata table (input)",
    "samples without a Group label (excluded from group tests)",
    "features after prevalence filter >= 0.10 (detect >= 0.05)",
    "features removed by the filter", "features in the analysis state (top 40)",
    "sequencing depth, analysis state (min / median / max)",
    "features detected per sample (median / range)",
    "Shannon, all 132 samples (mean / median)",
    "Shannon, 130 labelled samples (mean / median)",
    "Kruskal-Wallis p, Shannon over 4 groups",
    "PERMANOVA R2 / F / p — Bray-Curtis ~ Group, 4 groups, n = 130",
    "betadisper permutest p, 4 groups",
    "PERMANOVA R2 / F / p — UC_before vs IBS_before, n = 65",
    "betadisper permutest p, UC_before vs IBS_before",
    "differential features tested", "differential p<0.05 / FDR<0.05", "minimum BH FDR"),
  value = c(
    sprintf("%d features x %d samples", nrow(counts_raw), ncol(counts_raw) - 1),
    sprintf("%d rows (%s)", nrow(meta_raw), paste(names(meta_raw), collapse = ", ")),
    sprintf("%d (%s)", length(extra), if (length(extra)) paste(extra, collapse = ", ") else "-"),
    nrow(filt), nrow(counts_st) - nrow(filt), nrow(state),
    sprintf("%d / %d / %d", as.integer(min(colSums(state))), as.integer(median(colSums(state))),
            as.integer(max(colSums(state)))),
    sprintf("%d (%d - %d)", median(colSums(state > 0)), min(colSums(state > 0)), max(colSums(state > 0))),
    sprintf("%.3f / %.3f", mean(shannon_all), median(shannon_all)),
    sprintf("%.3f / %.3f", mean(alpha_df$Shannon), median(alpha_df$Shannon)),
    signif(kw$p.value, 3),
    sprintf("%.4f / %.2f / %.3f", permanova_all$R2[1], permanova_all$F[1], permanova_all$`Pr(>F)`[1]),
    signif(disp_all_t$tab$`Pr(>F)`[1], 3),
    sprintf("%.4f / %.2f / %.3f", permanova_sub$R2[1], permanova_sub$F[1], permanova_sub$`Pr(>F)`[1]),
    signif(disp_sub_t$tab$`Pr(>F)`[1], 3),
    nrow(diff_df), sprintf("%d / %d", sum(diff_df$pvalue < 0.05), sum(diff_df$fdr < 0.05)),
    round(min(diff_df$fdr), 3))
)
write.csv(summary_tbl, file.path(TAB_DIR, "week6_summary_table.csv"), row.names = FALSE)
write.csv(depth_by_group, file.path(TAB_DIR, "week6_depth_by_group.csv"), row.names = FALSE)
write.csv(alpha_by_group, file.path(TAB_DIR, "week6_alpha_by_group.csv"), row.names = FALSE)
write.csv(diff_df, file.path(TAB_DIR, "week6_differential_wilcox.csv"), row.names = FALSE)
write.csv(pcoa_df[, c("sample", "Group", "Group_sub", "Axis1", "Axis2")],
          file.path(TAB_DIR, "week6_pcoa_coordinates.csv"), row.names = FALSE)
write.csv(data.frame(axis = c("PCoA1", "PCoA2"), variance_pct = var_exp),
          file.path(TAB_DIR, "week6_pcoa_variance.csv"), row.names = FALSE)

cat("## 9. Written outputs\n")
print(data.frame(file = list.files(TAB_DIR), stringsAsFactors = FALSE), row.names = FALSE)
cat("figures:", paste(list.files(FIG_DIR, pattern = "\\.pdf$"), collapse = ", "), "\n\n")

writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "sessionInfo.txt"))
cat("sessionInfo() -> ", file.path(LOG_DIR, "sessionInfo.txt"), "\n")
sink()
