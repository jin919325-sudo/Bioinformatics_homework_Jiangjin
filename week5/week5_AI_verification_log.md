# Week 5 — AI verification log (Homework 1)

Agent used: **Hermes Agent** (Nous Research), model `deepseek-v4-flash`, running in a
terminal on the student's own Windows machine (R 4.6.1, DESeq2 1.52.0, apeglm 1.34.0,
Bioconductor 3.23, ggplot2 4.0.3). All numbers below come from the logs in
`verification/logs/`.

---

## 1. Prompt preserved (verbatim)

The assignment file `Week5_Homework_1_Instructions.md`, `Week5_Homework_Starter.R` and
`Week5_Homework_Verification_Checklist.md` were pasted in full as context, followed by
this instruction:

> Read the Week 5 instructions, starter and verification checklist. Write **one runnable
> R script** at `scripts/01_deseq2_analysis.R` that performs every required step in the
> stated order: import; assert that the count matrix contains only non-negative integers
> and that its columns are `identical()` to the metadata row names in the same order;
> make `condition` and `batch` factors with `control` as the reference level; build
> `DESeqDataSetFromMatrix` with `design = ~ batch + condition`; assert the design table
> is balanced and the model matrix full rank; pre-filter genes with >= 10 counts in
> >= 3 samples and report the retention; run `DESeq()`; **print `resultsNames(dds)` and
> verify the coefficient name instead of assuming it**; extract the treated-versus-control
> contrast; apply `apeglm` log2-fold-change shrinkage; export the **complete** result table
> including non-significant genes; save the fitted object and `sessionInfo()`; produce a
> PCA plot and a volcano plot as `week5_pca.png` and `week5_de_plot.png` at 300 dpi with
> readable labels. Thresholds are `padj < 0.05` **and** `|log2FoldChange| >= 1`. Do not
> transform the counts before DESeq2 and do not remove samples. Run it with `Rscript`,
> keep stdout in `logs/`, and then re-verify every number you report against that log.

## 2. What the AI generated

- `scripts/01_deseq2_analysis.R` — the full analysis (assertions, fit, shrinkage,
  figures, exports, sensitivity check, truth-key cross-check).
- `scripts/02_verification.R` — the sceptical re-checks described in section 4.
- The two figure files and the draft of `week5_interpretation.md`.

The AI also wrote the code for the two "extra" checks that the homework did not require
(batch-in/batch-out sensitivity, comparison against the provided annotation key); they
are reported because they are evidence, not because they were requested.

## 3. Package functions and arguments checked (three AI errors found and fixed)

| # | AI wrote | R error / symptom | Fix and evidence |
|---|---|---|---|
| 1 | `vst(dds, blind = FALSE)` | `less than 'nsub' rows, it is recommended to use varianceStabilizingTransformation directly` — the pre-filter leaves 989 genes, below the default `nsub = 1000` | `vst(dds, blind = FALSE, nsub = nrow(dds))`; the log prints `vst(): nsub = 989`. The transformation is still VST, using all genes instead of a 1,000-gene subset |
| 2 | `lfcShrink(dds, coef = ..., type = "apeglm", alpha = 0.05)` | `apeglm::apeglm(...) : 参数没有用(alpha = 0.05)` — apeglm 1.34.0 rejects `alpha` | LFC/lfcSE from `lfcShrink(type = "apeglm")`, hypothesis-test columns from `results(..., alpha = 0.05)`. The log also proves the choice is harmless here: `genes whose padj changed when alpha moved 0.1 -> 0.05: 0` (83 genes below 0.05 either way) |
| 3 | `keep <- rownames(res)[res$significant]` in the verification script | `rowMeans ... 下标出界` (subscript out of bounds) — the exported CSV has no row names, so `rownames()` returned `"1","2",...` | `keep <- res$gene_id[res$significant]`; the fixed run checks 60/60 genes |
| 4 | `cor(res_df$log2FoldChange, nb$log2FoldChange)` for the batch-in/out comparison | nonsense output (`r = 0.1818`, sign agreement 0.69) although both models called ~60 genes | `res_df` had already been **sorted by padj**, so the element-wise pairing compared different genes. Merged by `gene_id` instead: `r = 0.9977`, 58 shared genes, mean \|ΔLFC\| = 0.018 |
| 5 | sign agreement across all 989 genes against the annotation key | meaningless value (0.099) | the key contains `truth_log2FC == 0` for 900 genes and `sign(0)` matches nothing; the metric was restricted to the 98 genes with `truth != 0` → sign agreement 1.00 |

Also logged: `calcNormFactors()` is deprecated in edgeR >= 4 (`has been renamed to
normLibSizes`) — the verification script now calls `normLibSizes()` when it exists.

## 4. Independent verification of the biology and the numbers

| Check | How it was done (not via DESeq2's own output) | Result |
|---|---|---|
| Sample identity | the sample IDs were decoded from their own pattern (first letter `C`/`T` = condition, second letter = batch) and compared with the metadata columns | 12/12 match; design balanced 2/2/2 per batch |
| Coefficient direction | the sign of every significant gene's shrunken LFC was compared with the sign of the fold change of **library-size-normalised mean counts** computed directly from the CSV | 60/60 agree (r = 0.999); also 83/83 for the `padj < 0.05`-only set |
| Second method | the same filtered matrix and the same `~ batch + condition` design were re-run with **limma-voom** | 60 significant in DESeq2, 72 in voom, **all 60 DESeq2 calls are shared** (0 unique to DESeq2), r = 0.958 over all genes, 60/60 direction agreement |
| Reproducibility | the DESeq2 fit was rebuilt from scratch and compared with the saved table | max \|padj diff\| = 5.0e-16, max \|shrunken LFC diff\| = 4.9e-15, identical significant set = TRUE |
| Truth labels | the shipped `..._Gene_Annotation_Instructor_Key.csv` was joined **after** the analysis and never used to fit anything | 80/98 truth-DE genes detected at `padj < 0.05` (81.6%), 58/88 with the dual threshold (65.9%), 0 false positives among 891 truth-null genes, direction agreement 60/60 |

## 5. What the AI did **not** do / limitations stated honestly

- The AI did not invent any number: every figure caption, table and sentence in
  `week5_interpretation.md` is traceable to `verification/logs/01_deseq2_analysis.log`,
  `02_verification.log` or `week5_summary_stats.txt`.
- The dataset is a simulated 1,000-gene matrix with unannotated gene IDs (`Gene0001`…),
  so no pathway-level or mechanistic claim is made. A statistically significant gene is
  not claimed to be biologically important.
- No sample was removed on PCA position or library size, as the homework requires; the
  mean VST sample-to-sample distances are logged (21.5–22.3, i.e. no outlier sample).
- The 12 extra limma-voom calls are explained by shrinkage, not by disagreement: voom
  LFCs are unshrunken, so marginally strong genes exceed \|LFC\| >= 1 there.
