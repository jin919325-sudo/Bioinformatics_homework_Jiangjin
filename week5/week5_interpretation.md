# Week 5 — Interpretation of the treated-versus-control DESeq2 result

**Comparison.** Twelve samples (two conditions, three balanced batches) were modelled as
`~ batch + condition`; the treated-versus-control contrast used apeglm shrinkage with
padj<0.05 and |log2FC|>=1.

**QC.** Library sizes are homogeneous (112,619–151,748 counts, max/min 1.35). PC1
(24% variance) tracks condition (R²=0.99) and PC2 (9%) tracks batch (R²=0.49);
no sample was removed.

**Result.** Sixty genes pass both thresholds — 36 up and 24 down in treated (83 at padj<0.05 alone), strongest Gene0008 (+1.93) and Gene0100 (−1.90).

**Meaning.** The consistent 36-up/24-down split suggests a coordinated transcriptional
response to treatment rather than scattered single-gene noise.

**Limitation.** Twelve samples, 1,000 simulated genes, no pathway annotation or
functional validation: these genes are statistically, not biologically, established.

**AI use.** AI drafted the script and figure code (prompt logged below). I verified
sample order, design rank, coefficient direction and every number against
`logs/01_deseq2_analysis.log`, which records two AI errors found and fixed.

<!-- word count (Word 口径, 不含标题): 143 -->
