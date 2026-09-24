# Week 6 · 16S summary table

数据：EMP `tests/16S_level-7.csv`（470 特征 × 132 样本）+ `tests/16S_mapping.csv`（130 行）。
所有数字由 `week6_16s_diversity.R` 从原始文件重算，日志见 `logs/week6_16s_diversity.log.txt`。

## 1. 关键指标

| 项目 | 结果 |
|---|---|
| counts table (input) | 470 features x 132 samples |
| metadata table (input) | 130 rows (SampleID, Group, Group_sub) |
| samples without a Group label (excluded from group tests) | 2 (K_XYL_F_0009_03, K_XYL_F_0035_03) |
| features after prevalence filter >= 0.10 (detect >= 0.05) | 132 |
| features removed by the filter | 338 |
| features in the analysis state (top 40) | 40 |
| sequencing depth, analysis state (min / median / max) | 1624 / 9761 / 20033 |
| features detected per sample (median / range) | 17 (2 - 27) |
| Shannon, all 132 samples (mean / median) | 1.494 / 1.563 |
| Shannon, 130 labelled samples (mean / median) | 1.489 / 1.562 |
| Kruskal-Wallis p, Shannon over 4 groups | 0.279 |
| PERMANOVA R2 / F / p — Bray-Curtis ~ Group, 4 groups, n = 130 | 0.0338 / 1.47 / 0.088 |
| betadisper permutest p, 4 groups | 0.618 |
| PERMANOVA R2 / F / p — UC_before vs IBS_before, n = 65 | 0.0265 / 1.71 / 0.101 |
| betadisper permutest p, UC_before vs IBS_before | 0.336 |
| differential features tested | 40 |
| differential p<0.05 / FDR<0.05 | 3 / 0 |
| minimum BH FDR | 0.39 |

## 2. 测序深度（每个分组）

| Group | n | depth min | depth median | depth max |
|---|---|---|---|---|
| IBS_after | 36 | 1624 | 9584.5 | 17839 |
| IBS_before | 36 | 5393 | 10041 | 15862 |
| UC_after | 29 | 5340 | 10742 | 20033 |
| UC_before | 29 | 4210 | 9409 | 15624 |

## 3. Alpha 多样性（Shannon，分组）

| Group | n | Shannon mean | Shannon median | 检出特征数 median |
|---|---|---|---|---|
| IBS_before | 36 | 1.548 | 1.646 | 17.5 |
| IBS_after | 36 | 1.505 | 1.625 | 17.5 |
| UC_before | 29 | 1.501 | 1.570 | 15 |
| UC_after | 29 | 1.383 | 1.443 | 16 |

## 4. 差异丰度（Wilcoxon，UC_before vs IBS_before），按 raw p 排序前 8

| feature | p | BH FDR | median IBS_before | median UC_before | 偏高的一侧 |
|---|---|---|---|---|---|
| s__torques | 0.0181 | 0.390 | 0.0 | 0.0 | UC_before |
| s__mucosae | 0.0288 | 0.390 | 0.0 | 0.0 | UC_before |
| __ | 0.0293 | 0.390 | 2548.0 | 1719.0 | IBS_before |
| s__copri | 0.0655 | 0.465 | 19.5 | 0.0 | IBS_before |
| s__formicilis | 0.0668 | 0.465 | 101.0 | 0.0 | IBS_before |
| s__plebeius | 0.0745 | 0.465 | 0.0 | 0.0 | UC_before |
| s__bromii | 0.0906 | 0.465 | 0.0 | 0.0 | UC_before |
| s__longum | 0.0930 | 0.465 | 88.0 | 109.0 | UC_before |

## 5. 与 EMP 提交结果的一致性核对

| 核对项 | 结果 |
|---|---|
| Shannon per sample, this script vs EMP export (65 shared samples) | mean 1.527 vs 1.527 | max abs diff 0.0000 | Pearson r 1.0000 |
| Wilcoxon p / BH FDR per feature vs EMP export | n = 40 | max abs p diff 0.00000 | max abs FDR diff 0.00000 |
