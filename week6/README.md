# Week 6 homework · 16S 微生物组（补充交付包）

作业要求（课程仓库 `Week 6/Homework/Homework-Week6.docx` 与 `Week 6/PPT` 末页）：

1. 用 EMP `tests/` 里的 16S 数据在 EMP-Web 上完成完整分析流程，并通过 EMP-Web 系统提交最终结果；
2. 基于参数和分析结果给出科学假设，并说明用哪些参数支撑该假设。

第 1、2 条已经在 EMP-Web 上完成并同步（仓库 `EMP2026/Week_06/microbiome_16s/weekly/runs/`）。
这个文件夹是**补充交付物**：把同一套分析在 EMP 之外用 R 重新实现一遍，并补齐阅读材料
第 30 页评分表要求、但 EMP 提交里没有的内容（R 脚本、汇总表、PERMANOVA 与离散度、英文解读、
AI 验证日志、sessionInfo）。

## 文件清单

| 文件 | 对应评分项 | 说明 |
|---|---|---|
| `week6_16s_diversity.R` | R script（25） | 从干净 session 一次跑完：identity checks、每步过滤/变换后的维度与深度、Alpha、PCoA、PERMANOVA、betadisper、差异分析、出图、`sessionInfo()`。`Rscript week6_16s_diversity.R` |
| `week6_summary_table.md` / `tables/week6_summary_table.csv` | Summary table（15） | sample depth、保留特征数、Shannon、PERMANOVA（R²/F/p）、离散度（betadisper p），以及差异分析汇总 |
| `figures/week6_alpha_shannon.pdf` | Alpha-diversity figure（15） | 单样本点 + 分组箱线 + 症状分层（Group_sub）+ Kruskal-Wallis 检验标注 + 可读轴标签 |
| `figures/week6_pcoa_bray.pdf` | Bray-Curtis PCoA（15） | 排序图 + 组内椭圆，副标题写明 PERMANOVA R²/F/p 与离散度 p |
| `interpretation_EN.md` | Interpretation（20） | 142 词英文解读：效应、不确定性、离散度、局限性 |
| `AI_verification_log.md` + `logs/sessionInfo.txt` | AI verification log + sessionInfo（10） | prompt 摘要、用到的输出、核实/修正/驳回的条目、版本 |
| `figures/week6_dispersion.pdf` | 证据链「PERMANOVA + 离散度」 | betadisper 距离箱线图 + permutest p |
| `logs/week6_16s_diversity.log.txt` | 证据链「每步过滤后的深度与特征数」 | 脚本的完整 stdout：导入维度、首个样本/特征标识、过滤前后、condense 前后、深度、统计量 |
| `tables/*.csv` | 证据 | 分组深度、分组 Alpha、PCoA 坐标与轴方差、差异表、与 EMP 的一致性核对 |

## 分析状态（与 EMP 提交完全一致，便于对照）

原始计数 → 过滤 min prevalence 0.10 / min detect 0.05（470 → 132 特征）→ 按该 level-7 表的最后
一级分类字段合并、丢掉未注释桶、按总丰度保留 40 个特征 → **不抽平**（局限性里写明）。
样本层面：132 个样本里 130 个有 Group 标签，`K_XYL_F_0009_03`、`K_XYL_F_0035_03` 两个样本没有分组，
只参与描述性统计，不进入任何组间检验（脚本里显式记录）。

## 关键数字

- 保留特征 132 / 删除 338；分析状态 40 个特征；样本深度 min 1624 / median 9761.5 / max 20033
- Shannon（全部 132 样本）均值 1.494、中位 1.563
- PCoA1 27.6%、PCoA2 19.3%（与 EMP 的 27.5% / 19.2% 同一算法约定）
- PERMANOVA `adonis2(bray ~ Group)`，999 次置换，n = 130：**R² = 0.034，F = 1.47，p = 0.088**
- 离散度 `betadisper` + `permutest`：**p = 0.618**（组间离散度无差异）
- 差异分析（wilcox，UC_before 29 vs IBS_before 36）：40 个特征，**p<0.05 3 个，FDR<0.05 0 个，最小 FDR 0.39**
- 与 EMP 提交的一致性：逐样本 Shannon 最大绝对差 0.0000（r = 1.0000）；40 个特征的 p / FDR 最大绝对差 0.00000

## 复现

```bash
Rscript week6_16s_diversity.R          # R 4.6.1 + vegan/ggplot2/readr/dplyr/tidyr/stringr
```
输入路径写在脚本顶部（`DATA_DIR` → EMP 的 `tests/` 目录），输出写到 `figures/`、`tables/`、`logs/`。
