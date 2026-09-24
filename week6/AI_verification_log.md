# Week 6 · AI 使用与核验记录（AI verification log）

本次作业的 AI 交互由我在 EMP-Web 里的 AI 代理（Hermes agent）完成，下面记录：给它的 prompt、
用到的输出、我核实／修正／驳回的条目，以及它在执行过程中自行修正并被日志记录下来的问题。

## 1. Prompt 摘要（给 AI 的要求）

- 「照课程内置指南《Demo 操作指南：16S Microbiome（全模块流程）》的参数顺序，在 EMP-Web 上跑完
  16S 全流程：过滤 min prevalence 0.10 / min detect 0.05 → 16S Taxonomy 聚合 Genus Top 40 →
  Alpha（Shannon）→ 降维（Bray-Curtis PCoA）→ 可视化 → rclr → 差异分析 wilcox
  (UC_before vs IBS_before)；每一步的参数和返回都要留日志。」
- 「补上 rubric 里要求的证据：过滤前后的输入维度与样本深度、PERMANOVA 与离散度检验、
  100–150 词英文解读、AI 验证日志 + sessionInfo；代码要能从干净 session 直接跑。
  不要做指南以外的额外分析。」
- 「把课程页『解读与假设』的四段文字按我的口吻重写，不要模板腔。」

## 2. 用到的 AI 输出

- EMP-Web 上 8 个分析步骤的请求参数与返回：`logs/01–17*.json`（过滤 132/338、Top 40、Alpha、
  PCoA 坐标、rclr、差异表、快照列表）。
- 课程页「解读与假设」四段文字（已随 EMP 同步进仓库的 `teaching/report.md`）。
- 独立重算的 R 脚本 `week6_16s_diversity.R`、汇总表 `tables/week6_summary_table.csv`、
  三张图 `figures/*.pdf`。

## 3. 我核实／修正／驳回的条目

| # | 类型 | 内容 | 证据 |
|---|---|---|---|
| 1 | 核实并采纳 | AI 报「过滤后保留 132 特征、删除 338；Top 40 后 40 个特征；Shannon 均值 1.494 / 中位 1.563」。我用 `week6_16s_diversity.R` 从原始 CSV 独立重算，并逐样本与 EMP 导出的 `m16s_course_alpha.csv` 对齐：**最大绝对差 0.0000，Pearson r = 1.0000**；40 个特征的 wilcoxon p 与 BH FDR 与 EMP 导出表**最大绝对差 0.00000** | `tables/week6_verification_vs_EMP.csv` |
| 2 | 修正 | AI 第一版流程顺序为「Alpha → 降维 → rclr → 差异」。核对日志时发现 EMP 的差异分析会重置会话对象并把样本裁到参与比较的两组，先算的 Alpha/降维结果会从会话里消失。改为「差异 → Alpha → 降维 → rclr」重放后，三个结果表才同时保留在会话里 | `logs/11–17*.json`、预处理快照列表 |
| 3 | 驳回 | EMP 拒绝在全部 132 样本上做 wilcox（分组因子仍带着未选中的水平，报 `Column Group has beed deteced missing value`），AI 提议改用 Kruskal-Wallis 多组检验顶替。**不采纳**：指南与评分表都要求两组 wilcox 比较，最终按 EMP 默认的成对子集（UC_before 29 + IBS_before 36 = 65 样本）执行 | `logs/12_diff_wilcox.json` |
| 4 | 驳回 | AI 第一版「解读与假设」使用「参数链」「可证伪假设」这类模板化表达。**要求按我的口吻重写**，第二版（第一人称、短句、保留不确定性和局限性）才采用 | journal 两次更新：`logs/18_teaching_journal.json` → `logs/27_teaching_journal_humanized.json` |
| 5 | 补充 | AI 报告：EMP 没有独立的 PERMANOVA 接口，PCoA 散点图里的 PERMANOVA 标注也没有渲染出来，所以第一次提交里只有图形没有 PERMANOVA 数值。我要求补上，最后在 R 里用 `vegan::adonis2` + `vegan::betadisper` 计算并写进汇总表与图中 | `tables/week6_summary_table.csv`、`figures/week6_pcoa_bray.pdf`、`figures/week6_dispersion.pdf` |

## 4. 版本记录

- R 4.6.1（`vegan`、`ggplot2`、`readr`、`dplyr` 版本见 `logs/sessionInfo.txt`）
- EMP / EasyMultiProfiler 9.0.4（`manifest.json` 里的 `emp_version`）
- 脚本：`week6_16s_diversity.R`（一次 `Rscript` 跑完，无手工步骤）
