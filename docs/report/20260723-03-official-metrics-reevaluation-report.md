# FLOWR.ROOT 官方核心指标复评报告

**实验 ID:** `20260723-03-official-metrics-reevaluation`

**对应计划:** [20260723-03-official-metrics-reevaluation-plan.md](../plan/20260723-03-official-metrics-reevaluation-plan.md)

**状态:** 已完成

## 1. 结论

不重新采样的前提下, 本实验使用 FLOWR.ROOT 官方 validity、完整 dock PoseBusters、MMFF94s strain 和条件子结构匹配, 复评两组各 50 个冻结 SDF.

- Coordinate-only: validity `50/50`, 条件匹配 `50/50`, PB-valid `14/50`.
- Official inpainting: validity `49/50`, 条件匹配 `49/50`, PB-valid `24/50`.
- Official inpainting 的 strain 中位数为 `145.80 kcal/mol`, 低于 coordinate-only 的 `7486.80 kcal/mol`, 但仍高于 5 个 `G_good` 的 `40.37 kcal/mol`.

按官方生成质量口径, official inpainting 更容易产生完整 PB-valid、低 strain 的局部新结构. 但官方条件匹配只要求固定子结构能在输出中找到, 不检查原始固定原子映射、坐标或完整分子身份, 因而 `24/50` 不能解释为 24 次严格修回原配体.

## 2. 协议与完整性

- 输入为两组完全相同的 5 个 case x 10 个 seed, 共 100 个生成结果; 未重新运行模型或过滤样本.
- Preflight 确认 100 个结果 SDF、10 个基线 SDF 和 5 个 protein PDB 全部存在并可读取.
- 完整 PB-valid 使用仓库 `posebusters/config/dock.yml` 的全部 20 个布尔检查项.
- 环境为 PoseBusters `0.3.1`、RDKit `2025.03.6`; strain 使用 MMFF94s、最多 500 次优化.
- 正式流程完整执行两次, 两次主要计数和 strain 汇总一致. 5 个针对性测试均通过.

输入清单与配置哈希见 [input_manifest.json](../../outputs/20260723-03-official-metrics-reevaluation/input_manifest.json) 和 [preflight.json](../../outputs/20260723-03-official-metrics-reevaluation/preflight.json).

## 3. 总体结果

| 指标 | Coordinate-only | Official inpainting |
|---|---:|---:|
| 尝试数 | 50 | 50 |
| RDKit valid | 50 | 49 |
| Fully-connected valid | 50 | 49 |
| 官方条件匹配 | 50 | 49 |
| 有效产物中的条件匹配率 | 50/50 | 49/49 |
| 完整 PB-valid | 14 | 24 |
| Conditioned PB-valid | 14 | 24 |
| Strain 可计算 | 50 | 49 |
| Strain 均值, kcal/mol | 28739.87 | 209.78 |
| Strain 中位数, kcal/mol | 7486.80 | 145.80 |
| Strain Q1-Q3, kcal/mol | 66.84-40461.14 | 113.16-190.97 |

Official inpainting 唯一的无效产物是已知的 4F0S seed `2026072349`, 原因是氮原子显式价态为 4. 它保留在 50 次分母中, 条件匹配和 strain 记为不可计算, PB-valid 为失败.

## 4. 逐 Case 与配对结果

| Case | Coordinate PB-valid | Official PB-valid | Coordinate strain median | Official strain median |
|---|---:|---:|---:|---:|
| 6PVZ | 10/10 | 7/10 | 10.53 | 61.78 |
| 7DDL | 0/10 | 6/10 | 79026.41 | 184.24 |
| 3ROG | 3/10 | 2/10 | 70.41 | 165.65 |
| 4BV5 | 0/10 | 7/10 | 21245.89 | 131.31 |
| 4F0S | 1/10 | 2/10 | 20035.30 | 148.84 |

相同 case/seed 的完整 PB-valid 配对为:

| 结果 | 数量 |
|---|---:|
| 两组都通过 | 7 |
| 仅 coordinate-only 通过 | 7 |
| 仅 official inpainting 通过 | 17 |
| 两组都失败 | 19 |

## 5. 主要失败项

| PB 子项失败数 | Coordinate-only | Official inpainting |
|---|---:|---:|
| Bond lengths | 29 | 12 |
| Bond angles | 29 | 16 |
| Internal steric clash | 26 | 1 |
| Internal energy | 29 | 5 |
| Aromatic ring flatness | 8 | 1 |
| Minimum distance to protein | 10 | 3 |

Coordinate-only 的主要瓶颈仍是局部共价几何和应变. 它固定离散图, 但没有把可编辑环或支链作为刚性结构处理, 因而 7DDL、4BV5 和 4F0S 出现极高 strain. Official inpainting 可以重新生成原子和键, 因而整体几何更容易落入官方容差.

## 6. 与旧 native/strict 的关系

Coordinate-only 的 PB-valid `14/50` 与旧 native `14/50` 恰好完全重合. Official inpainting 则为:

| 旧 native / 官方 PB-valid | 数量 |
|---|---:|
| 两者都通过 | 13 |
| 仅官方 PB-valid 通过 | 11 |
| 仅旧 native 通过 | 1 |
| 两者都失败 | 25 |

差异来自两套标准本身:

- 官方条件检查做子结构搜索, 不要求原始固定原子逐一映射. 7DDL 的 10 个有效产物全部通过官方条件匹配, 其中 6 个完整 PB-valid; 这不会推翻旧评估中固定边界原子芳香性发生改变的事实.
- 官方 dock 配置对键长和键角采用 25% 容差, 对 internal clash 采用 30% 容差; 旧评估器直接调用模块默认值, 三项均为 20%.
- 官方完整 PB 额外检查芳香环/双键平面性和内部能量. 因此有 1 个旧 native 产物因为 internal energy 未通过而不是官方 PB-valid.

## 7. 基线与边界

5 个 `G_good` 全部 PB-valid, strain 中位数为 `40.37 kcal/mol`. 5 个 `G_bad` 全部因 protein clash 不是 PB-valid, strain 中位数为 `50.79 kcal/mol`, 说明人工扭转主要制造蛋白质冲突, 没有像 coordinate-only 的部分输出那样严重破坏内部几何.

本实验只有 5 个案例, 不进行显著性检验. 未计算 Vina、唯一性或多样性, 也未使用官方筛选后重采样. 官方 PB-valid 衡量生成质量, 不检查同分子身份; 旧 strict `7/50` 仍是判断 official inpainting 能否直接充当严格同分子坐标修复器的正式结果.

逐样本结果见 [runs.csv](../../outputs/20260723-03-official-metrics-reevaluation/runs.csv), 基线见 [baselines.csv](../../outputs/20260723-03-official-metrics-reevaluation/baselines.csv), 完整汇总见 [summary.json](../../outputs/20260723-03-official-metrics-reevaluation/summary.json).
