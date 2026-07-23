# SPINDR 全测试集 Oracle Repair Benchmark 报告

**实验 ID:** `20260723-04-full-test-oracle-repair-benchmark`

**对应计划:** [20260723-04-full-test-oracle-repair-benchmark-plan.md](../plan/20260723-04-full-test-oracle-repair-benchmark-plan.md)

**状态:** 已完成

## 1. 结论

本实验在 SPINDR 官方 test split 的 225 个复合物上自动构造中等强度扭转 clash. 其中 211 个复合物满足冻结的构造条件, 每个案例使用 10 个固定 seed, 公平比较 coordinate-only 与 official inpainting, 共完成 4,220 次 100-step 采样.

- 严格同分子修复: coordinate-only 为 `584/2110` (`27.68%`), official inpainting 为 `424/2110` (`20.09%`). 以复合物为簇的配对差值为 `+7.58` 个百分点, 95% bootstrap CI 为 `[+2.61,+12.70]`, coordinate-only 显著更高.
- 官方生成质量: coordinate-only 为 `606/2110` (`28.72%`), official inpainting 为 `1209/2110` (`57.30%`). 配对差值为 `-28.58` 个百分点, 95% CI 为 `[-35.40,-21.61]`, official inpainting 显著更高.
- Official inpainting 更擅长生成几何合理的局部新结构, 但只有 `691/2110` 保持完整 isomeric identity. Coordinate-only 保持全部非立体化学图, 但较大的可编辑区仍容易产生严重内部几何和 strain 问题.

因此, 两种方法回答的是不同问题. 若目标是严格修回同一个分子, 当前 coordinate-only 更合适; 若允许可编辑区重新设计, official inpainting 的官方生成质量明显更高. 两者都还不足以直接作为最终通用修复器.

## 2. 构造漏斗

| 阶段 | 复合物数 |
|---|---:|
| SPINDR 官方 test 条目 | 225 |
| G_good 基线 PB 失败 | 2 |
| 无可旋转键 | 8 |
| 无中等强度 clash 候选 | 3 |
| 候选未保持非蛋白 PB 项 | 1 |
| 最终 eligible | 211 |

211 个 `G_bad` 均只通过较小侧支链扭转构造, 固定原子坐标不变, 并满足:

- 新增 protein-ligand clash 为 1-4 对, 分布为 `83、64、39、25` 个案例.
- 最小相对 VDW 距离范围为 `0.5344-0.7448`, 中位数 `0.6273`.
- 可编辑原子数范围为 `1-18`, 中位数为 `6`.
- 191 个案例使用排序后的首个候选, 16 个使用第二个, 4 个使用第三个.
- 211 个 `G_good` 均通过完整 dock PoseBusters, 211 个 `G_bad` 均因人工 protein clash 不通过.

## 3. 运行完整性

- 两种方法均使用 `flowr_root_v2.2`, 100 个 flow 步, batch size 1, 相同 oracle mask 和相同 case/seed 清单.
- Coordinate-only 和 official inpainting 均为 `2110/2110 completed`, 无 OOM、模型错误、补抽或结果筛选.
- 4,220 个方法键全部唯一, 2,110 个 case/seed 均恰好包含两种方法, 10 个 rollout 索引各出现 422 次.
- Coordinate-only 推理中位耗时 `2.90 s`, 总计 `1.70 h`; official inpainting 中位耗时 `3.42 s`, 总计 `2.01 h`.
- 5-step 与 100-step GPU smoke 均通过. 正式评价完整执行两次, 两次主要计数和统计结论一致.
- 输出清单包含 4,231 个文件条目, 全部大小和 SHA-256 复核通过.

## 4. 总体结果

| 指标 | Coordinate-only | Official inpainting |
|---|---:|---:|
| 尝试/完成 | 2110/2110 | 2110/2110 |
| RDKit valid | 2110 | 2094 |
| Fully-connected valid | 2110 | 2093 |
| 官方条件匹配 | 2110 | 2087 |
| 完整 PB-valid | 606 | 1209 |
| 官方质量成功 | 606 | 1209 |
| Isomeric identity 相同 | 1736 | 691 |
| 固定原子身份保持 | 2110 | 1930 |
| 固定片段内部键保持 | 2110 | 2088 |
| 固定坐标保持 | 2110 | 2110 |
| Strict success | 584 | 424 |

相同 case/seed 的配对结果为:

| 终点 | 两者通过 | 仅 coordinate 通过 | 仅 official 通过 | 两者失败 |
|---|---:|---:|---:|---:|
| Strict success | 233 | 351 | 191 | 1335 |
| 官方质量成功 | 352 | 254 | 857 | 647 |

## 5. 失败漏斗

| 累计阶段 | Coordinate-only 剩余 | Official inpainting 剩余 |
|---|---:|---:|
| Attempted | 2110 | 2110 |
| Model completed | 2110 | 2110 |
| Valid | 2110 | 2094 |
| Fully-connected valid | 2110 | 2093 |
| Condition match | 2110 | 2087 |
| PB-valid | 606 | 1209 |
| 再要求 same molecule | 584 | 424 |
| 再要求固定原子、键和坐标 | 584 | 424 |

Coordinate-only 的主要 PB 失败项是 internal energy `1351`, bond angles `1347`, bond lengths `1295` 和 internal steric clash `1212`. Official inpainting 的主要失败项是 bond lengths `782` 和 bond angles `392`, 但其 internal energy 失败只有 `114`, internal steric clash 失败只有 `38`.

Strain 进一步显示两者的几何差异:

| Strain, kcal/mol | Coordinate-only | Official inpainting | G_good | G_bad |
|---|---:|---:|---:|---:|
| 可计算 | 2090 | 2066 | 209 | 209 |
| 中位数 | 15924.66 | 119.36 | 31.86 | 48.27 |
| Q1-Q3 | 93.68-123326.25 | 65.18-250.68 | 20.19-52.95 | 28.99-83.85 |

## 6. 可编辑区大小

| 可编辑原子数 | Cases | Coordinate strict | Official strict | Coordinate 官方质量 | Official 官方质量 |
|---|---:|---:|---:|---:|---:|
| 1-3 | 54 | 376/540 | 207/540 | 382/540 | 305/540 |
| 4-6 | 65 | 193/650 | 131/650 | 207/650 | 349/650 |
| 7-12 | 84 | 15/840 | 84/840 | 17/840 | 493/840 |
| 13+ | 8 | 0/80 | 2/80 | 0/80 | 62/80 |

Coordinate-only 的优势主要来自 1-6 个可编辑原子的局部修复. 当可编辑区达到 7 个以上时, 它的内部几何明显崩溃. Official inpainting 在较大可编辑区的官方质量优势很强, 但 strict 成功仍受到换分子和固定边界改变限制.

## 7. 分子身份解释

Coordinate-only 没有重新采样原子、键或电荷. 其 2,110 个结果的 canonical non-isomeric identity 全部与 `G_good` 相同, 但有 374 个结果因三维坐标变化导致立体化学身份不同. 其中 22 个原本已通过官方质量, 因此 strict 从 606 降为 584.

Official inpainting 中:

- 691 个结果保持完整 isomeric identity.
- 43 个保持非立体化学身份, 但立体身份变化.
- 1360 个连非立体化学身份也变化.
- 16 个无法完成 RDKit sanitization.

这与方法本身一致: official inpainting 会重新采样可编辑区的离散化学结构, 所以 PB-valid 不能自动解释为修回原来的分子.

## 8. 边界与下一步

本实验只覆盖可自动构造的中等强度扭转 clash, 使用 oracle mask, 不评价 Fixed Mask Head 或自然错误分布. 它也不是 FLOWR.ROOT 官方 benchmark 复现.

下一步应优先改进 strict coordinate repair 的几何保持能力, 尤其是 7 个以上可编辑原子、环系和边界共价几何. 在 oracle strict 成功率明显提高前, 不建议进入 Fixed Mask Head 的正式训练.

逐样本结果见 [runs.csv](../../outputs/20260723-04-full-test-oracle-repair-benchmark/runs.csv), 逐 case 结果见 [case_rates.csv](../../outputs/20260723-04-full-test-oracle-repair-benchmark/case_rates.csv), 基线见 [baselines.csv](../../outputs/20260723-04-full-test-oracle-repair-benchmark/baselines.csv), 完整统计见 [summary.json](../../outputs/20260723-04-full-test-oracle-repair-benchmark/summary.json), 文件哈希见 [output_manifest.json](../../outputs/20260723-04-full-test-oracle-repair-benchmark/output_manifest.json).
