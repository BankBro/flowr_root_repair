# Official Fragment Inpainting Control 实验报告

**实验 ID:** `20260723-02-official-fragment-inpainting-control`

**对应计划:** [20260723-02-official-fragment-inpainting-control-plan.md](../plan/20260723-02-official-fragment-inpainting-control-plan.md)

**Strict 结论:** `NO-GO`

## 1. 目的与结论

本实验在 `20260723-01` 完全相同的 5 个 `G_bad`、真值 `M_fixed` 和 50 个 seed 上, 使用 FLOWR.ROOT 官方 fragment inpainting 完整生成路径做公平对照.

共同评估器从两组最终 SDF 读取结果后得到:

- Coordinate-only: native `14/50`, strict `13/50`.
- Official fragment inpainting: native `14/50`, strict `7/50`.
- Native 总数相同, 差值为 `0` 个百分点. Strict 官方组少 6 次, 差值为 `-12` 个百分点.

官方组 strict `7/50`, 低于冻结 NO-GO 上限 `9/50`, 因此不能作为严格同分子坐标修复器的直接替代. 它在部分案例上能生成没有 clash 且内部几何合格的新分子, 更符合局部重设计方法的定位.

## 2. 执行完整性

- Plan commit: `3d5f81d`.
- 初始实现 commit: `34f37fc`.
- 正式采样实现 commit: `52e114f6fcf48f2797c55279efd2391ea3b5a424`.
- SDF 共同评估 commit: `5484f985c4f9fa4d44e5318328058a2dd39899c4`.
- Checkpoint SHA256: `b818f41dc12ffb6bc558bb0ad997055581e07cd9e49dcac1b794ed9993c46e4c`.
- 22 项测试通过, 包括 repair 单元测试、mask/index 映射、双终点、漏斗和官方生成取消回归.
- Preflight 通过. 冻结案例、旧 CSV、旧 summary 和 5 个 `G_bad` 哈希均一致; 共同评估器复现旧组 native `14/50`、strict `13/50`.
- 5 步和 100 步 GPU smoke 通过. 相同 seed 重跑的完整输出张量 SHA256 一致.
- 正式运行 50/50 completed, batch size 1, 无 OOM、模型异常、补抽或筛选. 50 个输出均保存为 SDF, 包括 1 个可构建但无法 sanitise 的失败分子.

正式采样使用官方 Gaussian 坐标 prior、均匀离散 prior、100 步 linear Euler、corrector 0 和 `final_inpaint=true`. 唯一方法适配是将冻结的真值 `M_fixed` 注入官方 `_build_inference_prior()`. 未修改 `flowr/gen/` 或 `flowr/models/` 的官方实现.

## 3. 公平对照结果

Native 允许可编辑区生成不同分子. Strict 还要求完整分子的 canonical isomeric SMILES 与 `G_bad` 相同.

| Case | Coordinate native | Coordinate strict | Official native | Official strict |
|---|---:|---:|---:|---:|
| 6PVZ | 10/10 | 10/10 | 6/10 | 6/10 |
| 7DDL | 0/10 | 0/10 | 0/10 | 0/10 |
| 3ROG | 3/10 | 3/10 | 1/10 | 1/10 |
| 4BV5 | 0/10 | 0/10 | 6/10 | 0/10 |
| 4F0S | 1/10 | 0/10 | 1/10 | 0/10 |
| **合计** | **14/50** | **13/50** | **14/50** | **7/50** |

相同 seed 的成对结果也表明两种方法不是简单的强弱替代:

| Endpoint | 两者成功 | 仅 coordinate | 仅 official | 两者失败 |
|---|---:|---:|---:|---:|
| Native | 6 | 8 | 8 | 28 |
| Strict | 6 | 7 | 1 | 36 |

案例数只有 5 个, 本实验不进行显著性检验, 也不将计数差异外推为数据集总体性能.

## 4. 逐层漏斗

每一层只保留同时通过此前所有条件的结果:

| 筛选层级 | Coordinate 剩余 | Official 剩余 | Official 本层淘汰 |
|---|---:|---:|---:|
| 固定尝试数 | 50 | 50 | 0 |
| 正常完成、有限、可 sanitise、单分子且固定片段保持 | 50 | 40 | 10 |
| 无 protein clash | 40 | 37 | 3 |
| 内部键长、键角和 internal clash 全部通过, 即 native | 14 | 14 | 23 |
| 完整分子身份和立体化学相同, 即 strict | 13 | 7 | 7 |

官方组 50/50 坐标有限, 49/50 可 sanitise 且为单连通分子. 50 次固定原子写文件前的原始最大漂移均小于 `1e-5 A`, 范围为 `7.53e-7` 至 `1.97e-6 A`, 均按预注册规则精确恢复后保存.

`completed_usable_fixed` 的 10 次淘汰来自:

- 7DDL 有 9/10 改变了固定边界原子 21 的芳香性. 坐标和固定区域内部键未漂移, 但生成区域与边界重新连接后, 固定原子的化学身份不再相同, 因而不能算固定片段完整保留.
- 4F0S 有 1/10 产生 N 原子显式价态不合法的分子, 无法通过 RDKit sanitise.

## 5. 方法差异的含义

官方组 47/50 输出没有 protein clash, 但仅 18/50 独立通过内部几何, 经过固定片段和前置条件逐层筛选后剩 14/50 native success. 它能重新生成原子类型、电荷和键, 因而在 4BV5 上得到 6 个合格的新局部结构, 这是 coordinate-only 没有做到的.

代价是分子身份更难保持. 官方组只有 23/50 输出与原分子 canonical isomeric SMILES 相同, 最终只有 7 个同时满足几何、clash、固定片段和同分子要求. 14 个 native success 中有 7 个其实是不同分子, 不能称为把原配体的坐标修好.

Coordinate-only 不改变模型张量中的离散图, 因而更贴近本课题的同分子坐标修复语义; 它的主要瓶颈仍是可编辑区域共价几何. 官方完整 inpainting 更像“保留一部分配体, 重新设计其余部分”, 不能用来证明当前 repair 路线已经解决几何问题.

## 6. SDF 评估审计

官方输出写文件前为 strict `6/50`, 从 SDF 重读后为 `7/50`. 变化来自 RDKit 根据 SDF 的 3D 构型恢复手性标记, 其中 3ROG seed `2026072323` 从不同分子改判为相同分子并通过 strict. Native 计数始终为 `14/50`.

最终结果以重读 SDF 的共同评估为准, 因为 coordinate-only 对照也从保存的 SDF 重评. 13 个涉及 `same_molecule` 的写前/读后差异完整记录在 `artifact_audit.json`, 没有隐藏或回写采样输出.

## 7. 决策与限制

本实验的 strict 结论为 `NO-GO`. 不采用官方 full fragment inpainting 作为首版严格坐标修复器, 也不据此开始 Fixed Mask Head 训练. 下一步仍应在 oracle mask 条件下研究能保持同一分子及局部共价几何的坐标表示、约束或采样方式.

主要限制包括: 只有 5 个人工扭转案例; 使用真值 mask; 未进行过滤、重采样或后处理优化; 官方组为 batch size 1, 旧组为 batch size 2; native 的固定片段要求会将固定边界芳香性改变视为失败. 这些都是本次受控实验的既定边界, 结论不代表 FLOWR.ROOT 在一般局部生成任务上的总体能力.

逐次结果见 [official_runs.csv](../../outputs/20260723-02-official-fragment-inpainting-control/official_runs.csv), 对照与漏斗见 [comparison.json](../../outputs/20260723-02-official-fragment-inpainting-control/comparison.json), 完整性审计见 [artifact_audit.json](../../outputs/20260723-02-official-fragment-inpainting-control/artifact_audit.json).
