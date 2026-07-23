# Oracle-mask Repair Pilot 实验计划

**实验 ID:** `20260723-01-oracle-mask-repair-pilot`

**上位设计:** [20260722-01-flowr-root-ligand-repair-v1-plan.md](20260722-01-flowr-root-ligand-repair-v1-plan.md)

**状态:** 已执行, 结论为 `CONDITIONAL`

## 1. 目标

使用冻结的 FLOWR.ROOT v2.2 checkpoint, 在 5 个验证集复合物上各运行 10 次真值掩码局部修复, 共 50 次. 本实验只判断 coordinate-only 修复器是否具备基本可行性, 不训练 Fixed Mask Head, 也不将结果外推为数据集总体上限或正式 benchmark.

## 2. 冻结案例

原子索引基于 checkpoint 使用的无氢 LMDB 表示. 角度遵循围绕有序轴的右手旋转, 旋转轴原子保持固定.

| Case | 有序旋转轴 | 角度 | `M_fixed=0` 重原子 | 预期新增 clash |
|---|---:|---:|---|---:|
| 6PVZ | 6 -> 5 | +30 deg | 4, 17, 19 | 1 |
| 7DDL | 20 -> 21 | -120 deg | 22-27 | 1 |
| 3ROG | 8 -> 9 | -120 deg | 0-3 | 2 |
| 4BV5 | 18 -> 16 | -90 deg | 7, 8, 13, 14, 15, 17, 19 | 3 |
| 4F0S | 0 -> 14 | -60 deg | 10-13, 15-18 | 4 |

完整机器可读定义保存在实验目录的 `cases.json` 中.

## 3. 修复协议

- `M_fixed=1` 表示固定, `M_fixed=0` 表示修复. 模型输入为 `fragment_mask = ligand_mask & M_fixed`, padding 为 0.
- 对可编辑原子采样 `z ~ N(0, I)`, 再整体平移, 使 prior 质心等于 `G_bad` 可编辑区域质心. prior 不读取 `G_good`.
- 固定坐标复制 `G_bad`. 原子类型、电荷、杂化和全部键始终复制输入化学图.
- 采用 100 步 linear Euler, corrector 为 0, 不使用 SDE 或额外坐标噪声, 保留 checkpoint 原有 self-conditioning.
- 连续时间中固定原子为 1, 可编辑原子为当前 `t`; 所有离散时间恒为 1.
- 每个积分步、self-conditioning 状态和最终预测后都恢复固定坐标与完整离散图.
- 正式 seed 为 `2026072300 + 10 * case_index + rollout_index`, case 和 rollout 均从 0 开始.
- 默认 batch size 为 2. OOM 只允许用相同 seed 降到 batch size 1 重试一次.

## 4. 成功标准

单次成功必须同时满足:

1. PoseBusters 蛋白质-配体 pairwise clash 为 0, 使用 VDW 半径、`radius_scale=1.0`、`clash_cutoff=0.75`.
2. PoseBusters 默认内部键长、键角和 internal clash 检查全部通过.
3. 写文件前的固定原子最大坐标漂移不超过 `1e-6 A`.
4. 原子身份、形式电荷、杂化、芳香性和键邻接/键类型与输入完全一致.

RMSD 只作诊断. 在口袋坐标系中直接报告 editable RMSD 和全配体 RMSD, 不做刚体对齐, 不设置 RMSD 成功阈值.

实验级 gate:

- **GO:** 至少 25/50 成功, 且至少 4 个样本各自达到 3/10.
- **NO-GO:** 最多 9/50 成功.
- **CONDITIONAL:** 其余情况.

CONDITIONAL 或 NO-GO 不在本实验内调参, 也不进入 Fixed Mask Head 训练. 后续调整必须新建 plan.

## 5. 执行流程

1. CPU preflight 检查 checkpoint 哈希、LMDB/raw ID 与原子顺序、冻结 corruption、拓扑、固定区域和 PoseBusters 基线.
2. 运行 repair 单元测试与原始生成路径回归测试.
3. 使用 6PVZ 和非正式 seed `424242/424243` 运行 5 步 GPU smoke 与 100 步 GPU 确认. Smoke 只要求无崩溃、无 NaN、约束保持和评估可运行, 不要求化学修复成功.
4. 冻结实现后运行 50 次正式采样. NaN、无效结构或模型输出异常计为失败, 不补抽 seed.
5. 汇总总体和逐样本成功率、四项条件失败分解、RMSD、耗时与显存, 再根据冻结 gate 得出结论.

如果基础设施错误在同 seed、batch size 1 下仍重复, 实验标记为未完成, 不输出 GO/NO-GO.

## 6. 代码与产物

- 可复用实现放在 `flowr/repair/`, 不修改官方 `_generate()` 默认行为.
- 实验入口、冻结案例和说明放在 `experiments/20260723-01-oracle-mask-repair-pilot/`.
- `G_good/G_bad`、全部 50 个最终 SDF、逐次 CSV 和汇总 JSON 放在同 ID 的 output 目录并提交 Git.
- 不提交 checkpoint、SPINDR 数据副本、采样轨迹、缓存和长期日志.
- checkpoint SHA256 固定为 `b818f41dc12ffb6bc558bb0ad997055581e07cd9e49dcac1b794ed9993c46e4c`.

## 7. 报告与交付

完成、失败或中止均创建同 slug 的 report. 报告记录实现 commit、完整性检查、基线、结果、gate、失败模式和局限. 计划、实现、结果与报告分批提交并推送到 `origin/main`, 且不纳入现有 `environment.yml` 与 `docs/paper/` 的无关改动.
