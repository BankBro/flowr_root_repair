# Oracle-mask Repair Pilot 实验报告

**实验 ID:** `20260723-01-oracle-mask-repair-pilot`

**对应计划:** [20260723-01-oracle-mask-repair-pilot-plan.md](../plan/20260723-01-oracle-mask-repair-pilot-plan.md)

**结论:** `CONDITIONAL`

## 1. 目的与结论

本实验使用真实 `M_fixed` 和冻结的 FLOWR.ROOT v2.2 checkpoint, 测试 coordinate-only 局部修复是否具备基本可行性. 5 个 SPINDR validation 案例各运行 10 个固定 seed, 共 50 次.

50 次均正常完成, 14/50 同时满足无 protein clash、内部几何合格、固定坐标不动和离散图不变, 成功率为 28%. 按预先冻结的 gate, 结果为 `CONDITIONAL`, 未达到进入 Fixed Mask Head 训练的 `GO` 标准.

## 2. 执行完整性

- 实现 commit: `bb75146025149443824edf59b48847949aaa484b`.
- checkpoint SHA256: `b818f41dc12ffb6bc558bb0ad997055581e07cd9e49dcac1b794ed9993c46e4c`.
- CPU preflight 通过: 5 个 case 均属于 68 个对象的官方 validation split, raw/LMDB 原子顺序一致, 人工扭转和 PoseBusters 基线符合冻结定义.
- 测试通过: 13 个 repair 单元测试和 3 个官方生成路径回归测试, 共 16 个.
- GPU smoke 通过: 6PVZ 的 5 步和 100 步 smoke 均无崩溃、NaN 或约束破坏.
- 正式运行: 100 步 linear Euler, batch size 2, 50/50 completed, 无 OOM、重试或补抽 seed.
- 硬约束: 50/50 固定原子写文件前漂移为 `0.0 A`, 50/50 冻结的离散图字段完全不变.

首次正式执行后发现, 中心坐标还原到物理坐标时仍保留 float32 量化误差, 导致部分固定原子被错误计为漂移. 该轮未作为正式结果. 在 commit `bb75146` 中加入物理坐标系的最终精确覆盖后, 使用相同 50 个 seed 完整重跑. 本报告只采用重跑结果, 模型预测、prior、积分设置和 gate 均未调整.

## 3. 正式结果

| Case | `G_bad` clash | 无 protein clash | 内部几何通过 | 固定/图保持 | 最终成功 |
|---|---:|---:|---:|---:|---:|
| 6PVZ | 1 | 10/10 | 10/10 | 10/10 | 10/10 |
| 7DDL | 1 | 10/10 | 0/10 | 10/10 | 0/10 |
| 3ROG | 2 | 3/10 | 9/10 | 10/10 | 3/10 |
| 4BV5 | 3 | 10/10 | 0/10 | 10/10 | 0/10 |
| 4F0S | 4 | 7/10 | 1/10 | 10/10 | 1/10 |
| **合计** | - | **40/50** | **20/50** | **50/50** | **14/50** |

实验级 gate 要求 GO 至少 25/50 成功, 且至少 4 个 case 各达到 3/10. 本实验只有 14/50 成功, 且只有 6PVZ 和 3ROG 达到逐 case 下限, 因此不能判定为 GO. 14/50 又高于 NO-GO 的 9/50 上限, 所以按协议判定为 `CONDITIONAL`.

诊断指标:

- editable RMSD: 中位数 `2.510 A`, 均值 `2.266 A`, 范围 `0.159-3.490 A`.
- 全配体 RMSD: 中位数 `1.396 A`, 均值 `1.184 A`, 范围 `0.062-2.150 A`.
- 单样本运行时间: 中位数 `1.476 s`, 均值 `1.490 s`.
- batch 峰值 GPU 显存: 中位数 `1248 MiB`, 最大 `1754 MiB`.

## 4. 失败模式与解释

- 40/50 已消除 protein clash, 说明真值掩码配合局部 FLOWR.ROOT 采样能够把可编辑区域移出蛋白质冲突.
- 30/50 未通过内部几何. 7DDL 和 4BV5 的 clash 全部消失, 但所有输出都破坏了键长或键角; 4F0S 也只有 1/10 保持完整内部几何.
- 3ROG 有 9/10 保持内部几何, 但仅 3/10 消除 protein clash, 表明不同案例的主要困难并不相同.
- 6PVZ 的 3 原子小区域达到 10/10, 是局部 coordinate-only 修复可行的正向证据, 但不能外推到更大或更复杂的可编辑区域.
- 固定坐标和冻结的离散图字段均为 50/50, 说明 `M_fixed` 方向、逐步硬覆盖和 coordinate-only 图约束实现可靠, 不是本轮主要瓶颈.

Post-hoc 产物检查发现, 50 个 SDF 重新读取后, 计划冻结的原子身份、电荷、杂化、芳香性和键字段仍为 50/50 一致, 但包含原子手性和键立体标记的严格签名有 18/50 发生变化, 其中 4BV5 为 8/10, 4F0S 为 10/10. 立体化学未包含在本轮预注册成功标准中, 因此不回改 gate, 但后续实验必须将其作为独立约束或评估项.

当前证据指向: 自由笛卡尔坐标 prior 和 Structure Head 输出缺少显式共价几何约束, 是扩大修复区域时的主要限制. 这是对 5 个受控案例的初步证据, 不是数据集总体修复率或正式 benchmark 结论.

## 5. 决策与下一步

本实验停止在 `CONDITIONAL`, 不在同一实验内调参, 也不开始 Fixed Mask Head 训练.

下一步应新建独立 plan, 继续使用 oracle mask, 优先验证能够保持局部共价几何的修复表示或约束. 在 oracle 修复达到预先定义的 GO 标准后, 再进入掩码预测训练.

完整逐次结果见 [runs.csv](../../outputs/20260723-01-oracle-mask-repair-pilot/runs.csv), 汇总见 [summary.json](../../outputs/20260723-01-oracle-mask-repair-pilot/summary.json), 50 个最终结构见 [repairs](../../outputs/20260723-01-oracle-mask-repair-pilot/repairs/).
