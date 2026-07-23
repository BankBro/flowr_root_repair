# FLOWR.ROOT 官方核心指标复评计划

**实验 ID:** `20260723-03-official-metrics-reevaluation`

**状态:** 已执行

## 1. 目的

直接复评前两次实验已经保存的两组各 50 个 SDF, 不重新采样模型. 本实验回答官方生成质量口径下两组产物分别表现如何, 不取代已经冻结的 native/strict 结论, 也不重新输出 `GO/NO-GO`.

## 2. 输入与边界

- Coordinate-only 组来自 `20260723-01-oracle-mask-repair-pilot` 的 50 个正式 SDF.
- Official inpainting 组来自 `20260723-02-official-fragment-inpainting-control` 的 50 个正式 SDF.
- 两组必须覆盖完全相同的 5 个 case 和 50 个 seed. 复评前记录全部输入哈希.
- 使用相同的 SPINDR prepared protein PDB、`G_bad` 和冻结 editable atom indices.
- 全程只使用 CPU, 不加载 checkpoint, 不重新运行 flow sampling.

## 3. 官方核心指标

每个产物计算:

1. FLOWR.ROOT `mol_is_valid` 的 validity 和 fully-connected validity.
2. `PoseBusters(config="dock")` 的全部官方检查列; 只有全部列通过才记为 PB-valid.
3. FLOWR.ROOT `evaluate_strain`, 使用加氢分子、MMFF94s 和默认 500 步优化.
4. FLOWR.ROOT `check_substructure_match`, 使用 `substructure_inpainting` 和冻结 editable atom indices.

条件子结构检查沿用官方语义: 只要求从参考配体得到的固定子结构仍能在输出中匹配, 不要求固定原子逐一映射、坐标相同或完整分子身份相同. 无效分子不计算条件匹配和 strain, 但仍保留在 50 次固定分母中.

## 4. 汇总与解释

- 报告总体和逐 case 的 validity、fully-connected validity、condition match、PB-valid 计数与比例.
- 条件匹配同时报告相对全部 50 次和相对 fully-connected valid 产物的分母.
- Strain 报告可计算数量、均值、标准差、中位数和四分位距, 不设置人为通过阈值.
- 按相同 case/seed 报告 PB-valid 与 conditioned PB-valid 的双方通过、单方通过和双方失败计数.
- 使用相同官方指标评估 5 个 `G_good/G_bad` 作为基线.
- `conditioned_pb_valid = fully_connected_valid and condition_match and pb_valid` 只表示官方工作流兼容产物率, 不命名为修复成功率.

## 5. 执行与产物

1. 单元测试 validity、条件匹配、PB 全列聚合、缺失 strain 和配对汇总.
2. Preflight 校验 100 个结果、5 个 `G_good/G_bad` 和 5 个蛋白质 PDB.
3. 对同一 case 的两组各一个样本运行 smoke.
4. 正式复评全部 100 个结果, 生成逐样本 CSV、基线 CSV、输入清单和汇总 JSON.
5. 创建同 slug 的 report, 更新 `docs/STATUS.md` 和 `docs/EXPERIMENT_LOG.md`.

主要入口放在 `experiments/20260723-03-official-metrics-reevaluation/run.py`, 结果放在 `outputs/20260723-03-official-metrics-reevaluation/`.
