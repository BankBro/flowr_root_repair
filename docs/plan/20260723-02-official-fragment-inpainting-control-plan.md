# Official Fragment Inpainting Control 实验计划

**实验 ID:** `20260723-02-official-fragment-inpainting-control`

**对照实验:** [20260723-01-oracle-mask-repair-pilot-plan.md](20260723-01-oracle-mask-repair-pilot-plan.md)

**状态:** 已执行, strict 结论为 `NO-GO`

## 1. 目标

在 `20260723-01` 完全相同的 5 个人工 clash 案例、真值 `M_fixed` 和 50 个正式 seed 上, 使用 FLOWR.ROOT 原始 fragment inpainting 生成流程进行修复. 本实验回答两个问题:

1. 官方完整生成方式能否得到化学有效、保留固定片段且消除 clash 的局部重设计结果.
2. 官方完整生成方式能否恢复为与输入相同的分子, 从而作为严格意义上的坐标修复器.

本实验不训练模型, 不修改官方生成算法, 也不把 5 个案例外推为正式 benchmark.

## 2. 公平对照边界

- 直接复用 `20260723-01` 冻结的 `cases.json`、`G_bad/G_good`、蛋白质 PDB 和 50 个 seed, 不重新制造缺陷.
- 坐标修复组直接复用其已冻结输出, 不重新采样.
- 两组均由本实验的共同评估器重新评价. 正式运行前必须复现坐标修复组的 native `14/50` 和 strict `13/50`.
- 官方组每个 seed 只生成一次, 不筛选、不补抽、不重采样. 无效输出计为失败并保留在分母中.
- 官方组采用 batch size 1, 且每次重置 Python、NumPy、PyTorch CPU/CUDA 随机数状态. 上一轮坐标修复组实际使用 batch size 2, 该差异写入报告.

## 3. 官方生成协议

- checkpoint 固定为 FLOWR.ROOT v2.2, SHA256 为 `b818f41dc12ffb6bc558bb0ad997055581e07cd9e49dcac1b794ed9993c46e4c`.
- 使用 `ComplexInterpolant._build_inference_prior(..., mode="fragment_inpainting", mask=M_fixed, is_local=True)` 构建官方先验.
- 唯一方法适配是注入冻结的真值 `M_fixed`, 替代官方推理时随机选择 fragment. `M_fixed=1` 表示固定原子.
- 使用官方 Gaussian 坐标 prior、均匀原子和键离散 prior, 连续和离散时间均从 0 开始.
- 使用原始 `model._generate()`、100 步 linear Euler、corrector 0, 不使用 SDE、guidance、优化或质子化.
- `sample_mol_sizes=false`, 保持与输入相同的重原子数.
- 显式设置 `final_inpaint=true`, 保证最终输出也恢复固定片段.
- fragment prior 的输出顺序是原始固定原子后接原始可编辑原子. 实验代码必须保存该索引映射, 评估时还原到原始原子语义.

固定原子写入物理坐标时可能出现 float32 量化误差. 先记录未修正最大漂移. 当其不超过 `1e-5 A` 时, 允许在保存和共同评估前精确恢复为 `G_bad` 的 float64 固定坐标; 超过阈值时固定片段判为失败, 不得用恢复操作掩盖模型错误.

## 4. 两类主要终点

### 4.1 Native Local Redesign Success

单次结果必须同时满足:

1. RDKit 可构建并 sanitise, 且只有一个连通分子.
2. 真值固定片段的原子身份、键关系和坐标得到保留.
3. PoseBusters 蛋白质-配体 pairwise clash 为 0.
4. PoseBusters 默认内部键长、键角和 internal clash 检查全部通过.

该指标允许可编辑区域生成不同但有效的原子、键和立体化学, 因此表示局部重设计能力.

### 4.2 Strict Same-Molecule Repair Success

除满足 native 条件外, 完整输出还必须与 `G_bad` 保持相同的分子图、形式电荷、芳香性及 canonical isomeric SMILES, 包括原子手性和键立体信息. 该指标表示严格的同分子坐标修复能力, 且不得依赖输出原子顺序相同.

RMSD 仅作诊断, 不设成功阈值.

## 5. Gate 与比较方法

Strict endpoint 沿用上一轮实验级 gate:

- **GO:** 至少 25/50 成功, 且至少 4 个案例各达到 3/10.
- **NO-GO:** 最多 9/50 成功.
- **CONDITIONAL:** 其余情况.

Native endpoint 作为并列主要描述结果, 但不单独决定能否进入同分子修复路线. 如果 native 表现良好而 strict 未通过, 结论只能是官方方法适合局部重设计, 不能称为成功修复同一个分子.

报告给出两组总体和逐案例计数、百分点差值以及逐层漏斗. 由于只有 5 个案例, 不进行显著性检验, 也不以一次额外成功宣称方法优越.

## 6. 执行流程

1. CPU 测试索引映射、固定片段检查、native/strict 口径和漏斗汇总.
2. Preflight 校验 checkpoint、案例和旧结果哈希, 确认源文件存在, 并用共同评估器复现旧组 native `14/50`、strict `13/50`.
3. 运行现有 repair 测试与官方生成路径回归测试.
4. 对 6PVZ 使用非正式 seed `424242/424243`, 分别运行 5 步和 100 步 GPU smoke. 重复一个 seed 验证确定性.
5. 按 `2026072300 + 10 * case_index + rollout_index` 执行 50 次正式采样, 支持按 case 和 seed 断点续跑.
6. 共同评估两组结果, 生成逐次表格、汇总、比较漏斗和报告.

若正式 seed 出现模型异常、NaN 或无效分子, 记为失败. 若 batch size 1 下重复发生 OOM 或基础设施错误而无法继续, 实验标记为未完成, 不输出 gate 结论.

## 7. 代码与产物

- 可复用的共同评估逻辑放在 `flowr/repair/`.
- 真值 mask 适配器、实验入口和测试辅助代码放在 `experiments/20260723-02-official-fragment-inpainting-control/`.
- 不修改 `flowr/gen/` 和 `flowr/models/` 中的官方默认实现.
- 正式输出保存在同 ID 的 output 目录, 包括 `official_runs.csv`、`summary.json`、比较漏斗以及可构建的逐样本 SDF.
- 输出分子评估不得读取 `G_good`; `G_good` 仅用于输入追溯和诊断 RMSD.
- 完成、失败或中止均创建对应 report, 并更新 `STATUS.md` 与 `EXPERIMENT_LOG.md`.
