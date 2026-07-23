# SPINDR 全测试集 Oracle Repair Benchmark 计划

**实验 ID:** `20260723-04-full-test-oracle-repair-benchmark`

**状态:** 已完成

## 1. 目的

对 SPINDR 官方 test split 的 225 个蛋白质-配体复合物做完整漏斗审计, 自动构造可复现的中等强度局部 clash, 再使用完全相同的 `G_bad`、oracle mask 和 seed, 比较冻结的 coordinate-only 与 official inpainting.

本实验区分两个问题:

1. 官方生成质量: 产物是否通过 FLOWR.ROOT 官方 validity、条件匹配和完整 dock PoseBusters.
2. 严格修复: 产物是否在官方质量合格的同时保持同一分子、固定片段和固定坐标.

## 2. 自动 Clash 构造

- 统一使用与 LMDB 对齐的 225 组预处理 pocket CIF/SDF, 将 CIF 链名确定性映射为单字符链名后交给 RDKit/PoseBusters.
- `G_good` 必须可读取、与 LMDB 原子顺序一致、通过完整 dock PoseBusters, 并至少包含一个 RDKit 可旋转键.
- 对每个可旋转键选择断键后的较小一侧支链, 枚举 `±30、±60、±90、±120、±150、180` 度扭转.
- 合格 `G_bad` 必须保持化学图和固定坐标, 产生 1-4 个 protein-ligand clash, 最小相对 VDW 距离位于 `[0.50, 0.75)`, 并保持所有非蛋白 PoseBusters 项通过.
- 候选依次按距 `0.625` 的偏差、角度绝对值、editable 原子数和稳定索引排序, 每个复合物只保留一个确定性 `G_bad`.
- 全部 225 个条目都记录漏斗状态. 无合格候选时保留失败原因, 不人工替换或放宽标准.

## 3. 正式运行

- 两种方法均使用 `flowr_root_v2.2` checkpoint、100 个 flow 步和 batch size 1, 不修改当前 prior、约束或官方生成流程.
- 每个 eligible case 每种方法运行 10 个固定 seed. `seed = 2026072300 + 10 × test_index + rollout_index`.
- 每个 `case × seed × method` 只采样一次. 无效分子、PB-fail 和模型错误均保留在固定分母中; 仅允许同 key 断点续跑未完成任务.
- 正式运行前依次完成 CPU 测试与 construction preflight、5-case/5-step GPU smoke 和一个中位 case 的 100-step smoke.
- 输入、checkpoint、配置和结果均记录 SHA-256, 两种方法必须覆盖完全相同的 case/seed 清单.

## 4. 评价与统计

官方质量终点定义为:

```text
fully_connected_valid
and official_condition_match
and complete_dock_pb_valid
```

主要终点 strict repair 在官方质量终点基础上进一步要求:

```text
canonical_isomeric_identity 与 G_good 相同
and 固定原子与固定键保持
and 映射后固定坐标漂移 <= 1e-5 Å
```

同时报告 validity、每个 PoseBusters 子项、MMFF94s strain、editable/all-atom RMSD、失败漏斗和逐 case 的 10-seed 成功率. 不加入 Vina、唯一性或多样性.

主要比较为 coordinate-only 减 official inpainting 的 strict case-level 成功率差. 使用复合物为簇的 10,000 次 paired bootstrap, 固定统计 seed `2026072304`. 95% CI 不跨 0 时才判断某方法更优, 跨 0 则记为证据不足. 官方质量终点作为次要比较.

## 5. 边界与产物

- 本实验只评价 oracle mask 下、可自动构造的中等强度扭转 clash, 不评价 Fixed Mask Head.
- 结果不代表自然 clash 分布, 也不声称复现 FLOWR.ROOT 官方 benchmark.
- plan、experiment、output 和 report 使用相同实验 ID 一一对应.
- 实验代码放在 `experiments/20260723-04-full-test-oracle-repair-benchmark/`, 输出放在同名 `outputs/` 目录.
