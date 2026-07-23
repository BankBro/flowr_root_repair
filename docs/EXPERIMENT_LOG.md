# 实验日志

本文件按时间顺序追加正式实验和会影响研究决策的阶段记录. 不重写历史条目; 需要纠正时追加更正说明.

正式实验记录见下文.

## 建议记录内容

每条记录至少说明日期, `experiment_id`, 对应 plan, 目的, 执行内容, 输出位置, 结果或失败原因以及下一步. 具体格式可以根据实验灵活调整.

## 2026-07-23: Oracle-mask repair pilot 开始

- `experiment_id`: `20260723-01-oracle-mask-repair-pilot`.
- plan: [20260723-01-oracle-mask-repair-pilot-plan.md](plan/20260723-01-oracle-mask-repair-pilot-plan.md).
- 目的: 在 5 个冻结的 SPINDR validation 人工 clash 案例上, 用 50 次真值掩码 rollout 验证 coordinate-only 修复器的基本可行性.
- 当前进度: 协议与案例已经冻结, experiment/output 骨架已建立, 即将实现 repair 内核并执行 preflight 和 smoke.
- 输出位置: `outputs/20260723-01-oracle-mask-repair-pilot/`.

## 2026-07-23: Oracle-mask repair pilot 完成

- `experiment_id`: `20260723-01-oracle-mask-repair-pilot`.
- 实现 commit: `bb75146025149443824edf59b48847949aaa484b`.
- 执行: CPU preflight、16 项测试、5/100 步 GPU smoke 和 5 case x 10 seed 正式采样均完成.
- 结果: 50/50 正常结束, 14/50 最终成功, 40/50 无 protein clash, 20/50 内部几何通过, 固定坐标和离散图均保持 50/50.
- gate: `CONDITIONAL`. 未达到 Fixed Mask Head 训练前要求的 oracle 修复 GO 标准.
- 主要限制: 大可编辑区域容易出现键长、键角或 internal clash 失败, 需要在新 plan 中验证几何保持方法.
- 报告: [20260723-01-oracle-mask-repair-pilot-report.md](report/20260723-01-oracle-mask-repair-pilot-report.md).
- 输出: `outputs/20260723-01-oracle-mask-repair-pilot/`.

## 2026-07-23: Official fragment inpainting control 开始

- `experiment_id`: `20260723-02-official-fragment-inpainting-control`.
- plan: [20260723-02-official-fragment-inpainting-control-plan.md](plan/20260723-02-official-fragment-inpainting-control-plan.md).
- 目的: 在完全相同的 5 个人工 clash 案例、真值 `M_fixed` 和 50 个 seed 上运行 FLOWR.ROOT 官方 fragment inpainting, 公平对照上一轮 coordinate-only 修复.
- 主要终点: 同时报告允许可编辑区换分子的 native local redesign success, 以及要求完整分子身份和立体化学不变的 strict same-molecule repair success.
- 当前进度: 协议已冻结, experiment/output 骨架已建立, 即将实现共同评估和官方 mask 适配器.
- 输出位置: `outputs/20260723-02-official-fragment-inpainting-control/`.

## 2026-07-23: Official fragment inpainting control 完成

- `experiment_id`: `20260723-02-official-fragment-inpainting-control`.
- 正式采样实现 commit: `52e114f6fcf48f2797c55279efd2391ea3b5a424`.
- 共同评估实现 commit: `5484f985c4f9fa4d44e5318328058a2dd39899c4`.
- 执行: 22 项测试、preflight、5/100 步 GPU smoke 和 5 case x 10 seed 正式采样均完成, 50/50 completed, 无 OOM、补抽或筛选.
- 结果: Coordinate-only 为 native `14/50`、strict `13/50`; 官方 fragment inpainting 为 native `14/50`、strict `7/50`.
- gate: 官方 strict 为 `NO-GO`. 官方方法可用于局部重设计, 但不能作为严格同分子坐标修复器的直接替代.
- 主要限制: 官方组 14 个 native success 中有 7 个是不同分子; 7DDL 还暴露了生成区重新连接后固定边界芳香性可能改变的问题.
- 报告: [20260723-02-official-fragment-inpainting-control-report.md](report/20260723-02-official-fragment-inpainting-control-report.md).
- 输出: `outputs/20260723-02-official-fragment-inpainting-control/`.

## 2026-07-23: Official metrics reevaluation 开始

- `experiment_id`: `20260723-03-official-metrics-reevaluation`.
- plan: [20260723-03-official-metrics-reevaluation-plan.md](plan/20260723-03-official-metrics-reevaluation-plan.md).
- 目的: 不重新采样, 使用 FLOWR.ROOT 官方 validity、完整 dock PoseBusters、strain 和条件子结构检查复评前两组各 50 个冻结结果.
- 解释边界: 官方条件匹配只要求固定子结构仍存在, 不等同于原子逐一映射或同分子坐标修复.
- 当前进度: 计划和实验骨架已建立, 即将执行测试、preflight、smoke 和正式 CPU 复评.
- 输出位置: `outputs/20260723-03-official-metrics-reevaluation/`.

## 2026-07-23: Official metrics reevaluation 完成

- `experiment_id`: `20260723-03-official-metrics-reevaluation`.
- 执行: 5 项测试、preflight、两样本 smoke 和两次完整 CPU 正式复评均通过; 两组各 50 个结果及 10 个基线全部保留.
- 结果: Coordinate-only validity `50/50`、条件匹配 `50/50`、PB-valid `14/50`; official inpainting validity `49/50`、条件匹配 `49/50`、PB-valid `24/50`.
- Strain: Coordinate-only 中位数 `7486.80 kcal/mol`, official inpainting 中位数 `145.80 kcal/mol`, `G_good` 基线中位数 `40.37 kcal/mol`.
- 解释: 官方条件匹配只做固定子结构搜索, 不要求原子逐一映射或同分子身份. 因此官方 PB-valid 表示生成质量, 不替代旧 strict 修复结论.
- 报告: [20260723-03-official-metrics-reevaluation-report.md](report/20260723-03-official-metrics-reevaluation-report.md).
- 输出: `outputs/20260723-03-official-metrics-reevaluation/`.

## 2026-07-23: Full-test oracle repair benchmark 开始

- `experiment_id`: `20260723-04-full-test-oracle-repair-benchmark`.
- plan: [20260723-04-full-test-oracle-repair-benchmark-plan.md](plan/20260723-04-full-test-oracle-repair-benchmark-plan.md).
- 目的: 在 SPINDR 官方 225-case test split 上自动构造可复现的中等强度扭转 clash, 使用 oracle mask 和固定 seed 完整比较 coordinate-only 与 official inpainting.
- 主要终点: strict same-molecule repair 为主要终点, 官方 validity、条件匹配和完整 dock PoseBusters 为次要终点.
- 输出位置: `outputs/20260723-04-full-test-oracle-repair-benchmark/`.

## 2026-07-23: Full-test oracle repair benchmark 完成

- `experiment_id`: `20260723-04-full-test-oracle-repair-benchmark`.
- 构造: 225 个测试条目中 211 个 eligible; 其余为 2 个基线 PB 失败、8 个无可旋转键、3 个无中等强度候选和 1 个无隔离 PB 候选.
- 执行: 两种方法各完成 211 case x 10 seed x 100 steps, 共 4,220 次采样, 无模型错误、补抽或筛选; 全量评价执行两次且主要计数一致.
- Strict: Coordinate-only `584/2110` (`27.68%`), official inpainting `424/2110` (`20.09%`); case-cluster 差值 `+7.58` 个百分点, 95% CI `[+2.61,+12.70]`.
- 官方质量: Coordinate-only `606/2110` (`28.72%`), official inpainting `1209/2110` (`57.30%`); official inpainting 显著更高.
- 解释: Coordinate-only 更适合小范围严格同分子修复, 但在 7 个以上可编辑原子时内部几何明显恶化; official inpainting 更适合允许换分子的局部重设计.
- 报告: [20260723-04-full-test-oracle-repair-benchmark-report.md](report/20260723-04-full-test-oracle-repair-benchmark-report.md).
- 输出: `outputs/20260723-04-full-test-oracle-repair-benchmark/`.
