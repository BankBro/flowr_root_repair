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
