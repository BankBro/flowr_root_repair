# 实验日志

本文件按时间顺序追加正式实验和会影响研究决策的阶段记录. 不重写历史条目; 需要纠正时追加更正说明.

当前尚无正式实验记录.

## 建议记录内容

每条记录至少说明日期, `experiment_id`, 对应 plan, 目的, 执行内容, 输出位置, 结果或失败原因以及下一步. 具体格式可以根据实验灵活调整.

## 2026-07-23: Oracle-mask repair pilot 开始

- `experiment_id`: `20260723-01-oracle-mask-repair-pilot`.
- plan: [20260723-01-oracle-mask-repair-pilot-plan.md](plan/20260723-01-oracle-mask-repair-pilot-plan.md).
- 目的: 在 5 个冻结的 SPINDR validation 人工 clash 案例上, 用 50 次真值掩码 rollout 验证 coordinate-only 修复器的基本可行性.
- 当前进度: 协议与案例已经冻结, experiment/output 骨架已建立, 即将实现 repair 内核并执行 preflight 和 smoke.
- 输出位置: `outputs/20260723-01-oracle-mask-repair-pilot/`.
