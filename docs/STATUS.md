# 项目状态

**更新时间:** 2026-07-23

## 当前阶段

- 首个 oracle-mask repair pilot 已完成, 冻结 gate 结论为 `CONDITIONAL`.
- 当前计划: [20260723-01-oracle-mask-repair-pilot-plan.md](plan/20260723-01-oracle-mask-repair-pilot-plan.md).
- 上位设计: [20260722-01-flowr-root-ligand-repair-v1-plan.md](plan/20260722-01-flowr-root-ligand-repair-v1-plan.md).
- 首版范围: 诊断蛋白质-配体 clash, 定位需要移动的配体原子, 并使用 FLOWR.ROOT 进行 coordinate-only 局部修复.
- 正式实验 50/50 正常完成, 14/50 成功. 固定坐标和离散图保持 50/50, 主要瓶颈是可编辑区域的内部几何.
- 实验报告: [20260723-01-oracle-mask-repair-pilot-report.md](report/20260723-01-oracle-mask-repair-pilot-report.md).

## 已确认边界

- FLOWR.ROOT 原有行为保持为默认路径, 修复模式显式启用.
- 首版只新增 Fixed Mask Head, 并复用 Structure Head 修复坐标.
- `M_fixed = 1` 表示固定, `M_fixed = 0` 表示修复, 与有效原子范围内的 `fragment_mask` 方向一致.
- 先验证真值掩码下的修复上限, 再训练掩码预测并进入自动修复.
- plan, experiment, report 一一对应, 当前首版计划作为一个包含多个子阶段的总实验.

## 下一步

1. 新建独立 plan, 在 oracle mask 条件下验证能够保持局部共价几何的修复表示或约束.
2. 保留本轮冻结协议作为对照, 不在已完成实验内调参.
3. Oracle 修复达到 GO 标准后, 再进入 Fixed Mask Head 训练.
