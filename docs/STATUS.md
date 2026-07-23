# 项目状态

**更新时间:** 2026-07-23

## 当前阶段

- 首版修复方案已经完成设计整理, 当前进入首个 oracle-mask repair pilot.
- 当前计划: [20260723-01-oracle-mask-repair-pilot-plan.md](plan/20260723-01-oracle-mask-repair-pilot-plan.md).
- 上位设计: [20260722-01-flowr-root-ligand-repair-v1-plan.md](plan/20260722-01-flowr-root-ligand-repair-v1-plan.md).
- 首版范围: 诊断蛋白质-配体 clash, 定位需要移动的配体原子, 并使用 FLOWR.ROOT 进行 coordinate-only 局部修复.
- 已创建 `20260723-01-oracle-mask-repair-pilot` 的 experiment 和 output 骨架, 正在实现 coordinate-only oracle 修复路径.

## 已确认边界

- FLOWR.ROOT 原有行为保持为默认路径, 修复模式显式启用.
- 首版只新增 Fixed Mask Head, 并复用 Structure Head 修复坐标.
- `M_fixed = 1` 表示固定, `M_fixed = 0` 表示修复, 与有效原子范围内的 `fragment_mask` 方向一致.
- 先验证真值掩码下的修复上限, 再训练掩码预测并进入自动修复.
- plan, experiment, report 一一对应, 当前首版计划作为一个包含多个子阶段的总实验.

## 下一步

1. 实现并测试可复用的 coordinate-only oracle 修复内核.
2. 完成 CPU preflight 与 5 步、100 步 GPU smoke.
3. 运行冻结的 5 个样本 x 10 个 seed 正式实验并形成 report.
