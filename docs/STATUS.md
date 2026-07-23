# 项目状态

**更新时间:** 2026-07-23

## 当前阶段

- 首个 oracle-mask repair pilot 已完成, 冻结 gate 结论为 `CONDITIONAL`.
- 当前计划: [20260723-02-official-fragment-inpainting-control-plan.md](plan/20260723-02-official-fragment-inpainting-control-plan.md).
- 上位设计: [20260722-01-flowr-root-ligand-repair-v1-plan.md](plan/20260722-01-flowr-root-ligand-repair-v1-plan.md).
- 当前工作: 在同一批人工 clash 输入和 seed 上执行 FLOWR.ROOT 官方 fragment inpainting, 与已冻结的 coordinate-only 结果做公平对照.
- 对照同时报告局部重设计成功和严格同分子修复成功, 避免把生成的新分子误计为原分子的坐标修复.

## 已确认边界

- FLOWR.ROOT 原有行为保持为默认路径, 修复模式显式启用.
- 首版只新增 Fixed Mask Head, 并复用 Structure Head 修复坐标.
- `M_fixed = 1` 表示固定, `M_fixed = 0` 表示修复, 与有效原子范围内的 `fragment_mask` 方向一致.
- 先验证真值掩码下的修复上限, 再训练掩码预测并进入自动修复.
- plan, experiment, report 一一对应, 当前首版计划作为一个包含多个子阶段的总实验.

## 下一步

1. 完成共同评估器和真值 fragment mask 适配器.
2. 通过旧结果复现检查与 GPU smoke 后运行 50 次官方采样.
3. 根据 native/strict 双终点判断官方方法适合作为局部重设计还是同分子修复基线.
