# 项目状态

**更新时间:** 2026-07-23

## 当前阶段

- SPINDR 全测试集 oracle repair benchmark 已完成, 225 个条目中 211 个可自动构造中等强度扭转 clash.
- 最近计划: [20260723-04-full-test-oracle-repair-benchmark-plan.md](plan/20260723-04-full-test-oracle-repair-benchmark-plan.md).
- 上位设计: [20260722-01-flowr-root-ligand-repair-v1-plan.md](plan/20260722-01-flowr-root-ligand-repair-v1-plan.md).
- 211 case x 10 seed 下, coordinate-only strict 为 `584/2110`, official inpainting strict 为 `424/2110`; case-cluster 95% CI 支持 coordinate-only 更高.
- 官方质量下, coordinate-only 为 `606/2110`, official inpainting 为 `1209/2110`; official inpainting 生成质量更好, 但大量结果改变完整分子身份.
- 实验报告: [20260723-04-full-test-oracle-repair-benchmark-report.md](report/20260723-04-full-test-oracle-repair-benchmark-report.md).

## 已确认边界

- FLOWR.ROOT 原有行为保持为默认路径, 修复模式显式启用.
- 首版只新增 Fixed Mask Head, 并复用 Structure Head 修复坐标.
- `M_fixed = 1` 表示固定, `M_fixed = 0` 表示修复, 与有效原子范围内的 `fragment_mask` 方向一致.
- 先验证真值掩码下的修复上限, 再训练掩码预测并进入自动修复.
- plan, experiment, report 一一对应, 当前首版计划作为一个包含多个子阶段的总实验.

## 下一步

1. 冻结本次 211-case benchmark 作为后续 oracle-mask 对照, 不在原实验内调参.
2. 新建独立 plan, 优先解决 7 个以上可编辑原子时的键长、键角、internal clash 和 strain 问题.
3. 同时约束立体化学身份和环/边界共价几何; oracle strict 修复明显提高后, 再进入 Fixed Mask Head 训练.
