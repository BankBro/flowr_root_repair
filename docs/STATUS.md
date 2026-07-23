# 项目状态

**更新时间:** 2026-07-23

## 当前阶段

- Official fragment inpainting control 已完成, strict gate 结论为 `NO-GO`.
- 最近计划: [20260723-02-official-fragment-inpainting-control-plan.md](plan/20260723-02-official-fragment-inpainting-control-plan.md).
- 上位设计: [20260722-01-flowr-root-ligand-repair-v1-plan.md](plan/20260722-01-flowr-root-ligand-repair-v1-plan.md).
- 同一评估器下, coordinate-only 为 native `14/50`、strict `13/50`; 官方完整 inpainting 为 native `14/50`、strict `7/50`.
- 官方方法能在部分案例生成有效的新局部结构, 但不能作为严格同分子坐标修复器的直接替代.
- 实验报告: [20260723-02-official-fragment-inpainting-control-report.md](report/20260723-02-official-fragment-inpainting-control-report.md).

## 已确认边界

- FLOWR.ROOT 原有行为保持为默认路径, 修复模式显式启用.
- 首版只新增 Fixed Mask Head, 并复用 Structure Head 修复坐标.
- `M_fixed = 1` 表示固定, `M_fixed = 0` 表示修复, 与有效原子范围内的 `fragment_mask` 方向一致.
- 先验证真值掩码下的修复上限, 再训练掩码预测并进入自动修复.
- plan, experiment, report 一一对应, 当前首版计划作为一个包含多个子阶段的总实验.

## 下一步

1. 保留两个已冻结实验作为后续 oracle-mask 对照, 不在原实验内调参.
2. 新建独立 plan, 验证保持同一分子和局部共价几何的坐标表示、约束或采样方式.
3. Oracle 修复达到 GO 标准后, 再进入 Fixed Mask Head 训练.
