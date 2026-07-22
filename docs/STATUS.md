# 项目状态

**更新时间:** 2026-07-22

## 当前阶段

- 首版修复方案已经完成设计整理, 尚未开始代码实现或正式实验.
- 当前计划: [20260722-01-flowr-root-ligand-repair-v1-plan.md](plan/20260722-01-flowr-root-ligand-repair-v1-plan.md).
- 首版范围: 诊断蛋白质-配体 clash, 定位需要移动的配体原子, 并使用 FLOWR.ROOT 进行 coordinate-only 局部修复.
- 仓库协作规则和顶层实验框架已经建立. 当前没有具体的 experiment 或 output 子目录.

## 已确认边界

- FLOWR.ROOT 原有行为保持为默认路径, 修复模式显式启用.
- 首版只新增 Fixed Mask Head, 并复用 Structure Head 修复坐标.
- `M_fixed = 1` 表示固定, `M_fixed = 0` 表示修复, 与有效原子范围内的 `fragment_mask` 方向一致.
- 先验证真值掩码下的修复上限, 再训练掩码预测并进入自动修复.
- plan, experiment, report 一一对应, 当前首版计划作为一个包含多个子阶段的总实验.

## 下一步

1. 单独确认项目运行环境和依赖可用性.
2. 准备开始实验时, 创建与当前计划对应的 experiment 和 output 目录.
3. 先进行轻量数据检查和 oracle-mask smoke test, 通过后再进入正式实验.
