# Official Fragment Inpainting Control

**实验 ID:** `20260723-02-official-fragment-inpainting-control`

**状态:** 执行中

**计划:** [../../docs/plan/20260723-02-official-fragment-inpainting-control-plan.md](../../docs/plan/20260723-02-official-fragment-inpainting-control-plan.md)

**输出:** [../../outputs/20260723-02-official-fragment-inpainting-control/](../../outputs/20260723-02-official-fragment-inpainting-control/)

本实验在上一轮完全相同的 5 个 `G_bad`、真值 `M_fixed` 和 50 个 seed 上, 使用 FLOWR.ROOT 官方 fragment inpainting 流程生成完整配体, 并以共同评估器公平比较 coordinate-only 修复与官方局部重设计.

主要入口和命令会随实现补充. 正式运行前必须依次通过 CPU 测试、preflight、5 步 smoke 和 100 步 smoke.

