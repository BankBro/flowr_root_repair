# Official Fragment Inpainting Control

**实验 ID:** `20260723-02-official-fragment-inpainting-control`

**状态:** 已完成, strict 结论为 `NO-GO`

**计划:** [../../docs/plan/20260723-02-official-fragment-inpainting-control-plan.md](../../docs/plan/20260723-02-official-fragment-inpainting-control-plan.md)

**输出:** [../../outputs/20260723-02-official-fragment-inpainting-control/](../../outputs/20260723-02-official-fragment-inpainting-control/)

**报告:** [../../docs/report/20260723-02-official-fragment-inpainting-control-report.md](../../docs/report/20260723-02-official-fragment-inpainting-control-report.md)

本实验在上一轮完全相同的 5 个 `G_bad`、真值 `M_fixed` 和 50 个 seed 上, 使用 FLOWR.ROOT 官方 fragment inpainting 流程生成完整配体, 并以共同评估器公平比较 coordinate-only 修复与官方局部重设计.

从仓库根目录运行:

```bash
export PYTHONPATH="$PWD"
conda run --no-capture-output -n flowr_root python \
  experiments/20260723-02-official-fragment-inpainting-control/run.py preflight
conda run --no-capture-output -n flowr_root python \
  experiments/20260723-02-official-fragment-inpainting-control/run.py smoke --steps 5
conda run --no-capture-output -n flowr_root python \
  experiments/20260723-02-official-fragment-inpainting-control/run.py smoke --steps 100
conda run --no-capture-output -n flowr_root python \
  experiments/20260723-02-official-fragment-inpainting-control/run.py formal
```

正式运行前必须依次通过 CPU 测试、preflight、5 步 smoke 和 100 步 smoke. `formal` 按 `system_id + seed` 断点续跑, 模型错误和无效分子均保留在 50 次固定分母中.

正式结果为 native `14/50`、strict `7/50`. Strict gate 为 `NO-GO`. Coordinate-only 对照经同一评估器重算为 native `14/50`、strict `13/50`.
