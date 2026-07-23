# Oracle-mask Repair Pilot

**实验 ID:** `20260723-01-oracle-mask-repair-pilot`

**状态:** 已完成, 结论为 `CONDITIONAL`

**计划:** [../../docs/plan/20260723-01-oracle-mask-repair-pilot-plan.md](../../docs/plan/20260723-01-oracle-mask-repair-pilot-plan.md)

**输出:** [../../outputs/20260723-01-oracle-mask-repair-pilot/](../../outputs/20260723-01-oracle-mask-repair-pilot/)

**报告:** [../../docs/report/20260723-01-oracle-mask-repair-pilot-report.md](../../docs/report/20260723-01-oracle-mask-repair-pilot-report.md)

本实验在 5 个冻结的 SPINDR validation 复合物上, 使用真值 `M_fixed` 测试 FLOWR.ROOT coordinate-only 局部修复. `cases.json` 是人工缺陷的唯一正式定义.

从仓库根目录运行:

```bash
export PYTHONPATH="$PWD"
conda run --no-capture-output -n flowr_root python \
  experiments/20260723-01-oracle-mask-repair-pilot/run.py preflight
conda run --no-capture-output -n flowr_root python \
  experiments/20260723-01-oracle-mask-repair-pilot/run.py smoke --steps 5
conda run --no-capture-output -n flowr_root python \
  experiments/20260723-01-oracle-mask-repair-pilot/run.py smoke --steps 100
conda run --no-capture-output -n flowr_root python \
  experiments/20260723-01-oracle-mask-repair-pilot/run.py formal
```

`run.py formal` 会按 `system_id + seed` 断点续跑. 正式协议固定为 5 个 case, 每个 case 10 个 seed, 共 50 次采样.

正式结果为 14/50 成功. 50/50 保持固定坐标与离散图, 但内部几何只通过 20/50, 因而未达到 GO 标准.
