# Full-Test Oracle Repair Benchmark

**实验 ID:** `20260723-04-full-test-oracle-repair-benchmark`

**状态:** 已完成

**计划:** [../../docs/plan/20260723-04-full-test-oracle-repair-benchmark-plan.md](../../docs/plan/20260723-04-full-test-oracle-repair-benchmark-plan.md)

**输出:** [../../outputs/20260723-04-full-test-oracle-repair-benchmark/](../../outputs/20260723-04-full-test-oracle-repair-benchmark/)

本实验在 SPINDR 官方 test split 上自动构造中等强度局部 clash, 使用 oracle mask 公平比较冻结的 coordinate-only 与 official inpainting.

`run.py` 提供 `preflight`、`construct`、`smoke`、`formal` 和 `evaluate` 命令. 正式运行前必须先完成测试、construction 和两级 GPU smoke.

```bash
PYTHONPATH=. conda run -n flowr_root python experiments/20260723-04-full-test-oracle-repair-benchmark/run.py preflight
PYTHONPATH=. conda run -n flowr_root python experiments/20260723-04-full-test-oracle-repair-benchmark/run.py construct
PYTHONPATH=. conda run -n flowr_root python experiments/20260723-04-full-test-oracle-repair-benchmark/run.py smoke --steps 5 --device cuda
PYTHONPATH=. conda run -n flowr_root python experiments/20260723-04-full-test-oracle-repair-benchmark/run.py smoke --steps 100 --device cuda
PYTHONPATH=. conda run -n flowr_root python experiments/20260723-04-full-test-oracle-repair-benchmark/run.py formal --method all --device cuda
PYTHONPATH=. conda run -n flowr_root python experiments/20260723-04-full-test-oracle-repair-benchmark/run.py evaluate
```

正式结果为 211 个 eligible case x 10 seed x 2 methods. 结论与完整统计见 [实验报告](../../docs/report/20260723-04-full-test-oracle-repair-benchmark-report.md).
