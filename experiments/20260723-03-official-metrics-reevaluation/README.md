# Official Metrics Reevaluation

**实验 ID:** `20260723-03-official-metrics-reevaluation`

**状态:** 已完成

**计划:** [../../docs/plan/20260723-03-official-metrics-reevaluation-plan.md](../../docs/plan/20260723-03-official-metrics-reevaluation-plan.md)

**输出:** [../../outputs/20260723-03-official-metrics-reevaluation/](../../outputs/20260723-03-official-metrics-reevaluation/)

**报告:** [../../docs/report/20260723-03-official-metrics-reevaluation-report.md](../../docs/report/20260723-03-official-metrics-reevaluation-report.md)

本实验不重新采样, 直接使用 FLOWR.ROOT 官方 validity、完整 dock PoseBusters、strain 和条件子结构检查, 复评前两次实验已经保存的两组各 50 个 SDF.

从仓库根目录运行:

```bash
export PYTHONPATH="$PWD"
conda run --no-capture-output -n flowr_root python \
  experiments/20260723-03-official-metrics-reevaluation/run.py preflight
conda run --no-capture-output -n flowr_root python \
  experiments/20260723-03-official-metrics-reevaluation/run.py smoke
conda run --no-capture-output -n flowr_root python \
  experiments/20260723-03-official-metrics-reevaluation/run.py formal
```

正式复评固定保留每组 50 次分母. 无效分子、PB 失败或 strain 不可计算均不得从原始分母中删除.

正式结果: Coordinate-only validity `50/50`、条件匹配 `50/50`、PB-valid `14/50`; official inpainting validity `49/50`、条件匹配 `49/50`、PB-valid `24/50`. Official 条件匹配只表示固定子结构仍能找到, 不表示同一个分子或原子逐一保持.
