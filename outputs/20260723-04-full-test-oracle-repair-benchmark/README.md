# Full-Test Oracle Repair Benchmark Outputs

**实验 ID:** `20260723-04-full-test-oracle-repair-benchmark`

**状态:** 已完成

**来源:** [../../experiments/20260723-04-full-test-oracle-repair-benchmark/](../../experiments/20260723-04-full-test-oracle-repair-benchmark/)

本目录保存 225-case 构造漏斗、211 个 eligible 输入、两种方法各 2,110 次正式结果、官方/strict 评价和汇总. 大规模 SDF、日志、缓存与可再生 smoke 结果默认忽略, 清单、CSV 和 JSON 保留.

主要入口:

- `construction.csv`: 全部 225 个测试条目的构造漏斗.
- `sampling_runs.csv`: 4,220 次正式采样库存.
- `runs.csv`: 逐次官方与 strict 评价.
- `case_rates.csv`: 211 个案例的 10-seed 成功率.
- `summary.json`: 总体统计、累计漏斗和配对 bootstrap.
- `output_manifest.json`: 正式产物的大小与 SHA-256.
- [实验报告](../../docs/report/20260723-04-full-test-oracle-repair-benchmark-report.md).
