# Official Metrics Reevaluation Outputs

**实验 ID:** `20260723-03-official-metrics-reevaluation`

**状态:** 已完成

**来源:** [../../experiments/20260723-03-official-metrics-reevaluation/](../../experiments/20260723-03-official-metrics-reevaluation/)

本目录保存两组冻结 SDF 的官方核心指标逐样本结果、基线、输入清单和汇总. 不复制源 SDF、checkpoint 或数据集.

主要产物:

- `runs.csv`: 两组共 100 个生成结果的逐样本官方指标.
- `baselines.csv`: 5 个 `G_good` 和 5 个 `G_bad` 的官方指标.
- `summary.json`: 总体、逐 case、配对和 PB 子项汇总.
- `input_manifest.json`: 所有输入 SDF 和 PDB 的路径、大小与 SHA-256.
- `preflight.json`: 输入完整性和评估配置审计.

正式结果: Coordinate-only PB-valid `14/50`, official inpainting PB-valid `24/50`. 两组 conditioned PB-valid 与 PB-valid 相同. 这表示官方生成质量口径下的可用产物率, 不表示严格同分子修复成功率.
