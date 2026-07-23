# Official Fragment Inpainting Control Outputs

**实验 ID:** `20260723-02-official-fragment-inpainting-control`

**状态:** 已完成, native `14/50`, strict `7/50`, strict gate 为 `NO-GO`

**来源:** [../../experiments/20260723-02-official-fragment-inpainting-control/](../../experiments/20260723-02-official-fragment-inpainting-control/)

本目录保存官方 fragment inpainting 的 50 个 SDF、逐次结果、共同评估汇总、两组比较漏斗和 SDF 重读审计. `coordinate_only_runs.csv` 是旧组由共同评估器得到的冻结对照, `official_runs.csv` 是官方组最终 artifact 口径结果.

正式结果为 native `14/50`、strict `7/50`. 官方组和 coordinate-only 组的 native 总数相同, 但 strict 比旧组的 `13/50` 少 6 次. Checkpoint、数据集副本、采样轨迹、缓存和长期日志不进入 Git.
