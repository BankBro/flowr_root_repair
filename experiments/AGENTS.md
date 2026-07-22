# experiments/AGENTS.md

## 实验目录规则

- 正式实验使用 `experiment_id = YYYYMMDD-NN-<experiment_name>`, 其中 `NN` 为两位数.
- 每个 experiment 必须对应一个 plan, 并在完成, 失败或中止后对应一个 report.
- 实验专用脚本, 命令, 配置和必要的中间内容放入 `experiments/<experiment_id>/`. 对应结果放入 `outputs/<experiment_id>/`.
- 每个实验目录必须有简短的 `README.md`, 写清目的, 对应 plan, 当前状态, 主要入口和输出位置.
- 目录结构可以按实验需要灵活组织. 不强制 metadata, schema, `configs/` 或统一子目录模板.
- 每个实验根目录按实际产物维护局部 `.gitignore`. 默认提交轻量, 可复用或支撑结论的内容, 忽略大文件, 缓存, 日志和可再生中间结果.
- 正式运行前先完成轻量检查或 smoke test. 用于研究决策的 smoke 结果应在实验 README, 日志或报告中留下可理解的记录.
- 实验状态变化后更新 `docs/STATUS.md`, 并在 `docs/EXPERIMENT_LOG.md` 追加记录.
