# Experiments

本目录保存单次实验专用的脚本, 命令, 配置和必要的中间内容.

正式实验使用 `experiments/<experiment_id>/`, 其中 `experiment_id = YYYYMMDD-NN-<experiment_name>`. 每个实验必须对应一个 plan, 并使用相同 ID 的 `outputs/<experiment_id>/` 保存结果.

实验目录结构可以灵活安排, 但必须包含简短的 `README.md`, 使读者能够看懂实验目的, 当前状态, 主要入口和输出位置. 每个实验还应按实际产物维护自己的 `.gitignore`.

正式实验开始前先完成轻量 CPU 或 GPU 检查以及必要的 smoke test. 当前尚未创建具体实验目录.
