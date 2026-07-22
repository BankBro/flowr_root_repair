# Outputs

本目录保存实验产生的结果, 使用 `outputs/<experiment_id>/` 与 `experiments/<experiment_id>/` 一一对应.

输出目录结构可以根据实验灵活安排, 但必须包含简短的 `README.md`, 说明来源实验, 主要内容, 关键结果和保留策略. 每个输出目录还应按实际产物维护自己的 `.gitignore`.

默认只向 Git 提交小型总结, 必要表格和关键复现信息. checkpoint, 大规模逐样本结果, 日志, 缓存和可再生中间产物应默认忽略. 当前尚未创建具体输出目录.
