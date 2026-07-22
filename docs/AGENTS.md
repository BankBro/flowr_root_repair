# docs/AGENTS.md

## 文档职责

- `plan/` 保存研究计划, 模型设计和实验方案.
- `report/` 保存完成, 失败或中止实验的结果, 限制和复盘.
- `STATUS.md` 只记录当前阶段, 当前计划和下一步, 保持简短.
- `EXPERIMENT_LOG.md` 按时间追加实验过程与阶段结论. 历史记录有误时追加更正, 不直接覆盖.

## 命名与对应关系

- plan 使用 `YYYYMMDD-NN-<slug>-plan.md`, report 使用 `YYYYMMDD-NN-<slug>-report.md`.
- `NN` 是当天的两位数顺序号. plan 和 report 的日期与编号独立分配, 对应同一实验时使用相同 slug, 并在正文中写明同一个 `experiment_id`.
- plan, experiment, report 一一对应. 一个 plan 可以规划同一总实验内部的多个子阶段, 不要求每个子阶段单独建立文档.

## 内容维护

- 计划开始指导实现后, 保留其当时的目标, 边界和验收标准. 实施偏差写入对应 report, 重大路线变化通过新 plan 处理.
- 结论必须区分已验证结果, 初步证据, 假设和待验证事项.
- 不伪造实验结果, 命令, 引用或完成状态. 失败和负结果同样需要记录.
- 模型和训练细节以当前 plan 为准, 不复制到 `AGENTS.md` 中形成第二份规格.
- 重要阶段变化后更新 `STATUS.md`. 正式实验开始, 完成, 失败或停止时, 在 `EXPERIMENT_LOG.md` 追加记录.
