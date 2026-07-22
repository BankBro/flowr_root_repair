# 基于 FLOWR.ROOT 的配体失败诊断与局部生成修复

- **版本:** 1.8, 分层精简版
- **基线代码:** FLOWR.ROOT commit `b2263e2516ad798d0119f8c4b531698860cc846e`
- **首版范围:** 蛋白质-配体 clash 的局部坐标修复

> 核心目标: 先找到需要移动的配体原子, 再利用 FLOWR.ROOT 只重建这些原子的坐标, 最后由外部验证器决定接受、拒绝或继续修复.

## 0. 一页总览

### 0.1 系统做什么

| 类别 | 内容 |
|---|---|
| 输入 | 蛋白质口袋、完整候选配体、可选的用户固定掩码 |
| 输出 | `M_fixed`、修复后的完整配体、clash 验证结果 |

首版只处理蛋白质-配体空间碰撞, 只允许修改坐标. 原子数量、顺序、类型、键、电荷和杂化状态全部保持不变.

### 0.2 主流程

```text
完整候选配体
-> diagnose()
-> M_fixed
-> 若全部有效原子均为 1, 直接返回
-> 在 M_fixed = 0 的区域构造局部坐标 prior
-> repair_step() 多步积分
-> 每步恢复固定坐标和全部离散属性
-> 外部验证
-> 接受、回滚或重新诊断
```

### 0.3 首版模型

```text
蛋白质口袋 -> FLOWR.ROOT Pocket Encoder
当前配体   -> FLOWR.ROOT Ligand Decoder
                         |
                         +-> 原有 Structure Head, 负责坐标修复
                         +-> 新增 Fixed Mask Head, 负责定位修复区域
```

首版只新增一个 Fixed Mask Head. Affinity Head 和 Confidence Head 保留但冻结. 其他预测头均后置.

### 0.4 唯一掩码

| 位置 | `ligand_mask` | `M_fixed` | 含义 |
|---|---:|---:|---|
| 有效固定原子 | 1 | 1 | 坐标保持不变 |
| 有效修复原子 | 1 | 0 | 坐标参与局部生成 |
| padding | 0 | 0 | 不解释修复语义 |

在有效原子范围内, `M_fixed` 与 FLOWR.ROOT 原有 `fragment_mask` 方向一致. padding 位置清零后可直接作为 `fragment_mask` 使用.

用户掩码也采用相同语义. 首版只做阈值化、padding 过滤和基本合法性检查, 不强制连通, 也不自动扩展芳香环或刚性基团. 初始训练优先使用单个连通修复区域.

### 0.5 推荐训练顺序

```text
数据与验证器
-> 真值掩码下验证冻结修复器
-> 必要时 LoRA 适配修复器
-> 冻结主干训练 Fixed Mask Head
-> predicted-mask 单轮自动修复
-> 掩码鲁棒性训练
-> 可选的交替联合适配
-> 外层多轮修复
```

其中真值掩码实验也称 oracle-mask experiment, 即直接使用人工 corruption 已知的 `M_fixed_gt`, 不让 Fixed Mask Head 预测.

## 1. 任务定义

### 1.1 关键对象

- `G_good`: 作为监督终点的参考合格配体构象, 不要求是唯一正确构象.
- `G_bad`: 人工制造 clash 后的失败构象.
- `G_0`: 修复区域从坐标 prior 采样得到的起始状态.
- `G_t`: flow 时间 `t` 上的局部带噪状态.
- `G_pred`: Structure Head 预测的最终完整配体.
- `M_fixed`: 1 固定、0 修复的二值掩码.

### 1.2 首版边界

首版采用 coordinate-only 硬覆盖:

- 只修改 `M_fixed = 0` 的有效原子坐标.
- `M_fixed = 1` 的坐标在每个积分步骤后恢复.
- 所有原子类型、键、电荷和杂化状态在每一步后恢复.
- 最终输出再次硬覆盖.
- clash 和 PoseBusters 只用于外部验证, 不直接反向传播.

以下内容后置:

- 多错误类型.
- 拓扑和键修复.
- 可修复性预测.
- 不依赖硬覆盖的部分 Flow Matching.
- 显式轨迹记忆和测试时训练.

## 2. 两条模型路径

### 2.1 诊断路径

接口:

```text
diagnose(pocket, complete_ligand) -> M_fixed
```

| 项目 | 设计 |
|---|---|
| 配体输入 | 完整、未加噪候选结构 |
| 坐标时间 | 所有有效原子均为 1 |
| 离散时间 | 所有有效原子均为 1 |
| self-conditioning | 不使用 |
| 正式输出 | `fixed_logits` 和硬化后的 `M_fixed` |
| Structure Head | 可以计算, 但结果忽略 |

诊断前向不能读取 `M_fixed_gt`. 如果旧接口强制要求 `fragment_mask`, 只能传与样本标签无关的全 0 占位值.

Fixed Mask Head 读取 Ligand Decoder 的逐原子不变特征 `invs_norm [B, N, d_inv]`. 最小实现见附录 D.

### 2.2 修复路径

接口:

```text
repair_step(pocket_repr, G_t, t, M_fixed, self_conditioning) -> G_pred
```

| 项目 | 设计 |
|---|---|
| 配体输入 | 局部带噪 `G_t` |
| 坐标时间 | 固定原子为 1, 修复原子为 `t` |
| 离散时间 | 所有有效原子均为 1 |
| 正式输出 | Structure Head 的完整坐标预测 |
| Fixed Mask Head | 结果忽略 |

单轮采样内 `M_fixed` 保持不变. Fixed Mask Head 不在内层 flow 时间步更新掩码. 下一轮掩码只能由最终完整结果重新进入 `diagnose()` 获得.

单轮修复中 Pocket Encoder 表示通常只计算一次并缓存. 一次 `repair_step()` 只重新估计终点和修复方向, 完整修复需要多个积分时间步.

### 2.3 训练路由

训练器在模型外区分两类 batch:

```text
diagnosis batch -> diagnose()    -> mask loss
repair batch    -> repair_step() -> local coordinate loss
```

训练路由不写入样本, 也不输入神经网络. 首版先分阶段训练; 稳定后如需共享层联合适配, 仍使用两个独立前向交替更新.

## 3. 数据设计

### 3.1 基础数据

从经过预处理的蛋白质-配体复合物开始:

```text
protein + G_good
```

预处理必须保证:

- 蛋白质和配体处于同一坐标系.
- altloc、显式氢和质子化流程统一.
- 沿用 FLOWR.ROOT 的口袋截取方式, 首版可使用 7 Å.
- RDKit 能读取和 sanitization 配体.
- 原子类型位于模型词表中.
- 原子顺序稳定并保存映射.

必须先按原始复合物划分训练、验证和测试集, 再生成 corruption. 同一原始复合物的所有衍生样本必须留在同一集合. 测试集还应尽量隔离蛋白质家族和配体 scaffold.

### 3.2 有 clash 的失败样本

1. 选择非环重原子单键.
2. 排除酰胺键、强共轭键和无有效重原子移动的末端键.
3. 确定键一侧的移动片段.
4. 枚举二面角并选择产生目标 clash 的构象.
5. 保证未移动区域坐标不变.
6. 移动片段标记为 `M_fixed_gt = 0`, 其余有效原子标记为 1.

直接发生 clash 的原子只用于分析. 训练标签使用完成修复所需的完整动作区域, 而不只标记碰撞原子.

### 3.3 诊断样本的三种类别

| 类别 | 输入 | `M_fixed_gt` | 目的 |
|---|---|---|---|
| 原始参考结构 | `G_good` | 全部有效原子为 1 | 学习正常时不触发修复 |
| clash 正样本 | `G_bad` | 动作区域为 0 | 学习定位修复区域 |
| 无 clash 扭转对照 | `G_torsion_no_clash` | 全部有效原子为 1 | 防止模型只识别“发生过扭转” |

无 clash 扭转对照至少应满足没有蛋白质-配体 clash、没有配体内部 clash, 且局部几何和应变处于合理范围. 全固定只表示按照首版 clash 标准不触发修复, 不表示该构象一定是实验真实姿势.

诊断训练初期在样本层面对三类数据近似平衡, 再根据部署分布和验证集校准.

### 3.4 两种样本契约

```text
诊断: (protein, G_candidate, M_fixed_gt)
修复: (protein, G_bad, G_good, M_fixed_gt)
```

诊断输入不包含 `M_fixed_gt`, 该字段只用于 mask loss. repair batch 只包含至少一个有效修复原子的 clash 样本. 全固定样本不进入局部坐标 loss.

## 4. 分阶段训练

### 阶段 0: 数据和验证器

不训练模型. 完成数据划分、corruption、困难负样本、原子映射、clash 检查、局部几何检查和可复现实验记录.

### 阶段 1: 真值掩码下验证冻结修复器

直接把 `M_fixed_gt` 交给冻结的 FLOWR.ROOT:

```text
G_bad + M_fixed_gt
-> coordinate-only 多步采样
-> 外部验证
```

要求:

- prior 只由 `G_bad`、固定区域和 pocket 构造.
- 所有离散属性保持输入值, 离散时间统一为 1.
- 每步和最终都执行硬覆盖.
- 随机 prior 可运行多个 rollout 估计成功率.

如果真值掩码下也无法稳定修复, 问题位于数据构造或修复器, 此时不训练 Fixed Mask Head.

### 阶段 1.5: 必要时适配修复器

仅在冻结修复器的真值掩码上限不足时进入:

1. 根据 `M_fixed_gt` 得到固定区域和修复区域.
2. 固定区域复制 `G_bad`, 修复区域从局部 prior 采样.
3. 为每个样本采样 `t`, 只对修复坐标构造 `G_t`.
4. 模型预测最终完整坐标 `G_pred`.
5. 只在修复区域计算坐标 loss.

参数策略:

- 冻结 Pocket Encoder.
- 优先对 Ligand Decoder 后部加入 LoRA.
- 必要时微调坐标相关 Structure Head.
- 冻结离散 Structure Head、Affinity Head 和 Confidence Head.
- 保留原有 self-conditioning.
- 不使用改变原子身份的 permutation alignment.

首版主线使用以 `G_bad` 修复区域为参考的局部各向同性 prior. 各向异性 prior 和 `G_bad + noise` 只作为消融. 训练和推理必须使用相同规则.

### 阶段 2: 训练 Fixed Mask Head

先冻结已验证的修复器和共享主干, 只训练 Fixed Mask Head:

1. 输入完整未加噪的 `G_candidate`.
2. 所有有效原子的坐标和离散时间设为 1.
3. 不输入 `M_fixed_gt`, 不使用 self-conditioning.
4. 计算逐原子加权 mask loss.

如果只训练新 Head 不足, 才小范围解冻 Decoder 最后若干层或加入 LoRA. 此时必须持续回放真值掩码修复评估, 防止共享层更新损害修复能力.

### 阶段 3: predicted-mask 单轮修复

```text
G_candidate
-> diagnose()
-> predicted M_fixed
-> 阈值化并清除 padding
-> 构造局部 prior
-> 单轮 repair
-> 外部验证
```

暂不通过阈值和采样过程端到端反向传播. 阈值根据 Repair Success Rate、正常样本误触发率和平均修复区域大小共同校准, 不只依据 mask F1.

### 阶段 4: 掩码鲁棒性

repair batch 逐步混入:

- 真值掩码.
- 多放开一跳邻居的扩大区域.
- 少选少量边界原子的欠覆盖区域.
- 同大小偏移区域.
- `detach` 后硬化的模型预测掩码.
- 后期再加入多个不连通区域.

预测掩码不反向传播到 Fixed Mask Head. 首版优先提高修复区域召回率, 同时用区域大小和正常样本误触发率约束过度修复.

### 阶段 5: 可选的交替联合适配

只有诊断器和修复器分别稳定后才考虑:

```text
diagnosis batch -> mask loss
repair batch    -> local coordinate loss
```

两个 batch 使用独立前向并交替更新. 共享层使用较小学习率, repair batch 继续以真值和受控扰动掩码为主. 如果真值掩码修复指标下降, 立即回到分阶段 checkpoint.

### 阶段 6: 外层循环

一轮修复后重新诊断. 验证通过且所有有效原子均固定时停止; 仍有修复区域时开始下一轮; 结果变差则回滚. 初始只做 1 轮, 稳定后扩展到 2 至 3 轮.

## 5. 损失函数

| 训练路径 | 主损失 | 计算范围 |
|---|---|---|
| Fixed Mask Head | 加权 `BCEWithLogitsLoss` | `ligand_mask = 1` 的原子 |
| 修复器 | 最终坐标 MSE | `ligand_mask = 1` 且 `M_fixed = 0` |

Mask loss 设计:

- 对标签为 0 的修复原子给予更高权重.
- 在样本层面平衡三种诊断类别.
- focal loss 和 soft Dice 只作为消融.
- 掩码大小正则只使用小权重.

Repair loss 设计:

- Structure Head 预测最终完整坐标 `G_pred`, 以 `G_good` 作为监督终点.
- 首版直接计算 `G_pred` 与 `G_good` 在修复区域的坐标 MSE, 不使用 velocity loss.
- 按每个样本的实际修复原子数归一化.
- 可选加入跨越固定区与修复区边界的键长和键角 loss.
- 不计算原子类型、键、电荷、杂化、Affinity 或 Confidence loss.
- 固定区域 RMSD 和属性一致性只作为强监控指标.

全固定样本不进入 repair loss. 首版不在同一次前向中直接相加 mask loss 和 repair loss.

## 6. 推理与验证

### 6.1 单轮推理

1. 编码蛋白质口袋.
2. 对完整候选配体运行 Fixed Mask Head.
3. 全部有效原子均固定时直接返回.
4. 只依据当前输入配体、固定区域和 pocket 构造局部 prior.
5. 使用同一个 `M_fixed` 运行整轮多步 FLOWR.ROOT.
6. 固定坐标时间为 1, 修复坐标时间为 `t`, 离散时间全部为 1.
7. 每步恢复固定坐标和全部离散属性.
8. 最终再次硬覆盖.
9. 运行外部验证并决定接受或回滚.

### 6.2 首版通过条件

- 目标蛋白质-配体 clash 消失.
- 没有产生新的蛋白质-配体或配体内部 clash.
- RDKit sanitization 通过.
- 化学图完全不变.
- 固定原子最大位移低于数值容差.
- 局部键长、键角和构象基本合理.

PoseBusters、strain energy 和关键相互作用保持作为辅助指标. Docking score 不作为唯一硬门槛.

## 7. 评价和实验门槛

### 7.1 核心指标

| 类别 | 指标 |
|---|---|
| 掩码 | 修复原子 AUPRC、precision、recall、F1、IoU |
| 安全性 | 正常样本和无 clash 对照的误触发率 |
| 修复 | Repair Success Rate、clash 消除率、新问题产生率 |
| 保持性 | 化学图保持率、固定区域 RMSD 和最大位移 |
| 端到端 | 真值掩码与预测掩码的成功率差距 |

掩码标签不唯一, 因此最终以修复成功、无新问题、固定区域保持和修改区域大小为主. 精确 mask IoU 不是唯一目标.

### 7.2 进入下一阶段的门槛

进入 Fixed Mask Head 训练前:

- 固定区域接近零位移.
- 离散化学图 100% 保持.
- prior 不读取 `G_good` 修复区域.
- 训练和推理的 prior 与时间构造一致.
- 真值掩码明显优于不修复和随机掩码.

进入交替联合适配前:

- Fixed Mask Head 能区分真实 clash 和无 clash 扭转对照.
- 正常样本误触发率处于可接受范围.
- 预测掩码与真值掩码的端到端差距已经量化.
- 共享层更新不降低真值掩码修复能力.

## 技术附录

<details>
<summary><strong>附录 A: 局部坐标 Flow 伪代码与时间语义</strong></summary>

修复训练使用 `(protein, G_bad, G_good, M_fixed_gt)`:

```python
fixed = ligand_mask & M_fixed
editable = ligand_mask & ~M_fixed

coords_0 = coords_bad.clone()
coords_0[editable] = sample_local_prior(
    reference_coords=coords_bad,
    editable=editable,
)

coords_t = coords_bad.clone()
coords_t[editable] = interpolate(
    coords_0[editable],
    coords_good[editable],
    time=t,
)

coord_times = t[:, None].expand_as(ligand_mask).float().clone()
coord_times[fixed] = 1.0
discrete_times = torch.ones_like(coord_times)
```

含义:

- 固定区域始终使用 `G_bad` 坐标, 并要求它与 `G_good` 相同.
- 修复区域从局部随机 prior 逐渐流向 `G_good`.
- 固定坐标已经完成, 所以坐标时间为 1.
- 所有离散属性已知且固定, 所以离散时间全部为 1.
- padding 由 `ligand_mask` 排除, 其时间值不参与模型计算和 loss.

`sample_local_prior()` 和 `interpolate()` 是概念性辅助函数, 实现时需映射到 FLOWR.ROOT interpolant.

prior 的中心、协方差和对齐只能读取 `G_bad`、固定区域和 pocket. 不使用会改变原子身份的 permutation alignment. 如做 SE(3) 增强, 必须对蛋白质、`G_bad` 和 `G_good` 应用同一刚体变换.

</details>

<details>
<summary><strong>附录 B: 样本记录示例</strong></summary>

```json
{
  "sample_id": "system_0001_corruption_03",
  "source_complex_id": "system_0001",
  "protein_path": "protein.pdb",
  "good_ligand_path": "good.sdf",
  "bad_ligand_path": "bad.sdf",
  "repair_atom_ids": [5, 6, 7, 8],
  "direct_clash_atom_ids": [7, 8],
  "rotatable_bond_atom_ids": [4, 5],
  "torsion_angle_deg": 120.0,
  "corruption_type": "torsion_clash",
  "has_target_clash": true,
  "diagnosis_label_kind": "clash_positive",
  "seed": 42
}
```

无 clash 扭转对照保存相同的旋转键、移动片段、角度和 seed, 并将 `has_target_clash` 设为 `false`. 所有样本保留来源复合物 ID, 用于审计数据划分和评价分母.

</details>

<details>
<summary><strong>附录 C: 后续候选模块</strong></summary>

这些模块不属于首版真实架构:

| 模块 | 何时考虑 | 作用 |
|---|---|---|
| Atom Error Type Head | 至少加入两种错误后 | 解释原子级错误类型 |
| Edge Error Head | 进入拓扑修复后 | 判断需要修改的键 |
| Repairability Head | 修复器冻结并获得 rollout 标签后 | 预测给定区域是否可被当前修复器修好 |

Atom Error Type Head 应优先采用多标签 sigmoid. Repairability 标签必须来自冻结修复器的重复 rollout 和外部验证, 不能用“掩码是否正确”代替.

</details>

<details>
<summary><strong>附录 D: 代码落地与测试清单</strong></summary>

Fixed Mask Head 的最小形式:

```python
self.fixed_mask_proj = torch.nn.Linear(d_inv, 1)
fixed_logits = self.fixed_mask_proj(invs_norm).squeeze(-1)
fixed_logits = fixed_logits.masked_fill(~ligand_mask, -1e9)

fixed_probability = torch.sigmoid(fixed_logits)
M_fixed = (fixed_probability >= fixed_threshold) & ligand_mask.bool()
```

建议新增:

```text
flowr/repair/
  corruption.py
  datasets.py
  interpolate.py
  masks.py
  validator.py
  oracle_runner.py
  records.py

scripts/repair/
  build_corruption_dataset.py
  run_oracle_repair.py
  evaluate_repair.py

tests/repair/
  test_masks.py
  test_forward_contracts.py
  test_atom_mapping.py
  test_coordinate_clamp.py
  test_corruption.py
  test_repair_interpolation.py
  test_training_routes.py
  test_no_target_leakage.py
  test_validator.py
```

主要修改点:

- `flowr/repair/interpolate.py`: coordinate-only prior、`G_t` 和两类时间.
- `flowr/repair/datasets.py`: diagnosis 和 repair 两种 batch.
- `flowr/models/mol_builder.py`: coordinate-only 固定状态和硬覆盖.
- `flowr/models/fm_pocket.py`: `diagnose()`、`repair_step()`、self-conditioning 和最终约束.
- `flowr/models/pocket.py`: Fixed Mask Head.
- `flowr/models/losses.py`: mask loss 和局部坐标 loss.
- `flowr/gen/generate.py`: 修复参数与整轮固定掩码.
- `flowr/gen/generate_from_pdb.py`: 修复命令行入口.

必须测试:

- `M_fixed` 与 `fragment_mask` 方向一致, padding 清零.
- 全固定掩码绕过.
- diagnosis 不输入 `M_fixed_gt`.
- repair prior 不读取 `G_good` 修复区域.
- 固定坐标时间为 1, 修复坐标时间为 `t`, 离散时间全部为 1.
- 局部坐标 loss 只覆盖修复原子并正确归一化.
- 单轮内 `M_fixed` 不变, 修复前向忽略 Fixed Mask Head.
- 每步和最终硬覆盖, 离散属性完全不变.
- SDF 往返读写后原子映射不变.
- validator 区分 clash、原始正常结构和无 clash 扭转对照.

推荐实施顺序:

1. 跑通官方 PDB 条件生成和手工 `substructure_inpainting`.
2. 完成数据划分、`M_fixed`、原子映射和验证器.
3. 生成 torsion clash 和无 clash 扭转对照.
4. 实现只依赖 `G_bad` 的 coordinate-only prior、插值和时间条件.
5. 实现每步及最终硬覆盖.
6. 在少量样本上跑通冻结修复器的真值掩码闭环.
7. 扩展样本并估计真值掩码上限.
8. 必要时训练 LoRA 修复器.
9. 增加并训练 Fixed Mask Head.
10. 跑通 predicted-mask 单轮修复并校准阈值.
11. 加入掩码鲁棒性, 最后评估交替训练和多轮修复.

</details>

<details>
<summary><strong>附录 E: 基线、消融和主要风险</strong></summary>

必要基线:

- 不修复.
- RDKit 局部优化.
- FLOWR.ROOT + 真值掩码.
- FLOWR.ROOT + 随机同大小修复区域.
- FLOWR.ROOT + 直接 clash 原子区域.
- FLOWR.ROOT + 预测掩码.

关键消融:

- 冻结修复器与 LoRA.
- 分阶段训练与后期交替适配.
- self-conditioning 开启与关闭.
- 直接 clash 标签与完整动作区域标签.
- 各向同性 prior、各向异性 prior 与 `G_bad + noise`.
- 不同掩码阈值和区域大小.
- 单连通与多连通区域.
- 单轮与多轮.
- 最终坐标 MSE 与 velocity loss.
- 硬覆盖与部分 Flow Matching.

主要风险与控制:

| 风险 | 控制方式 |
|---|---|
| 真值掩码也修不好 | 先解决数据或修复器, 不训练诊断头 |
| 模型只识别扭转 | 加入无 clash 扭转困难负样本 |
| prior 泄漏 `G_good` | prior 只允许读取 `G_bad`、固定区域和 pocket |
| 原子编号错位 | 保存映射并做 SDF round-trip 测试 |
| 早期预测掩码污染修复器 | 先用真值掩码, 预测掩码 detach 后逐步混入 |
| 两种 loss 梯度冲突 | 首版分阶段训练, 后期只做可回滚的交替适配 |
| 正常配体被误修 | 加入正常与困难负样本并校准阈值 |
| 多轮震荡 | 只接受改善结果, 保存历史最佳并限制轮数 |
| 硬覆盖边界不自然 | 监控锚点几何, 后续比较部分 Flow Matching |

</details>

## 参考资料

1. Lipman et al., *Flow Matching for Generative Modeling*, 2023.
2. Cremer et al., *FLOWR*, 2026.
3. Cremer et al., *FLOWR.ROOT*, 2026.
4. Dunn and Koes, *FlowMol3*, 2026.
5. Irwin et al., *SemlaFlow*, 2025.
