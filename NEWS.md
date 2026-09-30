# dtaMetaBridge 1.6.0

- SROC 图例改用 `grid` 实测每一行文字的渲染宽度，并按目标输出宽高、字号、符号列和内边距计算自适应边框；低尺寸无法容纳完整文字时明确报错，不再截断。

# dtaMetaBridge 1.5.9

- 清理横轴切换实现的格式检查警告；功能不变。

# dtaMetaBridge 1.5.8

- `plot_sroc()` 新增 `x_axis`：`"specificity"`（默认）或 `"fpr"`，分别显示 Specificity 与 1 - Specificity 横轴；图例位置按屏幕坐标定位，不受横轴方向影响。

# dtaMetaBridge 1.5.7

- SROC 图例全部使用常规字重。

# dtaMetaBridge 1.5.6

- 清理 SROC 图例实现的格式检查警告；功能不变。

# dtaMetaBridge 1.5.5

- SROC 图例框宽度按最长文字行自适应，高度按文字行数、字号和显示的轮廓条目自适应。

# dtaMetaBridge 1.5.4

- SROC 研究圈移至预测/置信轮廓及曲线下方，避免遮挡线条；新增 `study_label_col`，默认研究编号为黑色。

# dtaMetaBridge 1.5.3

- `plot_sroc()` 的汇总点改为 45° 旋转正方形（菱形）；新增 `study_col` 控制研究圈及编号颜色，`sroc_col` 仅控制 SROC 曲线颜色。

# dtaMetaBridge 1.5.2

- `fit_sroc()` 新增 `n_cores`（默认 6），以 R 标准库 `parallel` 并行执行频率学 AUC Bootstrap，并使用独立随机数流保证可复现；自动不超过可用物理核心数。

# dtaMetaBridge 1.5.1

- 频率学 `fit_sroc()` 现在默认执行 2,000 次研究层 Bootstrap，并返回 Rutter--Gatsonis 曲线 AUC 的 95% CI；需要快速预览时可设 `auc_boot = 0`。

# dtaMetaBridge 1.5.0

- `backend = "frequency"` 改为直接拟合 `dtametaTMB::fitRutterGatsonis()`；SROC 与 AUC 直接使用 Rutter--Gatsonis 的 \(\Lambda\) 与 \(\beta\) 参数，不再由 Reitsma 拟合结果转换。
- 频率学后端仅支持原生的 Rutter--Gatsonis 曲线（`sroc_type = 5`）；贝叶斯后端仍可选择 1--5 种曲线公式。
- 灵敏度、特异度和检验后概率的频率学不确定性由同一直接模型的参数协方差经 Delta 法计算。

# dtaMetaBridge 1.4.0

- 双森林图支持统一字体、表头与研究行字号、行距、列宽、独立横轴、颜色，以及异质性文字大小和位置。
- 灰色方块面积按各面板的 meta 随机效应权重缩放；缺失权重时警告并使用等大方块。完整 CI 画在方块上层，点估计用短竖线标记。
- SROC 默认使用 Reitsma REML 模型与 Rutter–Gatsonis 曲线；支持完整曲线外推、图例位置、字体、颜色、透明度及刻度设置。
- 后验概率支持多个检验前概率，并计算 PPV、NPV 及模拟区间。
- 为三个主要函数提供独立中文 R 帮助页和示例。
- 保留 QMD 与 Midas 风格可选方法；不声称与 Stata 数值完全一致，也不将输入数据视为已核验的论文原始数据。
