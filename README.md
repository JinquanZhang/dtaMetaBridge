# dtaMetaBridge

`dtaMetaBridge` 用于衔接配对的 `meta::metaprop()` 对象与标准诊断试验准确性
meta 分析，提供灵敏度/特异度双森林图、Reitsma 双变量模型、Rutter-Gatsonis
HSROC 曲线及依赖预检概率的 post-test probability。

## 使用教程

### 安装

```r
install.packages("remotes")       # 只需安装一次
remotes::install_github("JinquanZhang/dtaMetaBridge")
library(dtaMetaBridge)
```

### 数据格式

每一行是一项研究，必须有研究名称及四格表计数：`TP`（真阳性）、`FP`（假阳性）、
`FN`（假阴性）、`TN`（真阴性）。

```r
dta <- data.frame(
  study = c("Study A", "Study B", "Study C", "Study D", "Study E"),
  TP = c(40, 32, 48, 25, 55), FP = c(8, 15, 10, 20, 12),
  FN = c(10, 18, 7, 15, 9),  TN = c(72, 65, 80, 58, 90)
)
```

### 从 `meta` 对象开始：森林图、SROC 与 post-test probability

```r
library(meta)

# 1. 分别拟合灵敏度和特异度的单变量随机效应 GLMM。
meta_sens <- metaprop(TP, TP + FN, studlab = study, data = dta,
                      method = "GLMM", method.tau = "ML")
meta_spec <- metaprop(TN, TN + FP, studlab = study, data = dta,
                      method = "GLMM", method.tau = "ML")

# 2. 双森林图：菱形直接采用上述 meta 对象的随机效应汇总值。
plot_sensspec_forest_meta(meta_sens, meta_spec,
                           output_file = "forest.png")

# 3. 联合 Reitsma 模型与 Rutter-Gatsonis HSROC。
fit <- fit_bivariate_meta(meta_sens, meta_spec)
plot_sroc(fit)

# 4. 按预检概率计算阳性和阴性后的患病概率及 95% 不确定性区间。
posttest_probability(fit, prevalence = c(0.10, 0.30, 0.50))
```

### 示例图

下图由上面的示例数据和函数生成。左图的森林图汇总值来自两个独立的
`meta::metaprop()` 随机效应模型；右图的方块则来自联合 Reitsma 模型，故两组汇总
灵敏度/特异度数值可能略有差异，这是模型定义不同所致，并非计算不一致。

| 双森林图 | SROC 曲线 |
| --- | --- |
| ![双森林图](inst/figures/forest-example.png) | ![SROC 曲线](inst/figures/sroc-example.png) |

### 结果如何解释

- 森林图菱形：分别汇总灵敏度与特异度，适用于展示每个结局的异质性。
- SROC：使用 `mada::reitsma()` 的标准 Reitsma 双变量随机效应模型和
  Rutter-Gatsonis 参数化；曲线仅在研究实际观察到的假阳性率范围内绘制，以避免端点外推。
- `fit$metrics$auc`：HSROC 的跨研究区分能力汇总，不能替代单项研究中连续评分的 ROC AUC。
- post-test probability：区间反映汇总平均准确性的抽样不确定性，不是未来任一新场景的预测区间。

## 主要函数

| 函数 | 用途 |
| --- | --- |
| `dta_from_meta()` | 从匹配的灵敏度和特异度 `metaprop` 对象还原四格表。 |
| `plot_sensspec_forest_meta()` | 保留 `meta` 随机效应汇总菱形，绘制双森林图。 |
| `fit_bivariate_meta()` / `fit_bivariate_dta()` | 从 `meta` 对象或四格表拟合标准 Reitsma 双变量模型。 |
| `plot_sroc()` | 绘制 HSROC 曲线、置信轮廓、预测轮廓和研究点。 |
| `posttest_probability()` | 计算不同预检概率下的 PPV、NPV 及模拟法 95% 区间。 |

## 使用范围

本包用于可重复的诊断试验准确性 meta 分析。在合并前仍需审查研究定义、阈值及临床
适用性；本包不能替代研究方案或偏倚风险评价。
