# dtaMetaBridge

[中文](#中文教程) | [English](#english-tutorial)

`dtaMetaBridge` bridges paired `meta::metaprop()` objects and standard diagnostic
test-accuracy meta-analysis. It provides a paired sensitivity/specificity forest
plot, a Reitsma bivariate model with a Rutter-Gatsonis HSROC curve, and
prevalence-dependent post-test probabilities.

## 中文教程

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

| 双森林图 / Paired forest plot | SROC 曲线 / SROC curve |
| --- | --- |
| ![Paired forest plot](inst/figures/forest-example.png) | ![SROC plot](inst/figures/sroc-example.png) |

### 结果如何解释

- 森林图菱形：分别汇总灵敏度与特异度，适用于展示每个结局的异质性。
- SROC：使用 `mada::reitsma()` 的标准 Reitsma 双变量随机效应模型和
  Rutter-Gatsonis 参数化；曲线仅在研究实际观察到的假阳性率范围内绘制，以避免端点外推。
- `fit$metrics$auc`：HSROC 的跨研究区分能力汇总，不能替代单项研究中连续评分的 ROC AUC。
- post-test probability：区间反映汇总平均准确性的抽样不确定性，不是未来任一新场景的预测区间。

## English tutorial

### Install

```r
install.packages("remotes")       # once only
remotes::install_github("JinquanZhang/dtaMetaBridge")
library(dtaMetaBridge)
```

### Input data

Use one row per study with a study label and the 2x2 counts: `TP`, `FP`, `FN`,
and `TN`. The `dta` object in the Chinese example above is a complete runnable
example.

### Workflow from `meta` objects

```r
library(meta)

meta_sens <- metaprop(TP, TP + FN, studlab = study, data = dta,
                      method = "GLMM", method.tau = "ML")
meta_spec <- metaprop(TN, TN + FP, studlab = study, data = dta,
                      method = "GLMM", method.tau = "ML")

plot_sensspec_forest_meta(meta_sens, meta_spec, output_file = "forest.png")

fit <- fit_bivariate_meta(meta_sens, meta_spec)
plot_sroc(fit)
posttest_probability(fit, prevalence = c(0.10, 0.30, 0.50))
```

### Key points

- Forest-plot diamonds are separate univariate random-effects summaries retained
  by `meta`; the bivariate SROC summary point can therefore differ slightly.
- The SROC is the standard `mada::reitsma()` Rutter-Gatsonis HSROC curve. It is
  shown only over the observed false-positive-rate range.
- Post-test probabilities depend on the setting-specific pre-test probability.
  Their intervals quantify uncertainty in pooled mean accuracy, not prediction
  intervals for a new setting.

## Functions

| Function | Purpose |
| --- | --- |
| `dta_from_meta()` | Recover a 2x2 data frame from matched sensitivity and specificity `metaprop` objects. |
| `plot_sensspec_forest_meta()` | Create the paired forest plot while retaining the `meta` random-effects diamonds. |
| `fit_bivariate_meta()` / `fit_bivariate_dta()` | Fit the standard Reitsma bivariate model from `meta` objects or a 2x2 data frame. |
| `plot_sroc()` | Draw the HSROC curve, confidence contour, prediction contour, and study points. |
| `posttest_probability()` | Calculate prevalence-specific PPV and NPV with simulation-based 95% intervals. |

## Citation and scope

This package is intended for reproducible diagnostic-test-accuracy meta-analysis.
Check study definitions, thresholds, and clinical applicability before pooling;
the package does not replace a protocol or a risk-of-bias assessment.
