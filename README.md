# dtaMetaBridge

`dtaMetaBridge` 用于衔接配对的 `meta::metaprop()` 对象与标准诊断试验准确性
meta 分析，提供双森林图、双变量二项 GLMM、原始 QMD 的 SROC 函数和检验后概率。
0.3.0 默认使用 `HFphf_meta analysis.qmd` 中的 `fit_metandi()` / `plot_metandi()`。

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
  study = c("Parcha 2021", "Tada 2021", "Forsyth 2021", "Reddy 2022",
            "Ariyaratnam 2024", "Akerman 2025", "Wang 2023", "Nguyen 2025", "Rahi 2026"),
  TP = c(169,113,11,163,61,117,61,110,71), FP = c(29,5,16,3,11,42,2,36,27),
  FN = c(122,113,39,222,27,123,140,8,9), TN = c(350,173,11,99,21,214,116,126,85)
)
```

以上为原 QMD 的 H2 rule-in 输入，用于复现文件输出，不表示已重新核验各论文的四格表。

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

# 3. 默认采用原始 QMD 的拟合和绘图函数。
fit <- fit_bivariate_meta(meta_sens, meta_spec)
plot_sroc(fit)

# 也可以直接使用 QMD 中的原函数名：
# fit <- fit_metandi(dta)
# plot_metandi(fit)

# 4. 按预检概率计算阳性和阴性后的患病概率及 95% 不确定性区间。
posttest_probability(fit, prevalence = c(0.10, 0.30, 0.50))
```

### 示例图

下图由上面的示例数据和函数生成。左图的森林图汇总值来自两个独立的
`meta::metaprop()` 随机效应模型；SROC 的方块来自联合二项 GLMM，故两组汇总
灵敏度/特异度数值可能略有差异，这是模型定义不同所致，并非计算不一致。

下表左侧为生成代码，右侧为对应输出图。代码使用上文的 `dta`、`meta_sens`、
`meta_spec` 和 `fit` 对象。

<table>
<tr><th>绘图代码</th><th>示例图</th></tr>
<tr>
<td><pre><code>plot_sensspec_forest_meta(
  meta_sens, meta_spec,
  output_file = "forest-example.png",
  width = 10, res = 300
)</code></pre></td>
<td><img src="inst/figures/forest-example.png" alt="双森林图" width="600"></td>
</tr>
<tr>
<td><pre><code>plot_sroc(
  fit,
  output_file = "sroc-example.png",
  width = 6, height = 6, dpi = 300
)</code></pre></td>
<td><img src="inst/figures/sroc-example.png" alt="SROC 曲线" width="450"></td>
</tr>
</table>

### 异质性与列宽调整

使用 `plot_sensspec_forest_meta()` 时，图底部会自动显示灵敏度和特异度各自的
`I2`、`tau2`、Cochran Q 和 p 值。研究方块按相应的患病者或非患病者样本量缩放，
汇总菱形为红色。`column_widths` 是一个命名数值向量；数值是各列的相对宽度，可按
版面需要调整。下例加宽研究名称和两张森林图区，并指定与示例图相同的坐标刻度：

```r
plot_sensspec_forest_meta(
  meta_sens, meta_spec,
  output_file = "forest-wide.png",
  sens_axis = seq(0.2, 1, 0.2),
  spec_axis = seq(0.4, 1, 0.2),
  column_widths = c(
    study = 3.2, tp = 0.45, fp = 0.45, fn = 0.45, tn = 0.45,
    sens_text = 1.5, spec_text = 1.5,
    sens_plot = 1.2, spec_plot = 1.2
  )
)
```

### 结果如何解释

- 森林图菱形：分别汇总灵敏度与特异度，适用于展示每个结局的异质性。
- SROC：默认 `sroc_type = "qmd"`，使用所提供 QMD 的原始函数，保留研究编号、
  图例、汇总值与 AUC 模拟区间。`seed`、`n_mc` 和 `n_grid` 控制随机种子、模拟次数及网格。
  `midas`、`ruttergatsonis` 和 `naive` 仍可显式选择。
  `method`、`correction` 和 `correction_control` 仅用于两个 `mada` 后端。
- `fit$metrics$auc`：跨研究的 SROC 区分能力汇总，不能替代单项研究中连续评分的 ROC AUC。
- post-test probability：区间反映汇总平均准确性的抽样不确定性，不是未来任一新场景的预测区间。

## 主要函数

| 函数 | 用途 |
| --- | --- |
| `dta_from_meta()` | 从匹配的灵敏度和特异度 `metaprop` 对象还原四格表。 |
| `plot_sensspec_forest_meta()` | 保留 `meta` 随机效应汇总菱形，绘制双森林图。 |
| `fit_bivariate_meta()` / `fit_bivariate_dta()` | 从 `meta` 对象或四格表拟合双变量模型；默认原始 QMD 算法。 |
| `fit_metandi()` / `plot_metandi()` | 直接使用 QMD 中的原函数及其图例参数。 |
| `plot_sroc()` | 绘制 HSROC 曲线、置信轮廓、预测轮廓和研究点。 |
| `posttest_probability()` | 计算不同预检概率下的 PPV、NPV 及模拟法 95% 区间。 |

## 使用范围

### 原始 QMD 方法（默认，0.3.0）

原始算法保持不变：在 logit 特异度上，曲线斜率为
`sign(rho_re) * sigma_se / sigma_sp`。这不是 `mada` 的 naive 条件均值斜率，
也不是本包可选 Midas 后端的固定负斜率。相关性为正时，原函数仍可能画出反向曲线；
移植时不会自动更改其符号。
AUC 点估计为固定效应模拟所得面积的均值，区间未包含研究间方差参数的不确定性。
原函数执行 `set.seed(seed)`，会改变全局随机数状态。

```r
fit <- fit_metandi(dta, seed = 2026, n_mc = 3000, n_grid = 1000)
plot_metandi(fit, show_conf = TRUE, show_pred = TRUE, show_legend = TRUE)
# 原函数绘图样式也支持直接保存：
plot_sroc(fit, output_file = "sroc-qmd.png", width = 6, height = 6)
```

### Midas 可选方法与复现边界（`sroc_type = "midas"`）

依据本机 `midas.ado` 2.00（2008-12-21）的 SUMMARY ROC CURVE 段，令
`mu_se`、`mu_sp` 为两项 logit 均值，`v_se`、`v_sp` 为对应研究间方差：

```r
b <- (max(0.001, v_sp) / max(0.001, v_se))^0.25
a <- mu_se * b + mu_sp / b
sensitivity <- plogis((a - qlogis(specificity) / b) / b)
```

该曲线经过联合汇总点，以 500 个等距点在特异度 0–1 上作梯形积分得到 AUC。
这里的 Midas 公式与 `mada` 的 `naive` 条件均值曲线不同；此前将两者等同的说明已更正。
R 端使用 `lme4::glmer(nAGQ = 1)` 的 Laplace 估计；Midas 的自适应积分配置可能得到不同估计。
本次验证覆盖公式、端点、单调性、汇总点及数值积分，尚未完成 Stata 实际运行的数值对照，
因此不能宣称结果逐位一致。边界或不收敛信息见 `fit$diagnostics`。
旧版 Midas 使用绘图点数构造 AUC 区间，本包不把它作为统计不确定性区间，AUC 上下限保留 `NA`。
`pauc` 为观察到的特异度范围内的未标准化部分面积。

本包用于可重复的诊断试验准确性 meta 分析。在合并前仍需审查研究定义、阈值及临床
适用性；本包不能替代研究方案或偏倚风险评价。
