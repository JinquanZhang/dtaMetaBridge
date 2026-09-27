\name{fit_sroc}
\encoding{UTF-8}
\alias{fit_sroc}
\title{从四格表拟合双变量模型并选择五种 SROC 曲线}
\description{默认使用 \pkg{dtametaTMB} 的频率学 Reitsma 双变量模型；五种曲线公式与 \pkg{meta4diag} 的 \code{SROC()} 编号一致。}
\usage{
fit_sroc(data, backend = c("frequency", "bayes", "legacy"),
  sroc_type = 5L, study_col = "study", year_col = "Year",
  conf_level = 0.95, n_grid = 1000, auc_boot = 0L,
  posterior_samples = 2000L, seed = 2026)
}
\arguments{
\item{data}{数据框，必须有研究名称及 \code{TP}、\code{FP}、\code{FN}、\code{TN} 非负整数列。每行须代表独立研究。}
\item{backend}{\code{"frequency"}（默认）是基于二项似然的频率学双变量模型；\code{"bayes"} 是基于 INLA 的贝叶斯双变量模型；\code{"legacy"} 保留为旧的正态近似实现。}
\item{sroc_type}{1--5 的整数：1 回归线 1；2 主轴；3 Moses--Littenberg 型曲线；4 回归线 2；5 Rutter--Gatsonis（默认）。这是同一模型上的曲线公式选择，不是五种独立模型。}
\item{study_col}{研究名称列，默认 \code{study}。}
\item{year_col}{可选年份列，默认 \code{Year}；存在时附在显示标签后。设为 \code{NULL} 忽略年份。}
\item{conf_level}{灵敏度、特异度与轮廓的区间水平，默认 0.95。}
\item{n_grid}{观察到的 FPR 范围内曲线点数；为不小于 100 的整数。}
\item{auc_boot}{AUC 的研究层非参数 Bootstrap 次数。0（默认）只返回 AUC 点估计；通常建议正式分析设为 2000 或更多。}
\item{posterior_samples}{仅 \code{backend = "bayes"} 使用：用于非线性指标与 AUC 后验可信区间的抽样次数，默认 2000，至少为 20。}
\item{seed}{\code{auc_boot > 0} 时的随机种子；函数会恢复调用前的随机数状态。}
}
\details{
坐标统一为 \eqn{x=logit(FPR)}、\eqn{y=logit(sensitivity)}。令 \eqn{a,b,c} 分别为
灵敏度和 FPR 的随机效应方差、协方差，曲线写为 \eqn{y=mu_y+B(x-mu_x)}。
五个斜率依次为 \eqn{c/b}、主轴特征向量之比、\eqn{(a+c)/(b+c)}、\eqn{a/c} 和 \eqn{sqrt(a/b)}。
类型 5 在该参数化下对应 Rutter--Gatsonis HSROC 曲线。

\code{"frequency"} 直接拟合四格表的二项似然，返回汇总灵敏度、特异度及 Wald 95\% CI；指定 \code{auc_boot} 后，AUC 区间来自按研究重抽样、重新拟合模型的百分位 Bootstrap。
\code{"bayes"} 需要已安装并可加载 \pkg{INLA}，返回后验均值和 95\% CrI；AUC 的 CrI 由 \code{posterior_samples} 控制。
两种区间的概率解释不同，不能视为同一类不确定性区间。
\code{plot_sroc(fit, full_curve = TRUE)} 可画 FPR 0--1 的完整模型曲线；观察范围外为外推。

同一受试者队列的多个阈值不能按独立行送入本函数，即使研究名称不同也不行。请使用多阈值模型
（例如 \code{dtametaTMB::fitHoyer()} 或 \pkg{diagmeta}）保留阈值相关性。
}
\value{与 \code{plot_sroc()} 和 \code{posttest_probability()} 兼容的列表，包含拟合对象 \code{model}、\code{metrics}、\code{plot_data}、\code{curve_parameters} 及未修改的 \code{input_data}。}
\seealso{\code{\link{plot_sroc}}, \code{\link{posttest_probability}},
\url{https://rdrr.io/cran/dtametaTMB/man/fitReitsma.html},
\url{https://search.r-project.org/CRAN/refmans/meta4diag/html/AUC.html}}
\references{Arends, L. R. et al. (2008). Bivariate random effects meta-analysis of ROC curves. \emph{Medical Decision Making}, 28, 621--638.}
\examples{
\dontrun{
data("RF", package = "dtametaTMB")
fit <- fit_sroc(RF, backend = "frequency", year_col = "year", sroc_type = 5, auc_boot = 2000)
fit$metrics
plot_sroc(fit, full_curve = TRUE)

# 贝叶斯模型（需先安装并加载 INLA）
fit_bayes <- fit_sroc(RF, backend = "bayes", year_col = "year",
  sroc_type = 5, posterior_samples = 2000)
}
}
