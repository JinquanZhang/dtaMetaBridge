\name{fit_sroc}
\alias{fit_sroc}
\title{拟合诊断试验准确性的 SROC 模型}
\description{
对四格表诊断试验数据拟合直接频率学 Rutter--Gatsonis HSROC 模型，
或贝叶斯双变量模型，并计算汇总灵敏度、汇总特异度、AUC 及其区间。
}
\usage{
fit_sroc(data, backend = c("frequency", "bayes"), sroc_type = 5L,
  study_col = "study", year_col = "Year", conf_level = 0.95,
  n_grid = 1000L, auc_boot = 2000L, n_cores = 6L,
  posterior_samples = 2000L, seed = 2026L)
}
\arguments{
  \item{data}{含 TP、FP、FN、TN 四列的数据框；各列须为非负整数。}
  \item{backend}{\code{"frequency"}（默认）为直接 Rutter--Gatsonis
  频率学拟合；\code{"bayes"} 为 \pkg{meta4diag} 的贝叶斯双变量拟合。}
  \item{sroc_type}{SROC 公式编号。频率学后端仅支持 5（Rutter--Gatsonis）；
  贝叶斯后端支持 1--5。}
  \item{study_col}{研究标签列名。}
  \item{year_col}{年份列名；不存在时忽略。}
  \item{conf_level}{区间置信度或可信度，默认 0.95。}
  \item{n_grid}{绘制/积分 SROC 曲线的网格数。}
  \item{auc_boot}{频率学 AUC 的研究层 Bootstrap 次数；设为 0 不计算区间。}
  \item{n_cores}{频率学 Bootstrap 使用的最大核心数，默认 6。}
  \item{posterior_samples}{贝叶斯后端用于绘图和汇总的后验抽样数。}
  \item{seed}{随机数种子。}
}
\value{
返回类 \code{dta_sroc_fit} 对象，包含原始 \code{input_data}、拟合
\code{model}、\code{metrics}（灵敏度、特异度、AUC、部分 AUC 及区间）、
\code{plot_data}、\code{curve_parameters} 和 \code{random_effects}。
频率学区间为 Wald/Delta 或 Bootstrap 置信区间；贝叶斯区间为后验可信区间。
区间水平由 \code{conf_level} 决定，并保存在 \code{conf_level} 和
\code{interval_type} 元素中。
}
\details{
频率学后端直接调用 \code{dtametaTMB::fitRutterGatsonis()}。其曲线为
\deqn{\operatorname{logit}(Se)=\Lambda\exp(-\beta/2)-
\exp(-\beta)\operatorname{logit}(Sp).}
AUC 为该曲线在完整特异度区间上的数值积分；其默认 95\% CI 由研究层
Bootstrap 获得。贝叶斯后端依赖 \pkg{meta4diag} 和 INLA，公式 1--5
供不同 SROC 参数化选择。

每行必须代表独立研究。对于同一队列的多个阈值，应使用多阈值模型，
而不是将重复阈值行作为独立研究输入。
}
\references{
Rutter CM, Gatsonis CA (2001). A hierarchical regression approach to
meta-analysis of diagnostic test accuracy evaluations. \emph{Statistics in
Medicine}, 20, 2865--2884.
}
\examples{
\dontrun{
dat <- data.frame(study = c("A", "B", "C"),
  TP = c(52, 15, 20), FP = c(1, 4, 1),
  FN = c(19, 5, 4), TN = c(16, 29, 19))
fit <- fit_sroc(dat, backend = "frequency", auc_boot = 200)
fit$metrics
}
}
\keyword{models}
