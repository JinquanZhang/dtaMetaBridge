\name{posttest_probability}
\alias{posttest_probability}
\title{计算诊断检验后的阳性和阴性预测概率}
\description{
根据 \code{\link{fit_sroc}} 的汇总灵敏度、特异度及其协方差，计算给定
检验前概率下的阳性预测值（PPV）和阴性预测值（NPV）及不确定性区间。
}
\usage{
posttest_probability(fit, prevalence = 0.3, n_sims = 3000L, seed = 2026L)
}
\arguments{
  \item{fit}{由 \code{fit_sroc()} 返回的对象。}
  \item{prevalence}{单个数值或数值向量，表示检验前概率，必须位于 0 与 1 之间。}
  \item{n_sims}{模拟次数，默认 3000。}
  \item{seed}{随机数种子。}
}
\value{
数据框，每行对应一个 \code{prevalence}，含 \code{sensitivity}、
\code{specificity}、\code{ppv}、\code{ppv_lower}、\code{ppv_upper}、
\code{npv}、\code{npv_lower}、\code{npv_upper}。阴性检验后的患病概率为
\code{1 - npv}。
}
\details{
函数在 logit(灵敏度)--logit(特异度) 尺度上，按拟合模型的固定效应均值和协方差
矩阵进行正态模拟，再用贝叶斯公式计算 PPV 与 NPV。区间仅反映汇总灵敏度和特异度
的不确定性；不包含检验前概率的不确定性，也不是新研究的预测区间。频率学与贝叶斯
\code{fit_sroc()} 结果均可使用。
}
\examples{
\dontrun{
fit <- fit_sroc(dat, backend = "frequency", auc_boot = 200)
posttest_probability(fit, prevalence = c(0.1, 0.3, 0.5))
}
}
\keyword{models}
