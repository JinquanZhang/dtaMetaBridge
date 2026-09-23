\name{posttest_probability}
\encoding{UTF-8}
\alias{posttest_probability}
\title{计算检验后概率：阳性预测值与阴性预测值}
\description{结合双变量模型汇总灵敏度、特异度与指定检验前概率，用贝叶斯公式计算 PPV、NPV。}
\usage{
posttest_probability(fit, prevalence = 0.3, n_sims = 3000, seed = 2026)
}
\arguments{
\item{fit}{\code{fit_bivariate_meta} 或 \code{fit_bivariate_dta} 的结果；底层模型为 reitsma 或 glmerMod。}
\item{prevalence}{一个或多个检验前患病概率，必须严格介于 0 与 1；不是百分数。}
\item{n_sims}{区间计算的蒙特卡洛模拟次数，不小于 100 的整数；默认 3000。}
\item{seed}{随机种子，默认 2026；函数结束后恢复原随机数状态。}
}
\details{
PPV 是阳性后的患病概率；NPV 是阴性后的未患病概率。
阴性后的患病概率为 1-NPV，其区间端点按 1-npv_upr、1-npv_lwr 计算。
模拟从模型固定效应及其协方差的多元正态近似分布抽样，使用 2.5\% 与 97.5\% 分位数。
区间不包含检验前概率的不确定性，也不是新研究预测区间。
本函数不是重新计算“后验灵敏度、特异度”；返回的 sensitivity、specificity 是输入模型的汇总值。
}
\value{
数据框，每个检验前概率一行：
prevalence（检验前概率）、sensitivity（汇总灵敏度）、
specificity（汇总特异度）、ppv、ppv_lwr、ppv_upr、
npv、npv_lwr、npv_upr（预测值及其模拟 95\% 区间）。
}
\seealso{\code{\link{fit_bivariate_meta}}, \code{\link{plot_sroc}}}
\examples{
d <- data.frame(study = LETTERS[1:6],
  TP = c(35, 42, 28, 60, 45, 70), FN = c(15, 8, 22, 20, 15, 10),
  TN = c(80, 65, 90, 55, 75, 60), FP = c(20, 35, 10, 45, 25, 40))
fit <- fit_bivariate_dta(d)
posttest_probability(fit, prevalence = c(0.1, 0.3, 0.5),
  n_sims = 1000, seed = 2026)
}
