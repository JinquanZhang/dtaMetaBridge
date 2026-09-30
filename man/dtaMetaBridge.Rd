\name{dtaMetaBridge}
\alias{dtaMetaBridge}
\alias{dta_from_meta}
\title{诊断试验准确性 Meta 分析工具}
\description{
\pkg{dtaMetaBridge} 用于诊断试验准确性（DTA）Meta 分析：从配对的
\pkg{meta} 对象还原四格表、绘制灵敏度与特异度双森林图、拟合 SROC，
以及计算后验（检验后）概率。
}
\usage{
dta_from_meta(sensitivity_meta, specificity_meta)
}
\arguments{
  \item{sensitivity_meta}{由 \code{meta::metaprop()} 创建的灵敏度对象；
  事件数为 TP，样本量为 TP + FN。}
  \item{specificity_meta}{由 \code{meta::metaprop()} 创建的特异度对象；
  事件数为 TN，样本量为 TN + FP。研究标签必须能与
  \code{sensitivity_meta} 一一对应。}
}
\value{
\code{dta_from_meta()} 返回按研究排列的 \code{data.frame}，包含
\code{study}、\code{TP}、\code{FP}、\code{FN}、\code{TN} 五列。
}
\details{
用户可用的主要函数为：
\itemize{
  \item \code{dta_from_meta()}：将一对 \pkg{meta} 比例 Meta 分析对象转换为四格表；
  \item \code{plot_sensspec_forest_meta()}：绘制配对灵敏度/特异度森林图；
  \item \code{fit_sroc()}：拟合频率学 Rutter--Gatsonis 或贝叶斯双变量 SROC 模型；
  \item \code{plot_sroc()}：绘制 SROC、置信轮廓和预测轮廓；
  \item \code{posttest_probability()}：以拟合模型计算 PPV 和 NPV 的不确定性区间。
}
这些函数假定每行代表独立研究。若同一受试者群在多个阈值下重复报告，
不应把各阈值当作相互独立的研究；应采用专门的多阈值模型。
}
\seealso{\code{\link{fit_sroc}}, \code{\link{plot_sroc}},
\code{\link{plot_sensspec_forest_meta}}, \code{\link{posttest_probability}}}
\examples{
\dontrun{
library(meta)
sens <- metaprop(event = c(52, 15), n = c(71, 20),
                 studlab = c("A", "B"))
spec <- metaprop(event = c(16, 29), n = c(17, 33),
                 studlab = c("A", "B"))
dta_from_meta(sens, spec)
}
}
\keyword{models}
\keyword{dplot}
