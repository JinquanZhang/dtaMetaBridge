\name{plot_sensspec_forest_meta}
\encoding{UTF-8}
\alias{plot_sensspec_forest_meta}
\title{从 meta 对象绘制灵敏度与特异度双森林图}
\description{研究行、四格表、区间、汇总菱形和异质性信息并排显示；通过 \code{...} 调整布局。}
\usage{
plot_sensspec_forest_meta(sensitivity_meta, specificity_meta, ...,
  use_meta_summary = TRUE)
}
\arguments{
\item{sensitivity_meta}{灵敏度的 \code{meta::metaprop} 对象，event/n 为 TP/(TP+FN)。}
\item{specificity_meta}{特异度的配对对象，event/n 为 TN/(TN+FP)，按研究名对齐。}
\item{use_meta_summary}{默认 TRUE，菱形采用 meta 对象的随机效应汇总值；FALSE 使用合并四格表的粗比例，不是双变量模型汇总值。}
\item{...}{传给 \code{plot_sensspec_forest} 的绘图参数，完整列表见下节。不要重复传入由适配器生成的 data、summary_override、heterogeneity 或 study_weights。}
}
\section{全部可调绘图参数}{
\describe{
\item{output_file = NULL}{NULL 在当前设备绘图；文件输出支持 PNG、TIFF、TIF。}
\item{study_col = "study"}{重建数据的研究名列，通常保持默认。}
\item{sens_axis = seq(0, 1, 0.2)}{灵敏度面板刻度；最小值和最大值同时确定横轴范围。}
\item{spec_axis = seq(0, 1, 0.2)}{特异度面板刻度和范围。范围外的区间目前没有截断箭头，建议保留完整 0--1 范围。}
\item{column_widths}{正数命名向量，默认 c(study=2.8, tp=0.5, fp=0.5, fn=0.5, tn=0.5, sens_text=1.6, spec_text=1.6, sens_plot=1.2, spec_plot=1.2)，表示九列相对宽度。}
\item{width = 10, height = NULL, res = 300}{保存图像的宽、高（英寸）和分辨率（dpi）；height 为 NULL 时自动计算。}
\item{font_family = "serif"}{所有文字统一字体；可用 "sans"、"mono" 或设备支持的字体。}
\item{header_cex = 0.86}{表头相对字号。}
\item{study_cex = 0.72}{研究名称及四格表计数相对字号。}
\item{ci_text_cex = 0.68}{研究效应值及区间相对字号。}
\item{summary_cex = 0.78}{汇总行标签相对字号。}
\item{summary_count_cex = 0.75}{汇总行计数相对字号。}
\item{summary_ci_cex = 0.7}{汇总效应值及区间相对字号。}
\item{axis_cex = 0.68}{横轴刻度相对字号。以上 cex 均为正数倍数，不是磅值。}
\item{row_gap = NULL}{研究行距，占画布高度的比例；默认 min(0.055, 0.52/max(1, 研究数-1))，总跨度不得超过 0.52。}
\item{square_col = "grey70"}{研究方块颜色。}
\item{square_max_mm = 5}{每个面板最大方块边长，毫米。}
\item{ci_col = "black", ci_lwd = 1.2}{研究 CI 横线与点估计短竖线的颜色和线宽。}
\item{diamond_col = "#C00000"}{汇总菱形颜色。}
\item{heterogeneity_cex = 0.66}{异质性两行文字的相对字号。}
\item{heterogeneity_x = 0.015}{两行文字左端位置，0 为画布左侧、1 为右侧。}
\item{heterogeneity_y = NULL}{c(灵敏度, 特异度) 两行纵坐标，0 为底部、1 为顶部；NULL 自动置于汇总行下。自定义坐标不自动避让其他图形。}
}
}
\details{
方块面积分别与两个 meta 对象的 \code{w.random} 成正比；缺失或全部 NA 时警告并画等大方块，不改变模型估计。
先画方块，再画完整 CI 和点估计短竖线。
研究区间使用 Clopper--Pearson 精确二项区间，可能不同于输入 meta 对象所选的区间算法。
异质性来自对应 meta 对象，与 SROC 的双变量模型不是同一个模型；汇总灵敏度和特异度也不保证相同。
}
\value{绘制图形并不可见地返回输出文件名；当前设备绘图时返回 NULL。}
\seealso{\code{\link{plot_sensspec_forest}}, \code{\link{plot_sroc}}}
\examples{
d <- data.frame(study = c("A", "B", "C"),
  TP = c(8, 12, 16), FN = c(2, 4, 3),
  TN = c(17, 18, 20), FP = c(3, 2, 1))
se <- meta::metaprop(TP, TP + FN, studlab = study, data = d,
  method = "Inverse")
sp <- meta::metaprop(TN, TN + FP, studlab = study, data = d,
  method = "Inverse")
plot_sensspec_forest_meta(se, sp, font_family = "sans",
  header_cex = 0.9, study_cex = 0.75, square_col = "grey70",
  heterogeneity_cex = 0.65)
}
