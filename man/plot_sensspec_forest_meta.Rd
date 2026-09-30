\name{plot_sensspec_forest_meta}
\alias{plot_sensspec_forest_meta}
\title{绘制灵敏度和特异度双森林图}
\description{
以配对的 \pkg{meta} 比例 Meta 分析对象，或直接以诊断四格表，绘制灵敏度和
特异度的并列森林图。
}
\usage{
plot_sensspec_forest_meta(sensitivity_meta, specificity_meta = NULL,
  ..., use_meta_summary = TRUE)
}
\arguments{
  \item{sensitivity_meta}{灵敏度 \code{meta::metaprop} 对象；或含
  \code{study}、\code{TP}、\code{FP}、\code{FN}、\code{TN} 的数据框。}
  \item{specificity_meta}{对应特异度 \code{meta::metaprop} 对象；
  当第一个参数为四格表数据框时保持 \code{NULL}。}
  \item{...}{绘图参数。常用参数包括：\code{output_file}、\code{width}、
  \code{height}、\code{res}；两图横轴范围 \code{sens_axis}、\code{spec_axis}；
  列宽 \code{column_widths}；字体 \code{font_family}、\code{header_cex}、
  \code{study_cex}、\code{ci_text_cex}、\code{summary_cex}、\code{axis_cex}；
  行距 \code{row_gap}；方块 \code{square_col}、\code{square_max_mm}；
  CI 线 \code{ci_col}、\code{ci_lwd}；汇总菱形 \code{diamond_col}；
  异质性文字 \code{heterogeneity_cex}、\code{heterogeneity_x}、
  \code{heterogeneity_y}。}
  \item{use_meta_summary}{当输入为配对 \pkg{meta} 对象时，是否使用其
  随机效应汇总值及权重。默认 \code{TRUE}。}
}
\value{
以 \code{invisible()} 返回一个含绘图使用数据和汇总结果的列表；指定
\code{output_file} 时保存图片，否则绘制到当前图形设备。
}
\details{
研究行的灵敏度与特异度均由同一四格表计算。单研究 CI 使用精确二项分布区间；
若输入 \pkg{meta} 对象且 \code{use_meta_summary = TRUE}，汇总菱形、权重和
异质性取自各自的随机效应 Meta 分析。灰色方块面积与随机效应权重成比例，先绘制
方块，再绘制完整 CI 横线；点估计位置使用短竖线而非加号。

异质性为两个单变量比例 Meta 分析各自的 tau²、Q、P 和 I²，不能解释为 SROC
双变量模型的异质性。
}
\examples{
\dontrun{
dat <- data.frame(study = c("A", "B"),
  TP = c(52, 15), FP = c(1, 4), FN = c(19, 5), TN = c(16, 29))
plot_sensspec_forest_meta(dat, sens_axis = c(0, 1),
  spec_axis = c(0.4, 1), row_gap = 1.1)
}
}
\keyword{dplot}
