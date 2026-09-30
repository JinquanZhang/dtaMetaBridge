\name{plot_sroc}
\alias{plot_sroc}
\title{绘制 SROC 曲线、置信轮廓和预测轮廓}
\description{
绘制 \code{\link{fit_sroc}} 结果的 SROC 图，包括研究点、汇总点、
SROC 曲线、95\% 置信轮廓和 95\% 预测轮廓。
}
\usage{
plot_sroc(fit, show_confidence = TRUE, show_prediction = TRUE,
  output_file = NULL, width = 6, height = 5, dpi = 300,
  full_curve = FALSE, show_legend = TRUE, digits = 2,
  custom_se = NULL, custom_sp = NULL, custom_auc = NULL,
  legend_position = "bottomright", legend_text_size = 3.5,
  legend_bg = "#F8F9FA", font_family = "sans", base_size = 14,
  title = "SROC with Prediction & Confidence Contours",
  show_study_labels = TRUE, study_size = 4, study_label_size = 2.5,
  study_col = "#7F8C8D", study_label_col = "black",
  summary_size = 4.5, summary_col = "#C0392B",
  sroc_col = "#2C3E50", sroc_linewidth = 1.2,
  confidence_col = "#2980B9", confidence_alpha = 0.2,
  prediction_col = "#BDC3C7", prediction_alpha = 0.15,
  auc_digits = NULL, x_breaks = seq(0, 1, 0.2),
  y_breaks = seq(0, 1, 0.2), x_axis = c("specificity", "fpr"))
}
\arguments{
  \item{fit}{由 \code{fit_sroc()} 返回的对象。}
  \item{show_confidence,show_prediction}{是否显示相应的 95\% 轮廓。}
  \item{output_file}{可选输出文件路径；为空时返回 \code{ggplot} 对象。}
  \item{width,height,dpi}{导出图像的英寸宽、高与分辨率。}
  \item{full_curve}{是否外推并显示完整 FPR 0--1 的 SROC。默认仅显示
  数据支持范围，避免把外推部分误解为观测证据。}
  \item{show_legend,digits}{是否显示图例及汇总指标的小数位数。}
  \item{custom_se,custom_sp,custom_auc}{覆盖图例中相应指标文本的可选字符串。}
  \item{legend_position,legend_text_size,legend_bg}{图例位置、字号和背景色。}
  \item{font_family,base_size,title}{全图字体、基础字号和标题。}
  \item{show_study_labels,study_size,study_label_size,study_col,study_label_col}{研究点及标签的显示控制。研究点绘制在曲线和轮廓下方；标签默认黑色。}
  \item{summary_size,summary_col}{汇总点（45 度旋转正方形）的大小和颜色。}
  \item{sroc_col,sroc_linewidth}{SROC 曲线颜色和线宽。}
  \item{confidence_col,confidence_alpha,prediction_col,prediction_alpha}{两类轮廓的颜色和透明度。}
  \item{auc_digits}{AUC 单独的小数位数；为空时使用 \code{digits}。}
  \item{x_breaks,y_breaks}{横轴和纵轴刻度。}
  \item{x_axis}{\code{"specificity"} 为从 1 到 0 的特异度横轴；
  \code{"fpr"} 为从 0 到 1 的 1 - 特异度横轴。}
}
\value{一个 \code{ggplot} 对象；若指定 \code{output_file}，同时保存图像。}
\details{
频率学结果的 AUC 区间来自 \code{fit_sroc()} 的 Bootstrap；贝叶斯结果使用后验
可信区间。预测轮廓描述未来研究的潜在真实准确性变异，不包括新研究有限样本的
二项抽样误差。图例中的轮廓水平继承 \code{fit_sroc()} 的 \code{conf_level}。
图例框会按实际文字宽度和字号测量后自适应布局；若输出画布过小，
函数会报错而不截断文字。研究点为未调整的单研究灵敏度和特异度，轮廓反映模型
中的汇总点不确定性及研究间变异。
}
\examples{
\dontrun{
fit <- fit_sroc(dat, backend = "frequency", auc_boot = 200)
plot_sroc(fit, x_axis = "fpr", full_curve = TRUE,
  output_file = "sroc.png")
}
}
\keyword{dplot}
