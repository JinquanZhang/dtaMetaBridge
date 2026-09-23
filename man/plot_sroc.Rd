\name{plot_sroc}
\encoding{UTF-8}
\alias{plot_sroc}
\title{绘制 SROC 曲线及置信区域、预测区域}
\description{从双变量模型结果绘图，支持完整 Rutter--Gatsonis 曲线及外观调整。}
\usage{
plot_sroc(fit, show_confidence = TRUE, show_prediction = TRUE,
  output_file = NULL, width = 6, height = 5, dpi = 300, ...,
  full_curve = FALSE, show_legend = TRUE, digits = 2,
  custom_se = NULL, custom_sp = NULL, custom_auc = NULL,
  legend_position = "bottomright", legend_text_size = 3.5,
  legend_bg = "#F8F9FA", font_family = "sans", base_size = 14,
  title = "SROC with Prediction & Confidence Contours",
  show_study_labels = TRUE, study_size = 4, study_label_size = 2.5,
  summary_size = 4.5, summary_col = "#C0392B",
  sroc_col = "#2C3E50", sroc_linewidth = 1.2,
  confidence_col = "#2980B9", confidence_alpha = 0.2,
  prediction_col = "#BDC3C7", prediction_alpha = 0.15,
  auc_digits = NULL, x_breaks = seq(0, 1, 0.2),
  y_breaks = seq(0, 1, 0.2))
}
\arguments{
\item{fit}{\code{fit_bivariate_meta} 或 \code{fit_bivariate_dta} 返回的拟合对象。}
\item{show_confidence}{是否显示汇总点的 95\% 置信区域。}
\item{show_prediction}{是否显示 95\% 预测区域。}
\item{output_file}{保存路径；NULL 不保存。格式由 \code{ggplot2::ggsave} 支持情况决定。}
\item{width,height}{保存宽度、高度，英寸。}
\item{dpi}{保存分辨率。}
\item{...}{保留位置；未识别的额外参数会报错。其后的参数须使用名称传入。}
\item{full_curve}{TRUE 绘制 FPR 0--1 的完整曲线；仅支持 ruttergatsonis。FALSE 使用已存曲线，默认限制在观察范围。}
\item{show_legend}{是否显示图内注释式图例。}
\item{digits}{灵敏度、特异度的小数位数，0--10 的整数。}
\item{custom_se,custom_sp,custom_auc}{自定义图例显示文字；NULL 使用计算值，不修改模型。}
\item{legend_position}{四角位置 "bottomright"、"bottomleft"、"topright"、"topleft"，或 c(x,y) 指定图例框左下角。屏幕坐标 0 为左/下、1 为右/上，不受特异度反向坐标影响。}
\item{legend_text_size}{图例字号，毫米，正数；标题稍大。}
\item{legend_bg}{图例背景颜色，可用 "transparent"。}
\item{font_family}{所有文字字体，默认 "sans"。}
\item{base_size}{主题基础字号，磅；注释文字另行控制。}
\item{title}{标题；NULL 不显示。}
\item{show_study_labels}{是否在研究圆点内显示编号。}
\item{study_size}{研究点大小，毫米。}
\item{study_label_size}{研究编号字号，毫米。}
\item{summary_size}{汇总点大小，毫米。}
\item{summary_col}{汇总点颜色。}
\item{sroc_col}{SROC 线及研究编号颜色。}
\item{sroc_linewidth}{SROC 线宽，毫米。}
\item{confidence_col}{置信区域颜色。}
\item{confidence_alpha}{置信区域填充透明度，0--1。}
\item{prediction_col}{预测区域颜色。}
\item{prediction_alpha}{预测区域填充透明度，0--1。}
\item{auc_digits}{AUC 小数位数，0--10 整数；NULL 时 QMD 使用 digits，其他方法使用 3。}
\item{x_breaks,y_breaks}{横纵轴刻度，0--1 内不重复数值；不改变坐标范围和算法。}
}
\details{
横轴为特异度，由 1 递减到 0；纵轴为灵敏度。
默认拟合方法为 Reitsma 双变量 REML 加 Rutter--Gatsonis 曲线。
full_curve 只改变绘图，不改变汇总估计、AUC 或区域；观察范围外属于模型外推。
标准方法的 AUC 为完整曲线积分，未估计 AUC 置信区间。
图例宽度固定，高度随字号变化；框超出面板时报错，过长文字仍可能溢出。
置信区域描述汇总点不确定性，预测区域描述新研究可能的真实效应，不是曲线置信带。
}
\value{可继续用 ggplot2 修改或保存的 ggplot 对象。}
\seealso{\code{\link{fit_bivariate_meta}}, \code{\link{posttest_probability}}}
\examples{
d <- data.frame(study = LETTERS[1:6],
  TP = c(35, 42, 28, 60, 45, 70), FN = c(15, 8, 22, 20, 15, 10),
  TN = c(80, 65, 90, 55, 75, 60), FP = c(20, 35, 10, 45, 25, 40))
fit <- fit_bivariate_dta(d)
plot_sroc(fit, full_curve = TRUE, legend_position = "bottomright",
  show_study_labels = FALSE, sroc_col = "navy")
}
