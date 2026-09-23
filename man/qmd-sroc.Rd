\name{fit_metandi}
\alias{fit_metandi}
\alias{plot_metandi}
\title{Original SROC Functions from HFphf_meta analysis.qmd}
\description{The original user-provided QMD fitting and plotting functions, with namespace qualification and base pipes. Statistical formulas are retained.}
\usage{
fit_metandi(data, seed = 2026, n_mc = 3000, n_grid = 1000)
plot_metandi(model_obj, digits = 2, show_conf = TRUE, show_pred = TRUE,
  show_legend = TRUE, custom_se = NULL, custom_sp = NULL, custom_auc = NULL,
  legend_position = "bottomright", legend_text_size = 3.5,
  legend_bg = "#F8F9FA", font_family = "sans", base_size = 14,
  title = "SROC with Prediction & Confidence Contours",
  show_study_labels = TRUE, study_size = 4, study_label_size = 2.5,
  summary_size = 4.5, summary_col = "#C0392B",
  sroc_col = "#2C3E50", sroc_linewidth = 1.2,
  confidence_col = "#2980B9", confidence_alpha = 0.2,
  prediction_col = "#BDC3C7", prediction_alpha = 0.15,
  auc_digits = digits, x_breaks = seq(0, 1, 0.2), y_breaks = seq(0, 1, 0.2))
}
\arguments{
\item{legend_position}{One of bottomright, bottomleft, topright, topleft,
or numeric c(x, y) specifying the legend box's lower-left corner in screen
coordinates (0 = left/bottom, 1 = right/top), independent of the reversed
specificity axis. The box has fixed width 0.44 and a height scaled to the
text size. Out-of-panel boxes are rejected; long custom text may still
overflow. The legend is drawn as annotations, not a standard ggplot legend.}
\item{legend_text_size}{Positive legend text size in millimetres (heading is slightly larger).}
\item{legend_bg}{Legend background color; use "transparent" for no fill.}
\item{font_family}{Font family for axes, title, study labels and legend;
must be supported by the output graphics device.}
\item{base_size}{Positive base theme font size in points; annotation sizes are controlled separately.}
\item{title}{Plot title, or NULL to remove it.}
\item{show_study_labels}{Logical; show study-number labels inside study circles.}
\item{study_size,study_label_size,summary_size}{Positive sizes in millimetres
for study markers, study-number labels, and summary markers.}
\item{summary_col,sroc_col,confidence_col,prediction_col}{R color names or hexadecimal colors.}
\item{sroc_linewidth}{Positive SROC line width in millimetres, also applied to its legend key.}
\item{confidence_alpha,prediction_alpha}{Region fill opacity between 0 and 1.}
\item{auc_digits}{AUC decimal places, integer 0 to 10. In plot_sroc, NULL
uses digits for QMD fits and 3 otherwise. Missing AUC intervals are omitted,
not displayed as NA.}
\item{x_breaks,y_breaks}{Distinct numeric axis ticks between 0 and 1.
These do not change the axis range or the curve calculation.}
\item{data}{Four integer columns TP, FP, FN, TN and Study or study labels.}
\item{seed}{Random seed; changes global random state as in the original.}
\item{n_mc}{Number of fixed-effect simulations for AUC.}
\item{n_grid}{Number of grid points.}
\item{model_obj}{Output of fit_metandi.}
\item{digits}{Displayed decimal places.}
\item{show_conf,show_pred,show_legend}{Show contours and legend.}
\item{custom_se,custom_sp,custom_auc}{Optional legend text overrides; do not change estimates.}
}
\value{fit_metandi returns metrics, plot_data and a glmerMod model. plot_metandi returns a ggplot and prints summary values.}
\details{The curve uses sign(rho) times the random-effect SD ratio. It can reverse with positive sensitivity-specificity correlation. The AUC interval simulates only fixed effects, conditional on estimated random-effect covariance.}
