\name{dtaMetaBridge}
\alias{dtaMetaBridge}
\alias{dta_from_meta}
\alias{fit_forest_summary}
\alias{plot_sensspec_forest}
\alias{fit_bivariate_dta}
\alias{fit_bivariate_meta}
\title{Diagnostic Test Accuracy Meta-Analysis Bridging meta and SROC Models}
\description{
Tools to reconstruct a diagnostic 2x2 table from paired \code{meta::metaprop}
objects, produce paired forest plots, fit a Reitsma bivariate random-effects
model, draw an SROC curve, and calculate post-test probabilities.
}
\usage{
dta_from_meta(sensitivity_meta, specificity_meta)
fit_forest_summary(data, study_col = "study")
plot_sensspec_forest(data, output_file = NULL, study_col = "study",
  summary_override = NULL, heterogeneity = NULL, sens_axis = seq(0, 1, 0.2),
  spec_axis = seq(0, 1, 0.2), column_widths = c(study = 2.8, tp = 0.5,
  fp = 0.5, fn = 0.5, tn = 0.5, sens_text = 1.6, spec_text = 1.6,
  sens_plot = 1.2, spec_plot = 1.2), width = 10, height = NULL, res = 300,
  heterogeneity_cex = 0.66, heterogeneity_x = 0.015, heterogeneity_y = NULL,
  font_family = "serif", study_weights = NULL,
  header_cex = 0.86, study_cex = 0.72, ci_text_cex = 0.68,
  summary_cex = 0.78, summary_count_cex = 0.75, summary_ci_cex = 0.7,
  axis_cex = 0.68, row_gap = NULL, square_col = "grey70", square_max_mm = 5,
  ci_col = "black", ci_lwd = 1.2, diamond_col = "#C00000")
fit_bivariate_dta(data, study_col = "study", correction = 0.5,
  correction_control = c("single", "all", "none"),
  method = c("reml", "ml", "fixed"),
  sroc_type = c("ruttergatsonis", "qmd", "midas", "naive"), n_grid = 1000,
  seed = 2026, n_mc = 3000)
fit_bivariate_meta(sensitivity_meta, specificity_meta, ...)
}
\arguments{

\item{res}{Forest image resolution in dots per inch.}

  \item{header_cex,study_cex,ci_text_cex}{Positive relative sizes for column
  headers, study names/counts, and study confidence-interval text respectively.}
  \item{summary_cex,summary_count_cex,summary_ci_cex}{Positive relative sizes
  for the summary label, summary counts, and summary confidence-interval text.}
  \item{axis_cex}{Positive relative size for axis tick labels. All cex values
  are multipliers, not point sizes.}
  \item{row_gap}{Study-row spacing as a fraction of plot height. NULL uses
  min(0.055, 0.52 / max(1, number of studies - 1)). Must be positive and keep
  the total span of study rows at most 0.52. Summary and axes follow the last
  study; manually specified heterogeneity coordinates stay fixed. Larger
  fonts or squares can still overlap: inspect the plot or increase height.}
  \item{square_col,ci_col,diamond_col}{Valid R color names or hexadecimal colors
  for study squares, CI lines/point markers, and summary diamonds respectively.}
  \item{square_max_mm}{Positive maximum square side length in millimetres.
  Areas remain proportional to available weights within each panel.}
  \item{ci_lwd}{Positive line width for study confidence intervals and point markers.}
  \item{data}{A data frame with a study column and \code{TP}, \code{FP},
  \code{FN}, and \code{TN} count columns.}
  \item{sensitivity_meta,specificity_meta}{Paired random-effects
  \code{meta::metaprop} objects, where event/n encode TP/(TP + FN) and
  TN/(TN + FP), respectively.}
  \item{study_col}{Name of the study-label column.}
  \item{output_file}{Optional destination for a forest or SROC plot. Forest
  plots accept PNG, TIFF, or TIF; SROC output formats are supported by
  \code{ggplot2::ggsave}.}
  \item{summary_override}{Optional summary sensitivity and specificity values
  used for the forest-plot diamonds.}
  \item{heterogeneity}{Optional list with sensitivity and specificity elements,
  each containing \code{i2}, \code{tau2}, \code{q}, and \code{p}.}
  \item{study_weights}{For the data-frame interface, a list of non-negative
  numeric vectors named sensitivity and specificity, aligned to the data rows.
  Square areas are proportional to weights within each panel. The meta-object
  adapter supplies each object's \code{w.random}, aligned by study name;
  do not additionally pass this argument to that adapter.
  Missing or entirely NA weights (e.g. GLMM) produce a warning and equal-size
  squares, not an equal-weight model. Partial missing or invalid weights are rejected.
  Confidence intervals are drawn over squares, with a short vertical point marker.}
  \item{heterogeneity_cex}{Positive text-size multiplier for both heterogeneity lines;
  default 0.66. This is a relative size, not a size in points.}
  \item{font_family}{Font family for all SROC text and all forest-plot text, including headers,
  study labels, numbers, axis labels, and heterogeneity. Default \code{"serif"};
  \code{"sans"} and \code{"mono"} are also supported generic family names.
  Other fonts must be supported by the graphics device. Sizes and bold
  emphasis are unchanged. Also accepted by \code{plot_sensspec_forest_meta(...)}.}
  \item{heterogeneity_x}{Left-aligned text anchor for both heterogeneity lines,
  in normalized plot coordinates (0 = left, 1 = right). Default 0.015.}
  \item{heterogeneity_y}{Two normalized vertical coordinates, in sensitivity,
  specificity order (0 = bottom, 1 = top). Default NULL positions the lines
  automatically below the summary row. Custom positions do not reserve space
  or prevent overlap with other plot elements. These three options can also
  be passed through \code{plot_sensspec_forest_meta(...)}.}
  \item{sens_axis,spec_axis}{Axis tick locations for the two forest plots.}
  \item{column_widths}{Positive named relative widths for the study, count,
  text, and forest-plot columns.}
  \item{correction,correction_control}{Zero-cell continuity-correction controls
  passed to \code{mada::reitsma}. Omit for QMD and Midas backends.}
  \item{method}{Reitsma estimation method: \code{"reml"}, \code{"ml"}, or
  \code{"fixed"}. Omit for QMD and Midas backends.}
  \item{sroc_type}{Defaults to \code{"ruttergatsonis"}, the Reitsma model
  with Rutter-Gatsonis HSROC parameterisation. \code{"qmd"} selects the original QMD functions.
  \code{"midas"} selects a binomial GLMM with the
  Midas curve formula. \code{"ruttergatsonis"} uses the mada backend;
  \code{"naive"} is available as a conditional-mean curve.}
  \item{n_grid}{Number of plotting points for the SROC curve.}
  \item{n_mc}{Fixed-effect Monte Carlo simulations for the QMD AUC.}

  \item{width,height}{Saved plot width and height in inches. Forest height is automatic when NULL.}

  \item{seed}{Random seed for QMD Monte Carlo simulation.}
  \item{...}{Further arguments passed to fit_bivariate_dta.}
}
\value{
\code{dta_from_meta} returns a 2x2 data frame. The fit functions return a list
containing a \code{glmerMod} or \code{reitsma} model, metrics, and plotting data.
\code{posttest_probability} returns a data frame; plotting functions return a
ggplot object or, for forest plots, invisibly return the output filename.
}
\details{
The default uses mada Reitsma REML estimation and Rutter-Gatsonis HSROC
parameterisation for a model without covariates (Reitsma et al., 2005;
Rutter and Gatsonis, 2001; Harbord et al., 2007). Estimation uses the
logit-normal approximation, not the original Bayesian likelihood.
AUC integrates the fitted curve from FPR 0 to 1; pAUC is the unstandardised
integral over observed FPR. No AUC confidence interval is estimated.
The optional QMD curve slope is sign(rho) times the random-effect SD ratio; positive
sensitivity-specificity correlation can therefore yield a reverse curve.
Its simulated AUC interval conditions on the estimated random-effect covariance.
The optional Midas backend fits raw binomial counts with lme4 Laplace estimation and applies
the curve formula in Midas 2.00 (2008-12-21), with variances floored at 0.001.
Full AUC uses 500-point trapezoidal integration. Contours use the Midas
F(2, k-2) radius. AUC confidence limits are not estimated. Estimation can
differ from Stata adaptive quadrature; this is not a verified numerical clone.
The conditional-mean \code{"naive"} curve is available via
\code{sroc_type = "naive"}, but can run in a non-ROC direction when the
study-level covariance is negative. Midas shows the full curve (including
extrapolation); mada alternatives show the observed range. The SROC AUC is a model summary
and should not be equated with a continuous-score ROC AUC from an individual
study.
}
\examples{
\dontrun{
sens <- meta::metaprop(TP, TP + FN, studlab = study, data = dta,
                       method = "GLMM", method.tau = "ML")
spec <- meta::metaprop(TN, TN + FP, studlab = study, data = dta,
                       method = "GLMM", method.tau = "ML")
fit <- fit_bivariate_meta(sens, spec)
plot_sroc(fit)
posttest_probability(fit, c(0.1, 0.3))
}
}
