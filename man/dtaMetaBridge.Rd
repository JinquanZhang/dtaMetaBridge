\name{dtaMetaBridge}
\alias{dtaMetaBridge}
\alias{dta_from_meta}
\alias{fit_forest_summary}
\alias{plot_sensspec_forest}
\alias{plot_sensspec_forest_meta}
\alias{fit_bivariate_dta}
\alias{fit_bivariate_meta}
\alias{plot_sroc}
\alias{posttest_probability}
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
  sens_plot = 1.2, spec_plot = 1.2), width = 10, height = NULL, res = 300)
plot_sensspec_forest_meta(sensitivity_meta, specificity_meta, ...)
fit_bivariate_dta(data, study_col = "study", correction = 0.5,
  correction_control = "single", method = "reml", sroc_type = "ruttergatsonis", n_grid = 1000)
fit_bivariate_meta(sensitivity_meta, specificity_meta, ...)
plot_sroc(fit, show_confidence = TRUE, show_prediction = TRUE,
  output_file = NULL, width = 6, height = 5, dpi = 300)
posttest_probability(fit, prevalence = 0.3, n_sims = 3000, seed = 2026)
}
\arguments{
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
  \item{sens_axis,spec_axis}{Axis tick locations for the two forest plots.}
  \item{column_widths}{Positive named relative widths for the study, count,
  text, and forest-plot columns.}
  \item{correction,correction_control}{Zero-cell continuity-correction controls
  passed to \code{mada::reitsma}.}
  \item{method}{Reitsma estimation method: \code{"reml"}, \code{"ml"}, or
  \code{"fixed"}.}
  \item{sroc_type}{SROC parameterisation. Defaults to \code{"ruttergatsonis"};
  \code{"naive"} is available as a conditional-mean curve.}
  \item{n_grid}{Number of plotting points for the SROC curve.}
  \item{fit}{An object returned by \code{fit_bivariate_dta} or
  \code{fit_bivariate_meta}.}
  \item{width,height,dpi}{Dimensions and resolution used when saving an SROC
  plot with \code{output_file}.}
  \item{show_confidence,show_prediction}{Whether to display the respective
  95 percent regions.}
  \item{prevalence}{One or more pre-test probabilities strictly between zero
  and one.}
  \item{n_sims,seed}{Monte Carlo controls for post-test probability intervals.}
  \item{...}{Further arguments accepted by the underlying plotting function.}
}
\value{
\code{dta_from_meta} returns a 2x2 data frame. The fit functions return a list
containing a \code{mada::reitsma} model, metrics, and plotting data.
\code{posttest_probability} returns a data frame; plotting functions return a
ggplot object or, for forest plots, invisibly return the output filename.
}
\details{
By default the SROC uses the Rutter-Gatsonis parameterisation supplied by
\code{mada::sroc}. The conditional-mean \code{"naive"} curve is available via
\code{sroc_type = "naive"}, but can run in a non-ROC direction when the
study-level covariance is negative. The displayed range is limited to the
observed false-positive-rate range. The full SROC AUC remains a model summary
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
