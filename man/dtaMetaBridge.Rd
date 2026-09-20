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
\title{Diagnostic Test Accuracy Meta-Analysis Bridging meta and HSROC Models}
\description{
Tools to reconstruct a diagnostic 2x2 table from paired \code{meta::metaprop}
objects, produce paired forest plots, fit a Reitsma bivariate random-effects
model, draw a Rutter-Gatsonis HSROC curve, and calculate post-test probabilities.
}
\usage{
dta_from_meta(sensitivity_meta, specificity_meta)
fit_forest_summary(data, study_col = "study")
plot_sensspec_forest(data, output_file = NULL, study_col = "study", ...)
plot_sensspec_forest_meta(sensitivity_meta, specificity_meta, ...)
fit_bivariate_dta(data, study_col = "study", correction = 0.5,
  correction_control = "single", method = "reml", n_grid = 1000)
fit_bivariate_meta(sensitivity_meta, specificity_meta, ...)
plot_sroc(fit, show_confidence = TRUE, show_prediction = TRUE)
posttest_probability(fit, prevalence = 0.3, n_sims = 3000, seed = 2026)
}
\arguments{
  \item{data}{A data frame with a study column and \code{TP}, \code{FP},
  \code{FN}, and \code{TN} count columns.}
  \item{sensitivity_meta,specificity_meta}{Paired random-effects
  \code{meta::metaprop} objects, where event/n encode TP/(TP + FN) and
  TN/(TN + FP), respectively.}
  \item{study_col}{Name of the study-label column.}
  \item{output_file}{Optional PNG, TIFF, or TIF destination for the forest plot.}
  \item{correction,correction_control}{Zero-cell continuity-correction controls
  passed to \code{mada::reitsma}.}
  \item{method}{Reitsma estimation method: \code{"reml"}, \code{"ml"}, or
  \code{"fixed"}.}
  \item{n_grid}{Number of plotting points for the HSROC curve.}
  \item{fit}{An object returned by \code{fit_bivariate_dta} or
  \code{fit_bivariate_meta}.}
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
The HSROC curve uses the Rutter-Gatsonis parameterisation supplied by
\code{mada::sroc}, not a conditional-mean SROC. Its displayed range is limited
to the observed false-positive-rate range. The full HSROC AUC remains a model
summary and should not be equated with a continuous-score ROC AUC from an
individual study.
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
