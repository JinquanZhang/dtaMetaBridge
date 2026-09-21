\name{fit_metandi}
\alias{fit_metandi}
\alias{plot_metandi}
\title{Original SROC Functions from HFphf_meta analysis.qmd}
\description{The original user-provided QMD fitting and plotting functions, with namespace qualification and base pipes. Statistical formulas are retained.}
\usage{
fit_metandi(data, seed = 2026, n_mc = 3000, n_grid = 1000)
plot_metandi(model_obj, digits = 2, show_conf = TRUE, show_pred = TRUE,
  show_legend = TRUE, custom_se = NULL, custom_sp = NULL, custom_auc = NULL)
}
\arguments{
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
