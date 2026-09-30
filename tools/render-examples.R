# Run from the package root: Rscript tools/render-examples.R
devtools::load_all(".", quiet = TRUE)
dta <- utils::read.csv("example-dta.csv", check.names = FALSE)
sens <- meta::metaprop(TP, TP + FN, studlab = study, data = dta, method = "GLMM", method.tau = "ML")
spec <- meta::metaprop(TN, TN + FP, studlab = study, data = dta, method = "GLMM", method.tau = "ML")
fit <- fit_sroc(dta, backend = "frequency", auc_boot = 2000)
plot_sensspec_forest_meta(sens, spec, output_file = "inst/figures/forest-current.png")
plot_sroc(fit, output_file = "inst/figures/sroc-current.png", width = 6, height = 6,
  full_curve = TRUE)
print(fit$metrics)
