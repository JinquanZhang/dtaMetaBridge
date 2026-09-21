# Run from the package root: Rscript tools/render-examples.R
devtools::load_all(".", quiet = TRUE)
dta <- data.frame(study = paste("Study", LETTERS[1:5]),
  TP = c(40,32,48,25,55), FP = c(8,15,10,20,12),
  FN = c(10,18,7,15,9), TN = c(72,65,80,58,90))
sens <- meta::metaprop(TP, TP + FN, studlab = study, data = dta, method = "GLMM", method.tau = "ML")
spec <- meta::metaprop(TN, TN + FP, studlab = study, data = dta, method = "GLMM", method.tau = "ML")
fit <- fit_bivariate_meta(sens, spec)
plot_sensspec_forest_meta(sens, spec, output_file = "inst/figures/forest-example.png")
plot_sroc(fit, output_file = "inst/figures/sroc-example.png")
print(fit$metrics)
