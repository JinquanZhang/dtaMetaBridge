# Run from the package root: Rscript tools/render-examples.R
devtools::load_all(".", quiet = TRUE)
dta <- data.frame(study = c("Parcha 2021", "Tada 2021", "Forsyth 2021", "Reddy 2022", "Ariyaratnam 2024", "Akerman 2025", "Wang 2023", "Nguyen 2025", "Rahi 2026"),
  TP = c(169,113,11,163,61,117,61,110,71), FP = c(29,5,16,3,11,42,2,36,27),
  FN = c(122,113,39,222,27,123,140,8,9), TN = c(350,173,11,99,21,214,116,126,85))
sens <- meta::metaprop(TP, TP + FN, studlab = study, data = dta, method = "GLMM", method.tau = "ML")
spec <- meta::metaprop(TN, TN + FP, studlab = study, data = dta, method = "GLMM", method.tau = "ML")
fit <- fit_bivariate_meta(sens, spec)
plot_sensspec_forest_meta(sens, spec, output_file = "inst/figures/forest-example.png")
plot_sroc(fit, output_file = "inst/figures/sroc-example.png", width = 6, height = 6)
print(fit$metrics)
