testthat::test_that("2x2 inputs and meta-object adapter are validated", {
  toy_dta <- data.frame(study = c("A", "B", "C"), TP = c(8, 12, 16), FP = c(3, 2, 1), FN = c(2, 4, 3), TN = c(17, 18, 20))
  forest <- dtaMetaBridge:::.forest_data(toy_dta)
  testthat::expect_equal(nrow(forest$studies), 3)
  testthat::expect_equal(forest$summary$TP, 36)
  testthat::expect_true(all(forest$studies$sens > 0 & forest$studies$sens < 1))
  testthat::expect_error(dtaMetaBridge:::.assert_dta_data(transform(toy_dta, TP = -1)))

if (requireNamespace("meta", quietly = TRUE)) {
  sensitivity_meta <- meta::metaprop(TP, TP + FN, studlab = study, data = toy_dta, method = "GLMM", method.tau = "ML")
  specificity_meta <- meta::metaprop(TN, TN + FP, studlab = study, data = toy_dta, method = "GLMM", method.tau = "ML")
  reconstructed <- dta_from_meta(sensitivity_meta, specificity_meta)
  testthat::expect_identical(reconstructed[, c("TP", "FP", "FN", "TN")], toy_dta[, c("TP", "FP", "FN", "TN")])
  forest_summary <- fit_forest_summary(toy_dta)
  testthat::expect_true(is.finite(forest_summary$heterogeneity$sensitivity$i2))
  testthat::expect_true(is.finite(forest_summary$heterogeneity$sensitivity$q))
}
testthat::expect_error(plot_sensspec_forest(toy_dta, column_widths = c(study = 1)), "column_widths")
})

testthat::test_that("naive is the default SROC parameterisation", {
  toy_dta <- data.frame(
    study = LETTERS[1:5],
    TP = c(40, 32, 48, 25, 55), FP = c(8, 15, 10, 20, 12),
    FN = c(10, 18, 7, 15, 9), TN = c(72, 65, 80, 58, 90)
  )
  fit <- fit_bivariate_dta(toy_dta, n_grid = 100)
  testthat::expect_identical(fit$model_type, "mada::reitsma naive SROC")
  testthat::expect_true(is.finite(fit$metrics$auc[["est"]]))
  output_file <- tempfile(fileext = ".png")
  testthat::expect_s3_class(plot_sroc(fit, output_file = output_file), "ggplot")
  testthat::expect_true(file.exists(output_file))
})
