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

testthat::test_that("Midas curve is monotone, anchored and consistently integrated", {
  toy_dta <- data.frame(
    study = LETTERS[1:5],
    TP = c(40, 32, 48, 25, 55), FP = c(8, 15, 10, 20, 12),
    FN = c(10, 18, 7, 15, 9), TN = c(72, 65, 80, 58, 90)
  )
  fit <- fit_bivariate_dta(toy_dta, n_grid = 100, sroc_type = "midas")
  testthat::expect_identical(fit$backend, "midas")
  testthat::expect_true(all(diff(fit$plot_data$sroc$se) <= 0))
  testthat::expect_equal(fit$plot_data$sroc$se[c(1, 100)], c(1, 0))
  a <- fit$random_effects$alpha
  b <- fit$random_effects$beta
  f <- function(x) plogis((a - qlogis(x) / b) / b)
  testthat::expect_equal(f(fit$metrics$specificity[["est"]]), fit$metrics$sensitivity[["est"]])
  testthat::expect_equal(fit$metrics$auc[["est"]], integrate(f, 0, 1)$value, tolerance = 1e-3)
  post <- posttest_probability(fit, prevalence = .3)
  se <- fit$metrics$sensitivity[["est"]]
  sp <- fit$metrics$specificity[["est"]]
  testthat::expect_equal(post$ppv, unname(.3 * se / (.3 * se + .7 * (1 - sp))))
  testthat::expect_equal(post$npv, unname(.7 * sp / (.7 * sp + .3 * (1 - se))))
  old <- fit_bivariate_dta(toy_dta, sroc_type = "ruttergatsonis")
  testthat::expect_s3_class(old$model, "reitsma")
  testthat::expect_true(is.finite(fit$metrics$auc[["est"]]))
  output_file <- tempfile(fileext = ".png")
  testthat::expect_s3_class(plot_sroc(fit, output_file = output_file), "ggplot")
  testthat::expect_true(file.exists(output_file))
})
