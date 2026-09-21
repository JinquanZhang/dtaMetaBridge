testthat::test_that("explicit QMD adapter retains results and plot controls", {
  d <- data.frame(study = LETTERS[1:5], TP = c(40,32,48,25,55),
    FP = c(8,15,10,20,12), FN = c(10,18,7,15,9), TN = c(72,65,80,58,90))
  original <- fit_metandi(d, seed = 1, n_mc = 100, n_grid = 100)
  fit <- fit_bivariate_dta(d, seed = 1, n_mc = 100, n_grid = 100, sroc_type = "qmd")
  testthat::expect_identical(fit$backend, "qmd")
  testthat::expect_equal(fit$plot_data, original$plot_data)
  testthat::expect_equal(fit$metrics$auc, original$metrics$auc)
  testthat::expect_equal(fit$metrics$sensitivity, original$metrics$se)
  testthat::expect_s3_class(plot_sroc(fit, show_legend = FALSE), "ggplot")
  testthat::expect_true(all(is.finite(posttest_probability(fit)$ppv)))
})
