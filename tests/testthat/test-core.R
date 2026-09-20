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
}
})
