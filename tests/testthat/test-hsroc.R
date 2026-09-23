testthat::test_that("default HSROC agrees with mada and variance-ratio algebra", {
  d <- data.frame(study = LETTERS[1:5], TP = c(40,32,48,25,55),
    FP = c(8,15,10,20,12), FN = c(10,18,7,15,9), TN = c(72,65,80,58,90))
  fit <- fit_bivariate_dta(d)
  testthat::expect_identical(fit$model_type, "mada::reitsma ruttergatsonis SROC")
  reference <- mada::reitsma(d, method = "reml", correction = .5, correction.control = "single")
  testthat::expect_equal(stats::coef(fit$model), stats::coef(reference))
  mu <- as.numeric(stats::coef(reference))
  slope <- sqrt(reference$Psi[1,1] / reference$Psi[2,2])
  fpr <- 1 - fit$plot_data$sroc$sp
  analytic <- function(x) plogis(mu[1] + slope * (qlogis(x) - mu[2]))
  testthat::expect_equal(fit$plot_data$sroc$se, analytic(fpr), tolerance = 1e-10)
  testthat::expect_true(all(diff(fit$plot_data$sroc$se) >= 0))
  testthat::expect_equal(analytic(plogis(mu[2])), fit$metrics$sensitivity[["est"]])
  x <- seq(0, 1, length.out = 100001)
  y <- analytic(x)
  trap <- sum(diff(x) * (head(y,-1) + tail(y,-1)) / 2)
  testthat::expect_equal(fit$metrics$auc[["est"]], trap, tolerance = 1e-5)
  testthat::expect_true(all(is.na(fit$metrics$auc[c("lwr", "upr")])))
  original <- fit
  observed <- plot_sroc(fit)
  full <- plot_sroc(fit, full_curve = TRUE)
  curve_data <- function(p) Filter(function(layer) inherits(layer$geom, "GeomLine"), p$layers)[[1]]$data
  testthat::expect_equal(curve_data(observed), fit$plot_data$sroc)
  expanded <- curve_data(full)
  labels <- unlist(lapply(ggplot2::ggplot_build(full)$data, function(layer) layer$label))
  testthat::expect_true(any(grepl("AUC:", labels, fixed = TRUE)))
  testthat::expect_false(any(grepl("full curve", labels, fixed = TRUE)))
  testthat::expect_equal(range(expanded$sp), c(0, 1))
  testthat::expect_equal(expanded$se, analytic(1 - expanded$sp), tolerance = 1e-10)
  testthat::expect_true(all(is.finite(expanded$se)))
  testthat::expect_identical(fit, original)
  testthat::expect_error(plot_sroc(fit, full_curve = NA), "full_curve")
  unsupported <- fit
  unsupported$model_type <- "mada::reitsma naive SROC"
  testthat::expect_error(plot_sroc(unsupported, full_curve = TRUE), "ruttergatsonis")
})
