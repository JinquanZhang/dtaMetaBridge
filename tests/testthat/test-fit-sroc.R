test_that("five formulas agree with independent algebra and principal axis", {
  P <- matrix(c(.5, .2, .2, .8), 2)
  eig <- eigen(P)$vectors[, 1]
  expected <- c(.2/.8, eig[1]/eig[2], (.5+.2)/(.8+.2), .5/.2, sqrt(.5/.8))
  for (i in 1:5) expect_equal(.sroc_slope(P, i), expected[i], tolerance = 1e-12)
  P[1,2] <- P[2,1] <- -.2
  eig <- eigen(P)$vectors[,1]
  expect_equal(.sroc_slope(P, 2L), eig[1]/eig[2])
  expect_error(.sroc_slope(diag(c(.5,.8)), 2L), "zero covariance")
  expect_error(.sroc_slope(diag(c(.5,.8)), 4L), "zero covariance")
  expect_error(.sroc_slope(matrix(c(1,-.5,-.5,.5),2), 3L), "denominator")
  expect_equal(.sroc_values(c(0,.5,1), c(0,0), 0), rep(.5,3))
})

test_that("direct Rutter-Gatsonis fit retains rows, zero cells and original counts", {
  d <- data.frame(study=c("Briasoulis","Kravchenko","Ridouani","Ridouani"),
    Year=c("2023","2024","2018","2018"), TP=c(52,15,20,10),
    FP=c(1,4,1,0), FN=c(19,5,4,14), TN=c(16,29,19,20))
  testthat::skip_if_not_installed("dtametaTMB")
  f <- suppressWarnings(fit_sroc(d, backend="frequency", sroc_type=5, auc_boot=0))
  {
    expect_identical(f$input_data, d)
    expect_equal(nrow(f$plot_data$studies), 4L)
    expect_equal(anyDuplicated(f$plot_data$studies$Study), 0L)
    p <- f$curve_parameters
    expect_identical(f$model_type, "dtametaTMB::fitRutterGatsonis")
    expect_equal(.ruttergatsonis_values(c(.1, .5, .9), p$Lambda, p$beta),
      .sroc_values(c(.1, .5, .9), p$mu, p$slope), tolerance = 1e-12)
    expect_equal(.sroc_values(plogis(p$mu[2]), p$mu, p$slope), plogis(p$mu[1]))
    grid <- seq(0,1,length.out=100001)
    y <- .sroc_values(grid, p$mu, p$slope)
    if (p$slope > 0) {
      expect_true(all(diff(y) >= 0))
      expect_equal(sum(diff(grid)*(head(y,-1)+tail(y,-1))/2), unname(f$metrics$auc[1]), tolerance=1e-5)
    } else expect_true(is.na(f$metrics$auc[1]))
    expect_s3_class(plot_sroc(f, full_curve=TRUE, show_legend=FALSE), "ggplot")
  }
  expect_error(fit_sroc(d, backend="frequency", sroc_type=1), "requires sroc_type = 5")
  expect_error(fit_sroc(d, backend="frequency", sroc_type=6), "1 to 5")
  expect_error(fit_sroc(d, n_cores=0), "n_cores")
  expect_error(fit_sroc(d, n_grid=NA_real_), "n_grid")
  expect_error(fit_sroc(d, conf_level=1), "conf_level")
})

test_that("SROC legend box adapts to its text and visible entries", {
  testthat::skip_if_not_installed("dtametaTMB")
  d <- data.frame(study = LETTERS[1:6], TP = c(35,42,28,60,45,70),
    FN = c(15,8,22,20,15,10), TN = c(80,65,90,55,75,60), FP = c(20,35,10,45,25,40))
  f <- fit_sroc(d, auc_boot = 0)
  rect <- function(p) Filter(function(x) "xmin" %in% names(x), ggplot2::ggplot_build(p)$data)[[1]]
  standard <- rect(plot_sroc(f))
  long <- rect(plot_sroc(f, custom_auc = "AUC with a deliberately much longer explanatory label for sizing"))
  compact <- rect(plot_sroc(f, show_confidence = FALSE, show_prediction = FALSE))
  expect_gt(abs(long$xmax - long$xmin), abs(standard$xmax - standard$xmin))
  expect_gt(standard$ymax - standard$ymin, compact$ymax - compact$ymin)
})

test_that("SROC legend measures text against the export size", {
  testthat::skip_if_not_installed("dtametaTMB")
  d <- data.frame(study = LETTERS[1:6], TP = c(35,42,28,60,45,70),
    FN = c(15,8,22,20,15,10), TN = c(80,65,90,55,75,60), FP = c(20,35,10,45,25,40))
  f <- fit_sroc(d, auc_boot = 0)
  out <- tempfile(fileext = ".png")
  on.exit(unlink(out), add = TRUE)
  expect_s3_class(plot_sroc(f, output_file = out, width = 4, height = 4, dpi = 72), "ggplot")
  expect_true(file.exists(out))
  expect_error(plot_sroc(f, output_file = tempfile(fileext = ".png"), width = 3, height = 3, dpi = 72),
    "Legend does not fit")
})

test_that("SROC x axis switches between specificity and FPR", {
  testthat::skip_if_not_installed("dtametaTMB")
  d <- data.frame(study = LETTERS[1:6], TP = c(35,42,28,60,45,70),
    FN = c(15,8,22,20,15,10), TN = c(80,65,90,55,75,60), FP = c(20,35,10,45,25,40))
  f <- fit_sroc(d, auc_boot = 0)
  p_sp <- plot_sroc(f, show_legend = FALSE, x_axis = "specificity")
  p_fpr <- plot_sroc(f, show_legend = FALSE, x_axis = "fpr")
  expect_identical(p_sp$labels$x, "Specificity")
  expect_identical(p_fpr$labels$x, "1 - Specificity")
  b_sp <- ggplot2::ggplot_build(p_sp)
  b_fpr <- ggplot2::ggplot_build(p_fpr)
  expect_lt(max(abs(b_fpr$data[[1]]$x - (1 + b_sp$data[[1]]$x))), 1e-12)
})

test_that("SROC plot structure keeps study points beneath contours and labels the fitted level", {
  testthat::skip_if_not_installed("dtametaTMB")
  d <- data.frame(study = LETTERS[1:6], TP = c(35,42,28,60,45,70),
    FN = c(15,8,22,20,15,10), TN = c(80,65,90,55,75,60), FP = c(20,35,10,45,25,40))
  f <- fit_sroc(d, conf_level = .90, auc_boot = 0)
  b <- ggplot2::ggplot_build(plot_sroc(f))
  expect_identical(b$data[[1]]$shape[1], 21)
  expect_identical(b$data[[2]]$fill[1], "#BDC3C7")
  labels <- unlist(lapply(b$data, function(x) if ("label" %in% names(x)) x$label else character()))
  expect_true(any(grepl("90% Confidence Contour", labels, fixed = TRUE)))
  expect_true(any(grepl("90% Prediction Contour", labels, fixed = TRUE)))
})

test_that("Bayesian type-5 curve and AUC agree with the fitted meta4diag object", {
  testthat::skip_if_not_installed("meta4diag")
  testthat::skip_if_not_installed("INLA")
  d <- data.frame(study = LETTERS[1:6], TP = c(35,42,28,60,45,70),
    FN = c(15,8,22,20,15,10), TN = c(80,65,90,55,75,60), FP = c(20,35,10,45,25,40))
  f <- fit_sroc(d, backend = "bayes", sroc_type = 5, posterior_samples = 100, seed = 2026)
  x <- 1 - f$plot_data$sroc$sp
  expected_curve <- .sroc_values(x, f$curve_parameters$mu, f$curve_parameters$slope)
  expect_equal(f$plot_data$sroc$se, expected_curve, tolerance = 1e-12)
  reference_auc <- meta4diag::AUC(f$model, sroc.type = 5, est.type = "mean")
  expect_equal(unname(f$metrics$auc["est"]), unname(reference_auc["est"]), tolerance = 1e-12)
})
