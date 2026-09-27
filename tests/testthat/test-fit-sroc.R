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

test_that("user data retain duplicate rows, zero cells and original counts", {
  d <- data.frame(study=c("Briasoulis","Kravchenko","Ridouani","Ridouani"),
    Year=c("2023","2024","2018","2018"), TP=c(52,15,20,10),
    FP=c(1,4,1,0), FN=c(19,5,4,14), TN=c(16,29,19,20))
  testthat::skip_if_not_installed("dtametaTMB")
  fits <- lapply(1:5, function(i) suppressWarnings(fit_sroc(d, backend="frequency", sroc_type=i, auc_boot=0)))
  for (i in 1:5) {
    f <- fits[[i]]
    expect_identical(f$input_data, d)
    expect_equal(nrow(f$plot_data$studies), 4L)
    expect_equal(anyDuplicated(f$plot_data$studies$Study), 0L)
    p <- f$curve_parameters
    expect_equal(.sroc_values(plogis(p$mu[2]), p$mu, p$slope), plogis(p$mu[1]))
    grid <- seq(0,1,length.out=100001)
    y <- .sroc_values(grid, p$mu, p$slope)
    if (p$slope > 0) {
      expect_true(all(diff(y) >= 0))
      expect_equal(sum(diff(grid)*(head(y,-1)+tail(y,-1))/2), unname(f$metrics$auc[1]), tolerance=1e-5)
    } else expect_true(is.na(f$metrics$auc[1]))
    expect_s3_class(plot_sroc(f, full_curve=TRUE, show_legend=FALSE), "ggplot")
  }
  expect_error(fit_sroc(d, backend="frequency", sroc_type=6), "1 to 5")
  expect_error(fit_sroc(d, n_grid=NA_real_), "n_grid")
  expect_error(fit_sroc(d, conf_level=1), "conf_level")
})
