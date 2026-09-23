testthat::test_that("forest areas use aligned meta weights and CI is drawn over squares", {
  d <- data.frame(study = c("A", "B", "C"), TP = c(8, 12, 16),
                  FP = c(3, 2, 1), FN = c(2, 4, 3), TN = c(17, 18, 20))
  se <- meta::metaprop(TP, TP + FN, studlab = study, data = d, method = "Inverse")
  sp <- meta::metaprop(TN, TN + FP, studlab = study, data = d[3:1, ], method = "Inverse")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot_sensspec_forest_meta(se, sp)
  grobs <- as.list(grid::grid.grab()$children)
  rectangles <- Filter(function(g) inherits(g, "rect"), grobs)
  expected <- c(se$w.random / max(se$w.random), rev(sp$w.random) / max(sp$w.random))
  testthat::expect_equal(unname(vapply(rectangles, function(g) as.numeric(g$width)^2 / 25, numeric(1))), expected)

  grid::grid.newpage()
  dtaMetaBridge:::.draw_forest_panel(c(.4, .95), c(.2, .94), c(.6, .96),
    c(1, 4), c(.7, .5), .6, .4, .8, .3, c(.1, .9), c(0, 1))
  g <- as.list(grid::grid.grab()$children)
  for (i in 1:2) {
    box <- g[[2 + 3 * (i - 1)]]
    ci <- g[[3 + 3 * (i - 1)]]
    point <- g[[4 + 3 * (i - 1)]]
    testthat::expect_s3_class(box, "rect")
    testthat::expect_s3_class(ci, "segments")
    testthat::expect_equal(as.numeric(ci$x0), .1 + .8 * c(.2, .94)[i])
    testthat::expect_equal(as.numeric(ci$x1), .1 + .8 * c(.6, .96)[i])
    testthat::expect_equal(as.numeric(point$x0), .1 + .8 * c(.4, .95)[i])
    testthat::expect_equal(point$x0, point$x1)
  }
  testthat::expect_equal(as.numeric(g[[2]]$width)^2 / as.numeric(g[[5]]$width)^2, .25)

  se$w.random <- rep(NA, 3)
  testthat::expect_warning(plot_sensspec_forest_meta(se, sp), "sensitivity: random-effects study weights unavailable")
  rects <- Filter(function(g) inherits(g, "rect"), as.list(grid::grid.grab()$children))
  testthat::expect_equal(unname(vapply(rects[1:3], function(g) as.numeric(g$width), numeric(1))), rep(5, 3))
  for (bad in list(c(1, NA, 3), c(-1, 2, 3), c(0, 0, 0), c(1, 2))) {
    testthat::expect_error(plot_sensspec_forest(d, study_weights = list(sensitivity = bad, specificity = rep(1, 3))), "weights must")
  }
})
