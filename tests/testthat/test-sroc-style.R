testthat::test_that("SROC style controls and legend anchors are honored", {
  d <- data.frame(study=LETTERS[1:5], TP=c(40,32,48,25,55), FP=c(8,15,10,20,12),
                  FN=c(10,18,7,15,9), TN=c(72,65,80,58,90))
  fit <- fit_bivariate_dta(d)
  original <- fit
  p <- plot_sroc(fit, full_curve=TRUE, legend_position=c(.03,.04),
    legend_text_size=3, legend_bg="ivory", font_family="serif", base_size=12,
    title="Style test", show_study_labels=FALSE, study_size=3, study_label_size=2,
    summary_size=5, summary_col="red", sroc_col="navy", sroc_linewidth=2,
    confidence_col="blue", confidence_alpha=.1, prediction_col="grey50",
    prediction_alpha=.3, auc_digits=4, x_breaks=c(0,.5,1), y_breaks=c(0,.25,.5,1),
    custom_auc=NULL)
  built <- ggplot2::ggplot_build(p)$data
  testthat::expect_identical(fit, original)
  testthat::expect_equal(p$labels$title,"Style test")
  testthat::expect_equal(p$theme$text$family,"serif")
  testthat::expect_equal(p$theme$text$size,12)
  testthat::expect_equal(p$scales$get_scales("x")$breaks,c(0,.5,1))
  testthat::expect_equal(p$scales$get_scales("y")$breaks,c(0,.25,.5,1))
  line <- Filter(function(x) inherits(x$geom,"GeomLine"),p$layers)[[1]]
  testthat::expect_equal(line$aes_params$colour,"navy")
  testthat::expect_equal(line$aes_params$linewidth,2)
  testthat::expect_false(any(vapply(p$layers,function(x) !is.null(x$mapping$label),logical(1))))
  testthat::expect_equal(built[[1]]$alpha[1],.3)
  testthat::expect_equal(built[[3]]$alpha[1],.1)
  points <- Filter(function(x) inherits(x$geom,"GeomPoint"),p$layers)
  testthat::expect_equal(points[[1]]$aes_params$size,3)
  testthat::expect_equal(points[[2]]$aes_params$size,5)
  testthat::expect_equal(points[[2]]$aes_params$colour,"red")
  rect <- Filter(function(x) inherits(x$geom,"GeomRect"),p$layers)[[1]]
  testthat::expect_equal(rect$aes_params$fill,"ivory")
  testthat::expect_equal(rect$data$xmin,.97)
  testthat::expect_equal(rect$data$xmax,.53)
  testthat::expect_equal(rect$data$ymin,.04)
  labels <- unlist(lapply(built,function(x) x$label))
  testthat::expect_true(any(grepl(sprintf("AUC: %.4f",fit$metrics$auc[1]),labels,fixed=TRUE)))
  testthat::expect_false(any(grepl("NA|full curve",labels)))
  text <- Filter(function(x) inherits(x$geom,"GeomText"),p$layers)
  testthat::expect_true(all(vapply(text,function(x) identical(x$aes_params$family,"serif"),logical(1))))
  testthat::expect_equal(text[[2]]$aes_params$size,3)
  for (corner in c("bottomright","bottomleft","topright","topleft")) {
    q <- plot_sroc(fit,legend_position=corner)
    box <- Filter(function(x) inherits(x$geom,"GeomRect"),q$layers)[[1]]$data
    testthat::expect_equal(box$xmin,if(grepl("right",corner)) .46 else .98)
    if(grepl("top",corner)) testthat::expect_equal(box$ymax,.98)
  }
  q <- plot_sroc(fit,show_legend=FALSE,show_study_labels=TRUE,title=NULL)
  testthat::expect_false(any(vapply(q$layers,function(x) inherits(x$geom,"GeomRect"),logical(1))))
  testthat::expect_true(any(vapply(q$layers,function(x) !is.null(x$mapping$label),logical(1))))
  for (args in list(list(legend_position="center"),list(legend_position=c(.9,.9)),
    list(legend_text_size=0),list(confidence_alpha=2),list(prediction_col="bad-color"),
    list(auc_digits=-1),list(x_breaks=c(0,2)),list(show_study_labels=NA),list(unknown=1))) {
    testthat::expect_error(do.call(plot_sroc,c(list(fit=fit),args)))
  }
  qfit <- fit_bivariate_dta(d,sroc_type="qmd",n_mc=100)
  qp <- plot_sroc(qfit,auc_digits=4,custom_auc=NULL,legend_position="topleft")
  qlabels <- unlist(lapply(ggplot2::ggplot_build(qp)$data,function(x) x$label))
  testthat::expect_true(any(grepl(sprintf("AUC: %.4f (",qfit$metrics$auc[1]),qlabels,fixed=TRUE)))
})
