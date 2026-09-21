# Diagnostic-test meta-analysis helpers.
#
# Required packages by feature:
# - forest plots: grid, meta
# - bivariate SROC and post-test probabilities: dplyr, tidyr, lme4, MASS, ggplot2

.assert_dta_data <- function(data, study_col = "study") {
  required <- c(study_col, "TP", "FP", "FN", "TN")
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  counts <- as.matrix(data[, c("TP", "FP", "FN", "TN")])
  if (!is.numeric(counts) || any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts))) {
    stop("TP, FP, FN and TN must be finite, non-negative integer counts.", call. = FALSE)
  }
  if (any(data$TP + data$FN == 0) || any(data$TN + data$FP == 0)) {
    stop("Every study must contain at least one diseased and one non-diseased participant.", call. = FALSE)
  }
  invisible(data)
}

.binom_ci <- function(success, total, conf_level = 0.95) {
  ci <- stats::binom.test(success, total, conf.level = conf_level)$conf.int
  c(est = success / total, lwr = ci[1], upr = ci[2])
}

.forest_data <- function(data, study_col = "study", conf_level = 0.95, summary_override = NULL) {
  .assert_dta_data(data, study_col)
  studies <- data.frame(
    study = as.character(data[[study_col]]), TP = data$TP, FP = data$FP, FN = data$FN, TN = data$TN
  )
  sens <- t(vapply(seq_len(nrow(studies)), function(i) .binom_ci(studies$TP[i], studies$TP[i] + studies$FN[i], conf_level), numeric(3)))
  spec <- t(vapply(seq_len(nrow(studies)), function(i) .binom_ci(studies$TN[i], studies$TN[i] + studies$FP[i], conf_level), numeric(3)))
  studies <- cbind(studies, setNames(as.data.frame(sens), c("sens", "sens_lwr", "sens_upr")), setNames(as.data.frame(spec), c("spec", "spec_lwr", "spec_upr")))

  total <- colSums(studies[, c("TP", "FP", "FN", "TN")])
  total_sens <- .binom_ci(total[["TP"]], total[["TP"]] + total[["FN"]], conf_level)
  total_spec <- .binom_ci(total[["TN"]], total[["TN"]] + total[["FP"]], conf_level)
  summary <- data.frame(study = "Total (95% CI)", as.list(total), sens = total_sens[["est"]], sens_lwr = total_sens[["lwr"]], sens_upr = total_sens[["upr"]], spec = total_spec[["est"]], spec_lwr = total_spec[["lwr"]], spec_upr = total_spec[["upr"]])
  if (!is.null(summary_override)) {
    allowed <- intersect(names(summary_override), names(summary))
    unknown <- setdiff(names(summary_override), names(summary))
    if (length(unknown)) stop("Unknown summary fields: ", paste(unknown, collapse = ", "), call. = FALSE)
    for (name in allowed) summary[[name]] <- summary_override[[name]]
  }
  list(studies = studies, summary = summary)
}

#' Fit separate random-effects GLMMs for a sensitivity/specificity forest plot.
#'
#' The summary estimates are intentionally univariate. Use fit_bivariate_dta()
#' when the sensitivity-specificity correlation is part of the estimand.
fit_forest_summary <- function(data, study_col = "study") {
  .assert_dta_data(data, study_col)
  fit_prop <- function(event, n) meta::metaprop(event, n, studlab = data[[study_col]], method = "GLMM", method.tau = "ML", method.I2 = "Q")
  sens <- fit_prop(data$TP, data$TP + data$FN)
  spec <- fit_prop(data$TN, data$TN + data$FP)
  unpack <- function(x) c(est = meta::backtransf(x$TE.random, sm = x$sm), lwr = meta::backtransf(x$lower.random, sm = x$sm), upr = meta::backtransf(x$upper.random, sm = x$sm))
  heterogeneity <- function(x) {
    p <- x$pval.Q
    if (length(p) > 1 && "LRT" %in% names(p)) p <- p[["LRT"]] else p <- p[1]
    q <- x$Q
    if (length(q) > 1 && "LRT" %in% names(q)) q <- q[["LRT"]] else q <- q[1]
    list(i2 = x$I2, tau2 = x$tau2, q = q, p = p)
  }
  list(sensitivity = unpack(sens), specificity = unpack(spec), models = list(sensitivity = sens, specificity = spec), heterogeneity = list(sensitivity = heterogeneity(sens), specificity = heterogeneity(spec)))
}

.assert_metaprop <- function(object, label) {
  required <- c("event", "n", "studlab", "TE.random", "lower.random", "upper.random", "sm")
  if (!inherits(object, "meta") || !all(required %in% names(object))) {
    stop(label, " must be a meta::metaprop object with random-effects results.", call. = FALSE)
  }
  if (anyDuplicated(as.character(object$studlab))) {
    stop(label, " contains duplicate study labels; align the studies before using this adapter.", call. = FALSE)
  }
  invisible(object)
}

#' Reconstruct a diagnostic 2x2 data set from two meta::metaprop objects.
#'
#' @param sensitivity_meta A metaprop object with TP as `event` and TP + FN as `n`.
#' @param specificity_meta A metaprop object with TN as `event` and TN + FP as `n`.
#' @return A data frame with study, TP, FP, FN and TN columns.
dta_from_meta <- function(sensitivity_meta, specificity_meta) {
  .assert_metaprop(sensitivity_meta, "sensitivity_meta")
  .assert_metaprop(specificity_meta, "specificity_meta")
  sens_studies <- as.character(sensitivity_meta$studlab)
  spec_index <- match(sens_studies, as.character(specificity_meta$studlab))
  if (anyNA(spec_index) || length(spec_index) != length(specificity_meta$studlab)) {
    stop("The two meta objects must contain exactly the same uniquely labelled studies.", call. = FALSE)
  }
  data <- data.frame(
    study = sens_studies,
    TP = sensitivity_meta$event,
    FN = sensitivity_meta$n - sensitivity_meta$event,
    TN = specificity_meta$event[spec_index],
    FP = specificity_meta$n[spec_index] - specificity_meta$event[spec_index]
  )
  .assert_dta_data(data)
  data
}

.meta_summary_override <- function(sensitivity_meta, specificity_meta) {
  extract <- function(object, label) {
    result <- meta::backtransf(c(object$TE.random, object$lower.random, object$upper.random), sm = object$sm)
    if (any(!is.finite(result))) stop(label, " has no finite random-effects summary; refit it with random effects enabled.", call. = FALSE)
    result
  }
  sens <- extract(sensitivity_meta, "sensitivity_meta")
  spec <- extract(specificity_meta, "specificity_meta")
  list(sens = sens[1], sens_lwr = sens[2], sens_upr = sens[3], spec = spec[1], spec_lwr = spec[2], spec_upr = spec[3])
}

.meta_heterogeneity <- function(object) {
  p <- object$pval.Q
  if (length(p) > 1 && "LRT" %in% names(p)) p <- p[["LRT"]] else p <- p[1]
  q <- object$Q
  if (length(q) > 1 && "LRT" %in% names(q)) q <- q[["LRT"]] else q <- q[1]
  list(i2 = object$I2, tau2 = object$tau2, q = q, p = p)
}

#' Draw a two-panel forest plot directly from meta::metaprop objects.
#'
#' By default, the diamond is the exact random-effects summary stored in the
#' supplied meta objects; study rows are reconstructed from their event/n data.
plot_sensspec_forest_meta <- function(sensitivity_meta, specificity_meta, ..., use_meta_summary = TRUE) {
  data <- dta_from_meta(sensitivity_meta, specificity_meta)
  summary_override <- if (use_meta_summary) .meta_summary_override(sensitivity_meta, specificity_meta) else NULL
  heterogeneity <- list(sensitivity = .meta_heterogeneity(sensitivity_meta), specificity = .meta_heterogeneity(specificity_meta))
  plot_sensspec_forest(data, summary_override = summary_override, heterogeneity = heterogeneity, ...)
}

.draw_forest_panel <- function(values, lower, upper, weights, y, summary, summary_lower, summary_upper, summary_y, region, xlim) {
  x <- function(value) region[1] + (value - xlim[1]) / diff(xlim) * diff(region)
  grid::grid.segments(x0 = grid::unit(x(summary), "npc"), x1 = grid::unit(x(summary), "npc"), y0 = grid::unit(min(y) - .07, "npc"), y1 = grid::unit(max(y) + .04, "npc"), gp = grid::gpar(lty = 3, lwd = 1.2, col = "black"))
  square_size <- 1.8 + 3.2 * sqrt(weights / max(weights))
  for (i in seq_along(y)) {
    grid::grid.segments(x0 = grid::unit(x(lower[i]), "npc"), x1 = grid::unit(x(upper[i]), "npc"), y0 = grid::unit(y[i], "npc"), y1 = grid::unit(y[i], "npc"), gp = grid::gpar(lwd = 1.2, col = "black"))
    grid::grid.rect(x = grid::unit(x(values[i]), "npc"), y = grid::unit(y[i], "npc"), width = grid::unit(square_size[i], "mm"), height = grid::unit(square_size[i], "mm"), gp = grid::gpar(fill = "grey70", col = NA))
    grid::grid.points(x = grid::unit(x(values[i]), "npc"), y = grid::unit(y[i], "npc"), pch = 3, size = grid::unit(2.1, "mm"), gp = grid::gpar(col = "black", lwd = 1.1))
  }
  grid::grid.polygon(
    x = grid::unit(x(c(summary_lower, summary, summary_upper, summary)), "npc"),
    y = grid::unit(c(summary_y, summary_y + .022, summary_y, summary_y - .022), "npc"),
    gp = grid::gpar(fill = "#C00000", col = "#C00000")
  )
}

#' Draw a two-panel forest plot for sensitivity and specificity.
#'
#' @param summary_override Named list containing sens/sens_lwr/sens_upr and/or
#'   spec/spec_lwr/spec_upr, usually from fit_forest_summary().
#' @return Invisibly, the supplied output filename (or NULL when drawn to the active device).
plot_sensspec_forest <- function(data, output_file = NULL, study_col = "study", summary_override = NULL, heterogeneity = NULL, sens_axis = seq(0, 1, .2), spec_axis = seq(0, 1, .2), column_widths = c(study = 2.8, tp = .5, fp = .5, fn = .5, tn = .5, sens_text = 1.6, spec_text = 1.6, sens_plot = 1.2, spec_plot = 1.2), width = 10, height = NULL, res = 300) {
  forest <- .forest_data(data, study_col, summary_override = summary_override)
  studies <- forest$studies
  summary <- forest$summary
  n <- nrow(studies)
  required_columns <- c("study", "tp", "fp", "fn", "tn", "sens_text", "spec_text", "sens_plot", "spec_plot")
  if (!is.numeric(column_widths) || !identical(names(column_widths), required_columns) || any(!is.finite(column_widths) | column_widths <= 0)) {
    stop("column_widths must be a positive named vector with study, tp, fp, fn, tn, sens_text, spec_text, sens_plot and spec_plot.", call. = FALSE)
  }
  if (is.null(height)) height <- max(3.8, 1 + .28 * n)
  if (!is.null(output_file)) {
    extension <- tolower(tools::file_ext(output_file))
    if (extension == "png") grDevices::png(output_file, width = width, height = height, units = "in", res = res)
    else if (extension %in% c("tif", "tiff")) grDevices::tiff(output_file, width = width, height = height, units = "in", res = res, compression = "lzw")
    else stop("output_file must end in .png, .tif, or .tiff.", call. = FALSE)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  row_gap <- min(.055, .52 / max(1, n - 1))
  y <- seq(.79, .79 - row_gap * (n - 1), length.out = n)
  summary_y <- min(y) - .09
  heterogeneity_y <- c(summary_y - .075, summary_y - .115)
  text <- function(x, y, label, ...) grid::grid.text(label, x = grid::unit(x, "npc"), y = grid::unit(y, "npc"), ...)
  ci_label <- function(est, lwr, upr) sprintf("%.2f [%.2f; %.2f]", est, lwr, upr)
  format_heterogeneity <- function(x, label) {
    if (is.null(x)) return(NULL)
    p <- x$p
    p_label <- if (!is.finite(p)) "NA" else if (p < .001) "<0.001" else sprintf("%.3f", p)
    sprintf("Heterogeneity for %s: tau2 = %.2f, Q = %.2f, p %s, I2 = %.1f%%", label, x$tau2, x$q, if (p_label == "<0.001") "< 0.001" else paste0("= ", p_label), 100 * x$i2)
  }
  left <- .015
  widths <- .97 * column_widths / sum(column_widths)
  starts <- left + c(0, cumsum(widths)[-length(widths)])
  ends <- starts + widths
  centers <- (starts + ends) / 2
  names(starts) <- names(ends) <- names(centers) <- names(column_widths)
  grid::grid.newpage()
  headers <- c("Study", "TP", "FP", "FN", "TN", "Sensitivity", "Specificity", "Sensitivity", "Specificity")
  header_x <- c(starts[["study"]], centers[c("tp", "fp", "fn", "tn", "sens_text", "spec_text", "sens_plot", "spec_plot")])
  for (i in seq_along(headers)) text(header_x[i], .91, headers[i], gp = grid::gpar(fontface = "bold", fontfamily = "serif", cex = .86), just = if (i == 1) "left" else "centre")
  grid::grid.segments(x0 = grid::unit(.015, "npc"), x1 = grid::unit(.985, "npc"), y0 = grid::unit(.875, "npc"), y1 = grid::unit(.875, "npc"), gp = grid::gpar(lwd = 2.2))
  for (i in seq_len(n)) {
    row <- studies[i, ]
    text(starts[["study"]], y[i], row$study, just = "left", gp = grid::gpar(cex = .72))
    for (j in seq_along(c("TP", "FP", "FN", "TN"))) text(centers[[tolower(c("TP", "FP", "FN", "TN")[j])]], y[i], row[[c("TP", "FP", "FN", "TN")[j]]], gp = grid::gpar(cex = .72))
    text(centers[["sens_text"]], y[i], ci_label(row$sens, row$sens_lwr, row$sens_upr), gp = grid::gpar(cex = .68))
    text(centers[["spec_text"]], y[i], ci_label(row$spec, row$spec_lwr, row$spec_upr), gp = grid::gpar(cex = .68))
  }
  text(starts[["study"]], summary_y, summary$study, just = "left", gp = grid::gpar(fontface = "bold", fontfamily = "serif", cex = .78))
  for (j in seq_along(c("TP", "FP", "FN", "TN"))) text(centers[[tolower(c("TP", "FP", "FN", "TN")[j])]], summary_y, summary[[c("TP", "FP", "FN", "TN")[j]]], gp = grid::gpar(fontface = "bold", fontfamily = "serif", cex = .75))
  text(centers[["sens_text"]], summary_y, ci_label(summary$sens, summary$sens_lwr, summary$sens_upr), gp = grid::gpar(fontface = "bold", cex = .7))
  text(centers[["spec_text"]], summary_y, ci_label(summary$spec, summary$spec_lwr, summary$spec_upr), gp = grid::gpar(fontface = "bold", cex = .7))
  sens_region <- c(starts[["sens_plot"]] + .004, ends[["sens_plot"]] - .012)
  spec_region <- c(starts[["spec_plot"]] + .012, ends[["spec_plot"]] - .004)
  .draw_forest_panel(studies$sens, studies$sens_lwr, studies$sens_upr, studies$TP + studies$FN, y, summary$sens, summary$sens_lwr, summary$sens_upr, summary_y, sens_region, range(sens_axis))
  .draw_forest_panel(studies$spec, studies$spec_lwr, studies$spec_upr, studies$TN + studies$FP, y, summary$spec, summary$spec_lwr, summary$spec_upr, summary_y, spec_region, range(spec_axis))
  heterogeneity_label <- c(format_heterogeneity(heterogeneity$sensitivity, "sensitivity"), format_heterogeneity(heterogeneity$specificity, "specificity"))
  if (!is.null(heterogeneity)) for (i in seq_along(heterogeneity_label)) if (!is.null(heterogeneity_label[i])) text(starts[["study"]], heterogeneity_y[i], heterogeneity_label[i], just = "left", gp = grid::gpar(fontfamily = "serif", cex = .66, col = "black"))
  for (axis in list(list(ticks = sens_axis, region = sens_region), list(ticks = spec_axis, region = spec_region))) {
    x <- axis$region[1] + (axis$ticks - min(axis$ticks)) / diff(range(axis$ticks)) * diff(axis$region)
    axis_y <- summary_y - .05
    grid::grid.segments(x0 = grid::unit(axis$region[1], "npc"), x1 = grid::unit(axis$region[2], "npc"), y0 = grid::unit(axis_y, "npc"), y1 = grid::unit(axis_y, "npc"), gp = grid::gpar(lwd = 1.8))
    for (i in seq_along(x)) {
      grid::grid.segments(x0 = grid::unit(x[i], "npc"), x1 = grid::unit(x[i], "npc"), y0 = grid::unit(axis_y, "npc"), y1 = grid::unit(axis_y - .025, "npc"), gp = grid::gpar(lwd = 1.8))
      text(x[i], axis_y - .055, formatC(axis$ticks[i], format = "f", digits = 1), gp = grid::gpar(fontfamily = "serif", cex = .68))
    }
  }
  invisible(output_file)
}

#' Fit a bivariate diagnostic model and generate SROC plot data.
#'
#' By default, the SROC uses the Midas binomial model and curve formula. The
#' conditional-mean ("naive") curve is available explicitly, but can run in a
#' non-ROC direction when the study-level covariance is negative.
fit_bivariate_dta <- function(data, study_col = "study", correction = 0.5, correction_control = c("single", "all", "none"), method = c("reml", "ml", "fixed"), sroc_type = c("qmd", "midas", "ruttergatsonis", "naive"), n_grid = 1000, seed = 2026, n_mc = 3000) {
  .assert_dta_data(data, study_col)
  if (nrow(data) < 3) stop("At least three studies are required for the bivariate model.", call. = FALSE)
  sroc_type <- match.arg(sroc_type)
  if (sroc_type == "qmd") {
    if (!missing(method) || !missing(correction) || !missing(correction_control)) stop("method and correction arguments apply only to mada backends.", call. = FALSE)
    standardized <- data.frame(Study = data[[study_col]], data[, c("TP", "FP", "FN", "TN")])
    result <- fit_metandi(standardized, seed = seed, n_mc = n_mc, n_grid = n_grid)
    result$backend <- "qmd"
    result$model_type <- "Original QMD fit_metandi"
    result$metrics$sensitivity <- result$metrics$se
    result$metrics$specificity <- result$metrics$sp
    return(result)
  }
  if (sroc_type == "midas") {
    if (!missing(method) || !missing(correction) || !missing(correction_control)) stop("method and continuity-correction arguments apply only to the mada backends; omit them for midas.", call. = FALSE)
    return(.fit_midas(data, study_col, n_grid))
  }
  correction_control <- match.arg(correction_control)
  method <- match.arg(method)
  if (!is.numeric(correction) || length(correction) != 1 || correction < 0) stop("correction must be a single non-negative number.", call. = FALSE)
  if (n_grid < 100) stop("n_grid must be at least 100.", call. = FALSE)
  standardized <- data.frame(study = as.character(data[[study_col]]), TP = data$TP, FP = data$FP, FN = data$FN, TN = data$TN)
  fit <- mada::reitsma(standardized, TP = "TP", FN = "FN", FP = "FP", TN = "TN", correction = correction, correction.control = correction_control, method = method)
  mu <- stats::coef(fit)["(Intercept)", ]
  fixed_vcov <- as.matrix(stats::vcov(fit))
  z <- stats::qnorm(.975)
  se_fixed <- sqrt(diag(fixed_vcov))
  sensitivity <- c(est = unname(stats::plogis(mu[1])), lwr = unname(stats::plogis(mu[1] - z * se_fixed[1])), upr = unname(stats::plogis(mu[1] + z * se_fixed[1])))
  specificity <- c(est = unname(1 - stats::plogis(mu[2])), lwr = unname(1 - stats::plogis(mu[2] + z * se_fixed[2])), upr = unname(1 - stats::plogis(mu[2] - z * se_fixed[2])))
  to_roc_ellipse <- function(covariance) {
    points <- ellipse::ellipse(covariance, centre = mu, level = .95)
    data.frame(sp = 1 - stats::plogis(points[, 2]), se = stats::plogis(points[, 1]))
  }
  full_fpr_grid <- seq(.001, .999, length.out = n_grid)
  observed_fpr <- standardized$FP / (standardized$FP + standardized$TN)
  display_fpr_grid <- seq(max(.001, min(observed_fpr)), min(.999, max(observed_fpr)), length.out = n_grid)
  standard_sroc <- mada::sroc(fit, fpr = display_fpr_grid, type = sroc_type)
  auc_result <- mada::AUC(fit, fpr = full_fpr_grid, sroc.type = sroc_type)
  list(
    model = fit,
    model_type = paste("mada::reitsma", sroc_type, "SROC"),
    metrics = list(
      sensitivity = sensitivity,
      specificity = specificity,
      auc = c(est = auc_result[["AUC"]], lwr = NA_real_, upr = NA_real_),
      pauc = unname(auc_result[["pAUC"]])
    ),
    plot_data = list(confidence = to_roc_ellipse(fixed_vcov), prediction = to_roc_ellipse(fixed_vcov + fit$Psi), sroc = data.frame(sp = 1 - standard_sroc[, 1], se = standard_sroc[, 2]), studies = data.frame(Study = standardized$study, specificity = data$TN / (data$TN + data$FP), sensitivity = data$TP / (data$TP + data$FN))),
    random_effects = list(covariance = fit$Psi)
  )
}

#' Fit the bivariate SROC model directly from a sensitivity/specificity metaprop pair.
fit_bivariate_meta <- function(sensitivity_meta, specificity_meta, ...) {
  fit_bivariate_dta(dta_from_meta(sensitivity_meta, specificity_meta), ...)
}

#' Plot an SROC curve with confidence and prediction contours.
#'
#' @return A ggplot object, invisibly saved to output_file when supplied.
plot_sroc <- function(fit, show_confidence = TRUE, show_prediction = TRUE, output_file = NULL, width = 6, height = 5, dpi = 300, ...) {
  if (identical(fit$backend, "qmd") || !is.null(fit$plot_data$sroc_df)) {
    p <- plot_metandi(fit, show_conf = show_confidence, show_pred = show_prediction, ...)
    if (!is.null(output_file)) ggplot2::ggsave(output_file, plot = p, width = width, height = height, dpi = dpi)
    return(p)
  }
  pd <- fit$plot_data; mt <- fit$metrics
  p <- ggplot2::ggplot()
  if (show_prediction) p <- p + ggplot2::geom_polygon(data = pd$prediction, ggplot2::aes(x = sp, y = se), fill = "grey70", alpha = .2) + ggplot2::geom_path(data = pd$prediction, ggplot2::aes(x = sp, y = se), linetype = "dotted", colour = "grey45")
  if (show_confidence) p <- p + ggplot2::geom_polygon(data = pd$confidence, ggplot2::aes(x = sp, y = se), fill = "#2980B9", alpha = .18) + ggplot2::geom_path(data = pd$confidence, ggplot2::aes(x = sp, y = se), linetype = "dashed", colour = "#2980B9")
  p <- p +
    ggplot2::geom_line(data = pd$sroc, ggplot2::aes(x = sp, y = se), linewidth = 1.1, colour = "#2C3E50") +
    ggplot2::geom_point(data = pd$studies, ggplot2::aes(x = specificity, y = sensitivity), shape = 21, fill = "white", colour = "#7F8C8D", size = 3) +
    ggplot2::geom_point(ggplot2::aes(x = mt$specificity[["est"]], y = mt$sensitivity[["est"]]), shape = 15, size = 4, colour = "#C0392B") +
    ggplot2::scale_x_reverse(limits = c(1, 0), breaks = seq(0, 1, .2)) + ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) + ggplot2::coord_fixed() + ggplot2::labs(x = "Specificity", y = "Sensitivity", title = "SROC") + ggplot2::theme_classic()
  if (!is.null(output_file)) ggplot2::ggsave(output_file, plot = p, width = width, height = height, dpi = dpi)
  p
}

#' Calculate post-test probabilities from the pooled Reitsma-model effects.
#'
#' These intervals quantify uncertainty in pooled mean accuracy, not prediction
#' intervals for a future setting.
posttest_probability <- function(fit, prevalence = .3, n_sims = 3000, seed = 2026) {
  if (!is.numeric(prevalence) || any(!is.finite(prevalence)) || any(prevalence <= 0 | prevalence >= 1)) stop("prevalence must contain values strictly between 0 and 1.", call. = FALSE)
  if (length(n_sims) != 1 || n_sims < 100 || n_sims != floor(n_sims)) stop("n_sims must be an integer of at least 100.", call. = FALSE)
  model <- fit$model
  if (inherits(model, "glmerMod")) {
    beta <- unname(lme4::fixef(model)) * c(1, -1)
    transform <- diag(c(1, -1))
    vcov_beta <- transform %*% as.matrix(stats::vcov(model)) %*% transform
  } else if (inherits(model, "reitsma")) {
    beta <- stats::coef(model)["(Intercept)", ]; vcov_beta <- as.matrix(stats::vcov(model))
  } else stop("fit must be returned by fit_bivariate_dta() or fit_bivariate_meta().", call. = FALSE)
  sens <- stats::plogis(beta[1]); spec <- 1 - stats::plogis(beta[2])
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({ if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old_seed, envir = .GlobalEnv) }, add = TRUE)
  set.seed(seed)
  sims <- MASS::mvrnorm(n_sims, mu = beta, Sigma = vcov_beta)
  sim_sens <- stats::plogis(sims[, 1]); sim_spec <- 1 - stats::plogis(sims[, 2])
  result <- do.call(rbind, lapply(prevalence, function(p) {
    odds <- p / (1 - p)
    ppv <- odds * sens / (1 - spec) / (1 + odds * sens / (1 - spec))
    npv <- 1 - odds * (1 - sens) / spec / (1 + odds * (1 - sens) / spec)
    sim_ppv <- odds * sim_sens / (1 - sim_spec); sim_ppv <- sim_ppv / (1 + sim_ppv)
    sim_npv <- odds * (1 - sim_sens) / sim_spec; sim_npv <- 1 - sim_npv / (1 + sim_npv)
    data.frame(prevalence = p, sensitivity = sens, specificity = spec, ppv = ppv, ppv_lwr = unname(stats::quantile(sim_ppv, .025)), ppv_upr = unname(stats::quantile(sim_ppv, .975)), npv = npv, npv_lwr = unname(stats::quantile(sim_npv, .025)), npv_upr = unname(stats::quantile(sim_npv, .975)))
  }))
  row.names(result) <- NULL
  result
}
