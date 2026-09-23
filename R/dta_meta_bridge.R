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
  study_weights <- list(sensitivity = sensitivity_meta$w.random,
    specificity = specificity_meta$w.random[match(data$study, specificity_meta$studlab)])
  plot_sensspec_forest(data, summary_override = summary_override, heterogeneity = heterogeneity, study_weights = study_weights, ...)
}

.draw_forest_panel <- function(values, lower, upper, weights, y, summary, summary_lower, summary_upper, summary_y, region, xlim, square_col = "grey70", square_max_mm = 5, ci_col = "black", ci_lwd = 1.2, diamond_col = "#C00000") {
  x <- function(value) region[1] + (value - xlim[1]) / diff(xlim) * diff(region)
  grid::grid.segments(x0 = grid::unit(x(summary), "npc"), x1 = grid::unit(x(summary), "npc"), y0 = grid::unit(min(y) - .07, "npc"), y1 = grid::unit(max(y) + .02, "npc"), gp = grid::gpar(lty = 3, lwd = 1.2, col = "black"))
  square_size <- square_max_mm * sqrt(weights / max(weights))
  for (i in seq_along(y)) {
    grid::grid.rect(x = grid::unit(x(values[i]), "npc"), y = grid::unit(y[i], "npc"), width = grid::unit(square_size[i], "mm"), height = grid::unit(square_size[i], "mm"), gp = grid::gpar(fill = square_col, col = NA))
    grid::grid.segments(x0 = grid::unit(x(lower[i]), "npc"), x1 = grid::unit(x(upper[i]), "npc"), y0 = grid::unit(y[i], "npc"), y1 = grid::unit(y[i], "npc"), gp = grid::gpar(lwd = ci_lwd, col = ci_col))
    grid::grid.segments(x0 = grid::unit(x(values[i]), "npc"), x1 = grid::unit(x(values[i]), "npc"), y0 = grid::unit(y[i], "npc") - grid::unit(1, "mm"), y1 = grid::unit(y[i], "npc") + grid::unit(1, "mm"), gp = grid::gpar(col = ci_col, lwd = ci_lwd))
  }
  grid::grid.polygon(
    x = grid::unit(x(c(summary_lower, summary, summary_upper, summary)), "npc"),
    y = grid::unit(c(summary_y, summary_y + .022, summary_y, summary_y - .022), "npc"),
    gp = grid::gpar(fill = diamond_col, col = diamond_col)
  )
}

#' Draw a two-panel forest plot for sensitivity and specificity.
#'
#' @param summary_override Named list containing sens/sens_lwr/sens_upr and/or
#'   spec/spec_lwr/spec_upr, usually from fit_forest_summary().
#' @return Invisibly, the supplied output filename (or NULL when drawn to the active device).
plot_sensspec_forest <- function(data, output_file = NULL, study_col = "study", summary_override = NULL, heterogeneity = NULL, sens_axis = seq(0, 1, .2), spec_axis = seq(0, 1, .2), column_widths = c(study = 2.8, tp = .5, fp = .5, fn = .5, tn = .5, sens_text = 1.6, spec_text = 1.6, sens_plot = 1.2, spec_plot = 1.2), width = 10, height = NULL, res = 300, heterogeneity_cex = .66, heterogeneity_x = .015, heterogeneity_y = NULL, font_family = "serif", study_weights = NULL,
  header_cex = .86, study_cex = .72, ci_text_cex = .68,
  summary_cex = .78, summary_count_cex = .75, summary_ci_cex = .7,
  axis_cex = .68, row_gap = NULL, square_col = "grey70", square_max_mm = 5,
  ci_col = "black", ci_lwd = 1.2, diamond_col = "#C00000") {
  forest <- .forest_data(data, study_col, summary_override = summary_override)
  studies <- forest$studies
  summary <- forest$summary
  n <- nrow(studies)
  required_columns <- c("study", "tp", "fp", "fn", "tn", "sens_text", "spec_text", "sens_plot", "spec_plot")
  if (!is.numeric(column_widths) || !identical(names(column_widths), required_columns) || any(!is.finite(column_widths) | column_widths <= 0)) {
    stop("column_widths must be a positive named vector with study, tp, fp, fn, tn, sens_text, spec_text, sens_plot and spec_plot.", call. = FALSE)
  }
  if (is.null(height)) height <- max(3.8, 1 + .28 * n)
  if (!is.character(font_family) || length(font_family) != 1L || is.na(font_family) || !nzchar(trimws(font_family))) {
    stop("font_family must be one non-empty font family name supported by the graphics device.", call. = FALSE)
  }
  if (!is.numeric(heterogeneity_cex) || length(heterogeneity_cex) != 1L || !is.finite(heterogeneity_cex) || heterogeneity_cex <= 0) {
    stop("heterogeneity_cex must be one positive finite number.", call. = FALSE)
  }
  if (!is.numeric(heterogeneity_x) || length(heterogeneity_x) != 1L || !is.finite(heterogeneity_x) || heterogeneity_x < 0 || heterogeneity_x > 1) {
    stop("heterogeneity_x must be one finite number between 0 and 1.", call. = FALSE)
  }
  if (!is.null(heterogeneity_y) && (!is.numeric(heterogeneity_y) || length(heterogeneity_y) != 2L || any(!is.finite(heterogeneity_y)) || any(heterogeneity_y < 0 | heterogeneity_y > 1))) {
    stop("heterogeneity_y must be NULL or two finite numbers between 0 and 1 (sensitivity, specificity).", call. = FALSE)
  }
  sizes <- list(header_cex = header_cex, study_cex = study_cex, ci_text_cex = ci_text_cex,
    summary_cex = summary_cex, summary_count_cex = summary_count_cex, summary_ci_cex = summary_ci_cex,
    axis_cex = axis_cex, square_max_mm = square_max_mm, ci_lwd = ci_lwd)
  for (name in names(sizes)) {
    value <- sizes[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0)
      stop(name, " must be one positive finite number.", call. = FALSE)
  }
  colors <- list(square_col = square_col, ci_col = ci_col, diamond_col = diamond_col)
  for (name in names(colors)) {
    value <- colors[[name]]
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        inherits(try(grDevices::col2rgb(value), silent = TRUE), "try-error"))
      stop(name, " must be one valid R color.", call. = FALSE)
  }
  if (is.null(row_gap)) row_gap <- min(.055, .52 / max(1, n - 1))
  if (!is.numeric(row_gap) || length(row_gap) != 1L || !is.finite(row_gap) || row_gap <= 0 || row_gap * (n - 1) > .52 + 1e-12)
    stop("row_gap must be positive and fit the study rows within 0.52 of the plot height; reduce row_gap for more studies.", call. = FALSE)
  if (!is.null(study_weights) && (!is.list(study_weights) || !all(c("sensitivity", "specificity") %in% names(study_weights)))) {
    stop("study_weights must be NULL or a list with sensitivity and specificity weights.", call. = FALSE)
  }
  panel_weights <- lapply(c("sensitivity", "specificity"), function(label) {
    w <- study_weights[[label]]
    if (is.null(w) || (is.atomic(w) && length(w) == n && all(is.na(w)))) {
      warning(label, ": random-effects study weights unavailable (e.g. GLMM); using equal-size squares, not sample sizes or estimated weights.", call. = FALSE)
      return(rep(1, n))
    }
    if (!is.numeric(w) || length(w) != n || any(!is.finite(w)) || any(w < 0) || !any(w > 0)) {
      stop(label, ": weights must be finite, non-negative, aligned to all studies, with at least one positive value.", call. = FALSE)
    }
    w
  })
  if (!is.null(output_file)) {
    extension <- tolower(tools::file_ext(output_file))
    if (extension == "png") grDevices::png(output_file, width = width, height = height, units = "in", res = res)
    else if (extension %in% c("tif", "tiff")) grDevices::tiff(output_file, width = width, height = height, units = "in", res = res, compression = "lzw")
    else stop("output_file must end in .png, .tif, or .tiff.", call. = FALSE)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  y <- seq(.845, .845 - row_gap * (n - 1), length.out = n)
  summary_y <- min(y) - .09
  if (is.null(heterogeneity_y)) heterogeneity_y <- c(summary_y - .075, summary_y - .115)
  text <- function(x, y, label, gp = grid::gpar(), ...) {
    # ponytail: apply the common family here so no label can retain a different font.
    gp$fontfamily <- font_family
    grid::grid.text(label, x = grid::unit(x, "npc"), y = grid::unit(y, "npc"), gp = gp, ...)
  }
  ci_label <- function(est, lwr, upr) sprintf("%.2f [%.2f; %.2f]", est, lwr, upr)
  format_heterogeneity <- function(x, label) {
    if (is.null(x)) return(NULL)
    p <- x$p
    p_label <- if (!is.finite(p)) "NA" else if (p < .001) "<0.001" else sprintf("%.3f", p)
    sprintf("Heterogeneity for %s: \u03c4\u00b2 = %.2f, Q = %.2f, p %s, I\u00b2 = %.1f%%", label, x$tau2, x$q, if (p_label == "<0.001") "< 0.001" else paste0("= ", p_label), 100 * x$i2)
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
  for (i in seq_along(headers)) text(header_x[i], .91, headers[i], gp = grid::gpar(fontface = "bold", fontfamily = "serif", cex = header_cex), just = if (i == 1) "left" else "centre")
  grid::grid.segments(x0 = grid::unit(.015, "npc"), x1 = grid::unit(.985, "npc"), y0 = grid::unit(.875, "npc"), y1 = grid::unit(.875, "npc"), gp = grid::gpar(lwd = 2.2))
  for (i in seq_len(n)) {
    row <- studies[i, ]
    text(starts[["study"]], y[i], row$study, just = "left", gp = grid::gpar(cex = study_cex))
    for (j in seq_along(c("TP", "FP", "FN", "TN"))) text(centers[[tolower(c("TP", "FP", "FN", "TN")[j])]], y[i], row[[c("TP", "FP", "FN", "TN")[j]]], gp = grid::gpar(cex = study_cex))
    text(centers[["sens_text"]], y[i], ci_label(row$sens, row$sens_lwr, row$sens_upr), gp = grid::gpar(cex = ci_text_cex))
    text(centers[["spec_text"]], y[i], ci_label(row$spec, row$spec_lwr, row$spec_upr), gp = grid::gpar(cex = ci_text_cex))
  }
  text(starts[["study"]], summary_y, summary$study, just = "left", gp = grid::gpar(fontface = "bold", fontfamily = "serif", cex = summary_cex))
  for (j in seq_along(c("TP", "FP", "FN", "TN"))) text(centers[[tolower(c("TP", "FP", "FN", "TN")[j])]], summary_y, summary[[c("TP", "FP", "FN", "TN")[j]]], gp = grid::gpar(fontface = "bold", fontfamily = "serif", cex = summary_count_cex))
  text(centers[["sens_text"]], summary_y, ci_label(summary$sens, summary$sens_lwr, summary$sens_upr), gp = grid::gpar(fontface = "bold", cex = summary_ci_cex))
  text(centers[["spec_text"]], summary_y, ci_label(summary$spec, summary$spec_lwr, summary$spec_upr), gp = grid::gpar(fontface = "bold", cex = summary_ci_cex))
  sens_region <- c(starts[["sens_plot"]] + .004, ends[["sens_plot"]] - .012)
  spec_region <- c(starts[["spec_plot"]] + .012, ends[["spec_plot"]] - .004)
  .draw_forest_panel(studies$sens, studies$sens_lwr, studies$sens_upr, panel_weights[[1]], y, summary$sens, summary$sens_lwr, summary$sens_upr, summary_y, sens_region, range(sens_axis), square_col, square_max_mm, ci_col, ci_lwd, diamond_col)
  .draw_forest_panel(studies$spec, studies$spec_lwr, studies$spec_upr, panel_weights[[2]], y, summary$spec, summary$spec_lwr, summary$spec_upr, summary_y, spec_region, range(spec_axis), square_col, square_max_mm, ci_col, ci_lwd, diamond_col)
  heterogeneity_label <- c(format_heterogeneity(heterogeneity$sensitivity, "sensitivity"), format_heterogeneity(heterogeneity$specificity, "specificity"))
  if (!is.null(heterogeneity)) for (i in seq_along(heterogeneity_label)) if (!is.null(heterogeneity_label[i])) text(heterogeneity_x, heterogeneity_y[i], heterogeneity_label[i], just = "left", gp = grid::gpar(fontfamily = "serif", cex = heterogeneity_cex, col = "black"))
  for (axis in list(list(ticks = sens_axis, region = sens_region), list(ticks = spec_axis, region = spec_region))) {
    x <- axis$region[1] + (axis$ticks - min(axis$ticks)) / diff(range(axis$ticks)) * diff(axis$region)
    axis_y <- summary_y - .05
    grid::grid.segments(x0 = grid::unit(axis$region[1], "npc"), x1 = grid::unit(axis$region[2], "npc"), y0 = grid::unit(axis_y, "npc"), y1 = grid::unit(axis_y, "npc"), gp = grid::gpar(lwd = 1.8))
    for (i in seq_along(x)) {
      grid::grid.segments(x0 = grid::unit(x[i], "npc"), x1 = grid::unit(x[i], "npc"), y0 = grid::unit(axis_y, "npc"), y1 = grid::unit(axis_y - .025, "npc"), gp = grid::gpar(lwd = 1.8))
      text(x[i], axis_y - .055, formatC(axis$ticks[i], format = "f", digits = 1), gp = grid::gpar(fontfamily = "serif", cex = axis_cex))
    }
  }
  invisible(output_file)
}

#' Fit a bivariate diagnostic model and generate SROC plot data.
#'
#' By default, the SROC uses Reitsma REML and the Rutter-Gatsonis curve. The
#' conditional-mean ("naive") curve is available explicitly, but can run in a
#' non-ROC direction when the study-level covariance is negative.
fit_bivariate_dta <- function(data, study_col = "study", correction = 0.5, correction_control = c("single", "all", "none"), method = c("reml", "ml", "fixed"), sroc_type = c("ruttergatsonis", "qmd", "midas", "naive"), n_grid = 1000, seed = 2026, n_mc = 3000) {
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
  if (any(!is.finite(fit$Psi)) || any(diag(fit$Psi) <= 0)) stop("SROC requires positive between-study variances; inspect the fitted model or use a different model.", call. = FALSE)
  observed_fpr <- standardized$FP / (standardized$FP + standardized$TN)
  display_fpr_grid <- seq(max(.001, min(observed_fpr)), min(.999, max(observed_fpr)), length.out = n_grid)
  standard_sroc <- mada::sroc(fit, fpr = display_fpr_grid, type = sroc_type)
  curve_function <- mada::sroc(fit, type = sroc_type, return_function = TRUE)
  auc <- stats::integrate(curve_function, lower = 0, upper = 1, rel.tol = 1e-8)$value
  observed_range <- range(observed_fpr)
  pauc <- if (diff(observed_range) == 0) 0 else stats::integrate(curve_function, observed_range[1], observed_range[2], rel.tol = 1e-8)$value
  list(
    model = fit,
    model_type = paste("mada::reitsma", sroc_type, "SROC"),
    backend = "mada",
    metrics = list(
      sensitivity = sensitivity,
      specificity = specificity,
      auc = c(est = auc, lwr = NA_real_, upr = NA_real_),
      pauc = pauc
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
plot_sroc <- function(fit, show_confidence = TRUE, show_prediction = TRUE, output_file = NULL,
  width = 6, height = 5, dpi = 300, ..., full_curve = FALSE,
  show_legend = TRUE, digits = 2, custom_se = NULL, custom_sp = NULL, custom_auc = NULL,
  legend_position = "bottomright", legend_text_size = 3.5, legend_bg = "#F8F9FA",
  font_family = "sans", base_size = 14, title = "SROC with Prediction & Confidence Contours",
  show_study_labels = TRUE, study_size = 4, study_label_size = 2.5,
  summary_size = 4.5, summary_col = "#C0392B", sroc_col = "#2C3E50", sroc_linewidth = 1.2,
  confidence_col = "#2980B9", confidence_alpha = .2,
  prediction_col = "#BDC3C7", prediction_alpha = .15,
  auc_digits = NULL, x_breaks = seq(0, 1, .2), y_breaks = seq(0, 1, .2)) {
  if (length(list(...))) stop("Unknown plot_sroc arguments: ", paste(names(list(...)), collapse = ", "), call. = FALSE)
  if (!is.logical(full_curve) || length(full_curve) != 1L || is.na(full_curve))
    stop("full_curve must be TRUE or FALSE.", call. = FALSE)
  if (full_curve && !identical(fit$model_type, "mada::reitsma ruttergatsonis SROC"))
    stop("full_curve = TRUE currently supports the ruttergatsonis model only.", call. = FALSE)
  is_qmd <- identical(fit$backend, "qmd") || !is.null(fit$plot_data$sroc_df)
  if (is.null(auc_digits)) auc_digits <- if (is_qmd) digits else 3
  if (is_qmd) {
    adapted <- fit
  } else {
    pd <- fit$plot_data
    if (full_curve) {
      fpr <- seq(0, 1, length.out = max(2001L, nrow(pd$sroc)))
      curve <- mada::sroc(fit$model, type = "ruttergatsonis", return_function = TRUE)
      pd$sroc <- data.frame(sp = 1 - fpr, se = curve(fpr))
    }
    # ponytail: share one renderer across backends so every style control works.
    adapted <- list(metrics = list(se = fit$metrics$sensitivity, sp = fit$metrics$specificity, auc = fit$metrics$auc),
      plot_data = list(conf_df = pd$confidence, pred_df = pd$prediction, sroc_df = pd$sroc,
        study_df = data.frame(Study = pd$studies$Study, study_id = seq_len(nrow(pd$studies)),
          se = pd$studies$sensitivity, sp = pd$studies$specificity)))
  }
  p <- plot_metandi(adapted, digits = digits, show_conf = show_confidence, show_pred = show_prediction,
    show_legend = show_legend, custom_se = custom_se, custom_sp = custom_sp, custom_auc = custom_auc,
    legend_position = legend_position, legend_text_size = legend_text_size, legend_bg = legend_bg,
    font_family = font_family, base_size = base_size, title = title,
    show_study_labels = show_study_labels, study_size = study_size, study_label_size = study_label_size,
    summary_size = summary_size, summary_col = summary_col, sroc_col = sroc_col, sroc_linewidth = sroc_linewidth,
    confidence_col = confidence_col, confidence_alpha = confidence_alpha,
    prediction_col = prediction_col, prediction_alpha = prediction_alpha,
    auc_digits = auc_digits, x_breaks = x_breaks, y_breaks = y_breaks)
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
