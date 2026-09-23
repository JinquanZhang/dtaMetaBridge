# Extracted from HFphf_meta analysis.qmd: fit_metandi and plot_metandi.
# Changes: namespace qualification and base pipe only; original statistical formulas retained.
# ponytail: declare data-mask column names for package checks, without changing formulas.
utils::globalVariables(c("TP", "FP", "FN", "TN", "cell", "status", "count", "Study", "sp", "se", "study_id"))
fit_metandi <- function(data, seed = 2026, n_mc = 3000, n_grid = 1000) {
  required_cols <- c("TP", "FP", "FN", "TN")
  if (!all(required_cols %in% names(data))) {
    stop("Input data must contain columns: TP, FP, FN, TN.")
  }

  study_col <- dplyr::case_when(
    "Study" %in% names(data) ~ "Study",
    "study" %in% names(data) ~ "study",
    TRUE ~ NA_character_
  )
  if (is.na(study_col)) {
    stop("Input data must contain a study label column named `Study` or `study`.")
  }

  data_std <- data |>
    dplyr::rename(Study = !!rlang::sym(study_col))

  long_data <- data_std |>
    tidyr::pivot_longer(cols = c(TP, FP, FN, TN), names_to = "cell", values_to = "count") |>
    dplyr::mutate(
      type = dplyr::case_when(
        cell %in% c("TP", "FN") ~ "Sensitivity",
        cell %in% c("TN", "FP") ~ "Specificity"
      ),
      status = dplyr::case_when(
        cell %in% c("TP", "TN") ~ "Success",
        cell %in% c("FN", "FP") ~ "Failure"
      )
    ) |>
    dplyr::select(-cell) |>
    tidyr::pivot_wider(names_from = status, values_from = count) |>
    dplyr::mutate(Study = as.factor(Study))

  suppressMessages({
    glmm_fit <- lme4::glmer(
      cbind(Success, Failure) ~ 0 + type + (0 + type | Study),
      family = binomial,
      data = long_data
    )
  })

  fixef_est <- lme4::fixef(glmm_fit)
  mu_se <- unname(fixef_est["typeSensitivity"])
  mu_sp <- unname(fixef_est["typeSpecificity"])

  V <- stats::vcov(glmm_fit)
  s_se <- sqrt(V["typeSensitivity", "typeSensitivity"])
  s_sp <- sqrt(V["typeSpecificity", "typeSpecificity"])
  r_fixed <- V["typeSensitivity", "typeSpecificity"] / (s_se * s_sp)

  RE <- lme4::VarCorr(glmm_fit)
  re_study <- as.matrix(RE$Study)
  sigma_se <- attr(RE$Study, "stddev")["typeSensitivity"]
  sigma_sp <- attr(RE$Study, "stddev")["typeSpecificity"]
  rho_re <- attr(RE$Study, "correlation")[1, 2]

  V_pred <- V + re_study
  s_se_pred <- sqrt(V_pred[1, 1])
  s_sp_pred <- sqrt(V_pred[2, 2])
  r_pred <- V_pred[1, 2] / (s_se_pred * s_sp_pred)

  c_const <- sqrt(2 * stats::qf(0.95, df1 = 2, df2 = length(unique(long_data$Study)) - 2))
  z <- stats::qnorm(0.975)

  res <- list()
  res$se <- c(est = plogis(mu_se), lwr = plogis(mu_se - z * s_se), upr = plogis(mu_se + z * s_se))
  res$sp <- c(est = plogis(mu_sp), lwr = plogis(mu_sp - z * s_sp), upr = plogis(mu_sp + z * s_sp))

  set.seed(seed)
  sim_params <- MASS::mvrnorm(n = n_mc, mu = fixef_est, Sigma = V)
  sp_mc <- seq(1e-4, 1 - 1e-4, length.out = n_grid)
  logit_sp_mc <- qlogis(sp_mc)
  x_ord <- (1 - sp_mc)[order(1 - sp_mc)]

  auc_sims <- apply(sim_params, 1, function(p) {
    logit_se <- p[1] + (sigma_se / sigma_sp) * sign(rho_re) * (logit_sp_mc - p[2])
    y_ord <- plogis(logit_se)[order(1 - sp_mc)]
    sum(diff(x_ord) * (y_ord[-1] + y_ord[-length(y_ord)]) / 2)
  })
  res$auc <- c(est = mean(auc_sims), lwr = quantile(auc_sims, 0.025), upr = quantile(auc_sims, 0.975))

  t <- seq(0, 2 * pi, length.out = n_grid)
  plot_data <- list(
    conf_df = tibble::tibble(
      sp = plogis(mu_sp + s_sp * c_const * cos(t + acos(r_fixed))),
      se = plogis(mu_se + s_se * c_const * cos(t))
    ),
    pred_df = tibble::tibble(
      sp = plogis(mu_sp + s_sp_pred * c_const * cos(t + acos(r_pred))),
      se = plogis(mu_se + s_se_pred * c_const * cos(t))
    ),
    sroc_df = tibble::tibble(
      sp = sp_mc,
      se = plogis(mu_se + (sigma_se / sigma_sp) * sign(rho_re) * (logit_sp_mc - mu_sp))
    ),
    study_df = data_std |>
      dplyr::transmute(
        Study,
        study_id = seq_along(Study),
        se = TP / (TP + FN),
        sp = TN / (TN + FP)
      )
  )

  list(metrics = res, plot_data = plot_data, model = glmm_fit)
}

plot_metandi <- function(model_obj, digits = 2,
                         show_conf = TRUE,       # show the confidence contour
                         show_pred = TRUE,       # show the prediction contour
                         show_legend = TRUE,     # show the plot legend
                         custom_se = NULL,       # override summary sensitivity text
                         custom_sp = NULL,       # override summary specificity text
                         custom_auc = NULL, legend_position = "bottomright", legend_text_size = 3.5,
                         legend_bg = "#F8F9FA", font_family = "sans", base_size = 14,
                         title = "SROC with Prediction & Confidence Contours",
                         show_study_labels = TRUE, study_size = 4, study_label_size = 2.5,
                         summary_size = 4.5, summary_col = "#C0392B",
                         sroc_col = "#2C3E50", sroc_linewidth = 1.2,
                         confidence_col = "#2980B9", confidence_alpha = .2,
                         prediction_col = "#BDC3C7", prediction_alpha = .15,
                         auc_digits = digits, x_breaks = seq(0, 1, .2),
                         y_breaks = seq(0, 1, .2)) {
  for (name in c("show_conf", "show_pred", "show_legend", "show_study_labels")) {
    value <- get(name)
    if (!is.logical(value) || length(value) != 1L || is.na(value))
      stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  for (name in c("legend_text_size", "base_size", "study_size", "study_label_size",
                 "summary_size", "sroc_linewidth")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0)
      stop(name, " must be a positive finite number.", call. = FALSE)
  }
  for (name in c("digits", "auc_digits")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 0 || value > 10 || value != floor(value))
      stop(name, " must be an integer from 0 to 10.", call. = FALSE)
  }
  for (name in c("confidence_alpha", "prediction_alpha")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 0 || value > 1)
      stop(name, " must be between 0 and 1.", call. = FALSE)
  }
  for (name in c("legend_bg", "summary_col", "sroc_col", "confidence_col", "prediction_col")) {
    value <- get(name)
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        inherits(try(grDevices::col2rgb(value), silent = TRUE), "try-error"))
      stop(name, " must be one valid R color.", call. = FALSE)
  }
  if (!is.character(font_family) || length(font_family) != 1L || is.na(font_family) || !nzchar(trimws(font_family)))
    stop("font_family must be one non-empty string.", call. = FALSE)
  for (name in c("title", "custom_se", "custom_sp", "custom_auc")) {
    value <- get(name)
    if (!is.null(value) && (!is.character(value) || length(value) != 1L || is.na(value)))
      stop(name, " must be NULL or one string.", call. = FALSE)
  }
  for (name in c("x_breaks", "y_breaks")) {
    value <- get(name)
    if (!is.numeric(value) || !length(value) || any(!is.finite(value)) || any(value < 0 | value > 1) || anyDuplicated(value))
      stop(name, " must contain distinct finite values between 0 and 1.", call. = FALSE)
  }
  corners <- c("bottomright", "bottomleft", "topright", "topleft")
  if (!(is.character(legend_position) && length(legend_position) == 1L && !is.na(legend_position) && legend_position %in% corners) &&
      !(is.numeric(legend_position) && length(legend_position) == 2L && all(is.finite(legend_position)) && all(legend_position >= 0 & legend_position <= 1)))
    stop("legend_position must be a corner name or c(x, y) in screen coordinates from 0 to 1.", call. = FALSE)
  
  # 1. Extract model outputs
  pd <- model_obj$plot_data
  mt <- model_obj$metrics
  
  # 2. Format summary text
  fmt <- function(x) sprintf(paste0("%.", digits, "f"), x)
  
  str_se  <- if(!is.null(custom_se)) custom_se else sprintf("Sens: %s (%s-%s)", fmt(mt$se[1]), fmt(mt$se[2]), fmt(mt$se[3]))
  str_sp  <- if(!is.null(custom_sp)) custom_sp else sprintf("Spec: %s (%s-%s)", fmt(mt$sp[1]), fmt(mt$sp[2]), fmt(mt$sp[3]))
  fmt_auc <- function(x) sprintf(paste0("%.", auc_digits, "f"), x)
  str_auc <- if (!is.null(custom_auc)) custom_auc else if (all(is.finite(mt$auc[2:3]))) {
    sprintf("AUC: %s (%s-%s)", fmt_auc(mt$auc[1]), fmt_auc(mt$auc[2]), fmt_auc(mt$auc[3]))
  } else sprintf("AUC: %s", fmt_auc(mt$auc[1]))
  
  # 3. Define color palette
  col_sroc <- sroc_col   # SROC curve
  col_conf <- confidence_col   # confidence contour
  col_pred <- prediction_col   # prediction contour
  col_obs  <- "#7F8C8D"   # observed studies
  col_sum  <- summary_col   # summary point
  bg_box   <- legend_bg   # legend background
  
  # 4. Build the SROC plot
  p <- ggplot2::ggplot()
  
  # 4.1 Prediction contour
  if (show_pred) {
    p <- p + 
      ggplot2::geom_polygon(data = pd$pred_df, ggplot2::aes(x = sp, y = se), fill = col_pred, alpha = prediction_alpha) +
      ggplot2::geom_path(data = pd$pred_df, ggplot2::aes(x = sp, y = se), color = col_pred, linetype = "dotted", linewidth = 0.8)
  }
  
  # 4.2 Confidence contour
  if (show_conf) {
    p <- p + 
      ggplot2::geom_polygon(data = pd$conf_df, ggplot2::aes(x = sp, y = se), fill = col_conf, alpha = confidence_alpha) +
      ggplot2::geom_path(data = pd$conf_df, ggplot2::aes(x = sp, y = se), color = col_conf, linetype = "dashed", linewidth = 0.8)
  }
  
  # 4.3 SROC curve, observed studies, and summary point
  p <- p +
    ggplot2::geom_line(data = pd$sroc_df, ggplot2::aes(x = sp, y = se), color = col_sroc, linewidth = sroc_linewidth) +
    
    ggplot2::geom_point(data = pd$study_df, ggplot2::aes(x = sp, y = se), 
               shape = 21, size = study_size, fill = ggplot2::alpha("white", 0.6), color = col_obs, stroke = 1.2) +
    (if (show_study_labels) ggplot2::geom_text(data = pd$study_df, ggplot2::aes(x = sp, y = se, label = study_id),
              size = study_label_size, family = font_family, color = col_sroc) else NULL) +
    
    ggplot2::geom_point(ggplot2::aes(x = mt$sp[1], y = mt$se[1]), 
               shape = 15, size = summary_size, color = col_sum) +
    
    # Reverse the x-axis so higher specificity is on the left
    ggplot2::scale_x_reverse(limits = c(1, 0), expand = c(0, 0), breaks = x_breaks) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = y_breaks) +
    ggplot2::coord_fixed(ratio = 1) +
    
    # Add axis labels and title
    ggplot2::labs(x = "Specificity", y = "Sensitivity", 
         title = title) +
    ggplot2::theme_classic(base_size = base_size, base_family = font_family) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold", margin = ggplot2::margin(t = 10, r = 10)),
      axis.text  = ggplot2::element_text(color = "black"),
      axis.line  = ggplot2::element_line(linewidth = 0.8),
      axis.ticks = ggplot2::element_line(linewidth = 0.8),
      axis.ticks.length = grid::unit(0.2, "cm"),
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.margin = ggplot2::margin(15, 15, 15, 15)
    )
  
  # 5. Add legend annotations
  if (show_legend) {
    x0 <- 0.44       # left edge of the legend box
    x_text <- 0.38   # label anchor
    x_line <- 0.42   # line or point anchor
    # Stack legend entries from bottom to top
    y_current <- 0.07
    y_pos <- list()
    
    if (show_pred) { y_pos$pred <- y_current; y_current <- y_current + 0.07 }
    if (show_conf) { y_pos$conf <- y_current; y_current <- y_current + 0.09 }
    
    y_pos$sroc <- y_current; y_current <- y_current + 0.11
    y_pos$sum <- y_current; y_current <- y_current + 0.09
    y_pos$obs <- y_current
    
    # Compute legend box bounds
    box_ymax <- y_pos$obs + 0.04
    box_ymin <- 0.03
    scale <- legend_text_size / 3.5
    box_height <- (box_ymax - box_ymin) * scale
    anchor <- if (is.character(legend_position)) {
      c(if (grepl("right", legend_position)) .54 else .02,
        if (grepl("top", legend_position)) .98 - box_height else .03)
    } else legend_position
    if (anchor[1] + .44 > 1 || anchor[2] + box_height > 1 || any(anchor < 0))
      stop("Legend does not fit: reduce legend_text_size or move legend_position.", call. = FALSE)
    dx <- .54 - anchor[1]
    x0 <- x0 + dx
    x_text <- x_text + dx
    x_line <- x_line + dx
    y_pos <- lapply(y_pos, function(y) anchor[2] + (y - .03) * scale)
    box_ymin <- anchor[2]
    box_ymax <- anchor[2] + box_height
    
    p <- p +
      ggplot2::annotate("rect", xmin = x0+0.02, xmax = 0.02 + dx, ymin = box_ymin, ymax = box_ymax,
               fill = bg_box, color = "gray80", linewidth = 0.5) +
      
      # Observed study marker
      ggplot2::annotate("point", x = x_line, y = y_pos$obs, shape = 21, size = study_size, color = col_obs, fill = "white", stroke = 1.2) +
      ggplot2::annotate("text", x = x_text, y = y_pos$obs, label = "Observed Data", hjust = 0, size = legend_text_size * 4 / 3.5, family = font_family, fontface = "bold") +
      
      ggplot2::annotate("point", x = x_line, y = y_pos$sum, shape = 15, size = summary_size, color = col_sum) +
      ggplot2::annotate("text", x = x_text, y = y_pos$sum, 
               label = sprintf("Summary Point\n%s\n%s", str_se, str_sp), 
               hjust = 0, size = legend_text_size, family = font_family, lineheight = 1.1) +
      
      ggplot2::annotate("segment", x = x_line+0.02, xend = x_line-0.02, y = y_pos$sroc, yend = y_pos$sroc, 
               color = col_sroc, linewidth = sroc_linewidth) +
      ggplot2::annotate("text", x = x_text, y = y_pos$sroc, 
               label = sprintf("SROC Curve\n%s", str_auc), 
               hjust = 0, size = legend_text_size, family = font_family, lineheight = 1.1)
    
    # Add contour labels if requested
    if (show_conf) {
      p <- p + 
        ggplot2::annotate("segment", x = x_line+0.02, xend = x_line-0.02, y = y_pos$conf, yend = y_pos$conf, 
                 color = col_conf, linewidth = 0.8, linetype = "dashed") +
        ggplot2::annotate("text", x = x_text, y = y_pos$conf, label = "95% Confidence Contour", hjust = 0, size = legend_text_size, family = font_family)
    }
    
    if (show_pred) {
      p <- p + 
        ggplot2::annotate("segment", x = x_line+0.02, xend = x_line-0.02, y = y_pos$pred, yend = y_pos$pred, 
                 color = col_pred, linewidth = 0.8, linetype = "dotted") +
        ggplot2::annotate("text", x = x_text, y = y_pos$pred, label = "95% Prediction Contour", hjust = 0, size = legend_text_size, family = font_family)
    }
  }
  
  # 6. Print summary results to the console
  cat(sprintf("\n=== Bivariate Meta-Analysis Results ===\n"))
  cat(sprintf("%s\n%s\n%s\n", str_se, str_sp, str_auc))
  cat(sprintf("============================================\n"))
  
  return(p)
}
