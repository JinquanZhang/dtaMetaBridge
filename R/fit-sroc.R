# Formula numbering follows meta4diag::SROC.
.sroc_slope <- function(Psi, type) {
  a <- Psi[1, 1]; b <- Psi[2, 2]; c <- Psi[1, 2]
  if (any(!is.finite(Psi)) || a <= 0 || b <= 0)
    stop("SROC requires positive finite between-study variances.", call. = FALSE)
  tol <- 1e-10 * sqrt(a * b)
  if (type %in% c(2L, 4L) && abs(c) <= tol)
    stop("This SROC formula is undefined or unstable at zero covariance; use type 1 or 5.", call. = FALSE)
  if (type == 3L && abs(b + c) <= 1e-10 * max(a, b))
    stop("Type 3 has a zero or near-zero denominator.", call. = FALSE)
  # ponytail: the major axis is the eigenvector of the largest eigenvalue.
  slope <- switch(type, c / b,
    if (a >= b) (a - b + sqrt((a - b)^2 + 4 * c^2)) / (2 * c)
    else 2 * c / (b - a + sqrt((a - b)^2 + 4 * c^2)),
    (a + c) / (b + c), a / c, sqrt(a / b))
  if (!is.finite(slope)) stop("Non-finite SROC slope.", call. = FALSE)
  slope
}

.sroc_values <- function(fpr, mu, slope) {
  if (slope == 0) return(rep(stats::plogis(mu[1]), length(fpr)))
  stats::plogis(mu[1] + slope * (stats::qlogis(fpr) - mu[2]))
}

#' Fit a bivariate model and select one of five SROC formulas.
.dtameta_fit <- function(data, study_col, conf_level) {
  if (!requireNamespace("dtametaTMB", quietly = TRUE))
    stop("backend = 'dtametaTMB' requires the dtametaTMB package.", call. = FALSE)
  data$study <- data[[study_col]]
  dtametaTMB::fitReitsma(data = data, TP = TP, FP = FP, FN = FN, TN = TN,
    study = study, conflevel = conf_level)
}

.dtameta_parameters <- function(model) {
  e <- model$estimates
  mu <- c(unname(e["mu_A.sens", "Estimate"]), -unname(e["mu_B.spec", "Estimate"]))
  psi <- matrix(c(e["sigma2_A.sens", "Estimate"], -e["sigma_AB", "Estimate"],
                  -e["sigma_AB", "Estimate"], e["sigma2_B.spec", "Estimate"]), 2, 2)
  fixed <- as.matrix(model$vcov)[1:2, 1:2, drop = FALSE]
  fixed[2, ] <- -fixed[2, ]; fixed[, 2] <- -fixed[, 2]
  list(mu = mu, psi = psi, fixed = fixed)
}

.fit_sroc_dtameta <- function(data, study_col, sroc_type, conf_level, n_grid,
  auc_boot, seed) {
  model <- .dtameta_fit(data, study_col, conf_level)
  par <- .dtameta_parameters(model)
  slope <- .sroc_slope(par$psi, sroc_type)
  rho <- par$psi[1, 2] / sqrt(par$psi[1, 1] * par$psi[2, 2])
  if (abs(rho) > .999)
    warning("Between-study correlation is near its boundary; SROC formulas may coincide and estimates may be unstable.", call. = FALSE)
  if (slope <= 0)
    warning("Selected formula is flat or decreasing versus FPR; AUC is omitted. Consider type 5 and inspect model suitability.", call. = FALSE)
  curve <- function(x) .sroc_values(x, par$mu, slope)
  fpr <- data$FP / (data$FP + data$TN)
  auc <- if (slope > 0) stats::integrate(curve, 0, 1, rel.tol = 1e-8)$value else NA_real_
  auc_ci <- c(lwr = NA_real_, upr = NA_real_)
  if (auc_boot > 0L && slope > 0) {
    # Do not leave a modelling helper with a changed global random-number state.
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
    boot_auc <- vapply(seq_len(auc_boot), function(i) {
      d <- data[sample.int(nrow(data), nrow(data), replace = TRUE), , drop = FALSE]
      d[[study_col]] <- make.unique(as.character(d[[study_col]]), sep = " #")
      tryCatch({
        p <- .dtameta_parameters(.dtameta_fit(d, study_col, conf_level))
        b <- .sroc_slope(p$psi, sroc_type)
        if (b <= 0) return(NA_real_)
        stats::integrate(function(x) .sroc_values(x, p$mu, b), 0, 1, rel.tol = 1e-7)$value
      }, error = function(e) NA_real_)
    }, numeric(1))
    boot_auc <- boot_auc[is.finite(boot_auc)]
    if (length(boot_auc) < max(20L, auc_boot * .75))
      warning("Too many Bootstrap refits failed; AUC interval is unavailable.", call. = FALSE)
    else auc_ci <- stats::quantile(boot_auc, c(.025, .975), names = FALSE)
  }
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  sens <- c(est = stats::plogis(par$mu[1]),
    lwr = stats::plogis(par$mu[1] - z * sqrt(par$fixed[1, 1])),
    upr = stats::plogis(par$mu[1] + z * sqrt(par$fixed[1, 1])))
  spec <- c(est = stats::plogis(-par$mu[2]),
    lwr = stats::plogis(-par$mu[2] - z * sqrt(par$fixed[2, 2])),
    upr = stats::plogis(-par$mu[2] + z * sqrt(par$fixed[2, 2])))
  ellipse_roc <- function(v) {
    e <- ellipse::ellipse(v, centre = par$mu, level = conf_level)
    data.frame(sp = 1 - stats::plogis(e[, 2]), se = stats::plogis(e[, 1]))
  }
  grid <- seq(max(.001, min(fpr)), min(.999, max(fpr)), length.out = n_grid)
  list(model = model, backend = "frequency", model_type = "dtametaTMB::fitReitsma",
    sroc_type = sroc_type, curve_parameters = list(mu = par$mu, slope = slope),
    input_data = data, interval_type = "95% CI",
    metrics = list(sensitivity = sens, specificity = spec,
      auc = c(est = auc, lwr = unname(auc_ci[1]), upr = unname(auc_ci[2])),
      pauc = if (slope > 0 && diff(range(fpr)) > 0) stats::integrate(curve, min(fpr), max(fpr))$value else NA_real_),
    plot_data = list(confidence = ellipse_roc(par$fixed),
      prediction = ellipse_roc(par$fixed + par$psi),
      sroc = data.frame(sp = 1 - grid, se = curve(grid)),
      studies = data.frame(Study = as.character(data[[study_col]]),
        specificity = data$TN / (data$TN + data$FP), sensitivity = data$TP / (data$TP + data$FN))),
    random_effects = list(covariance = par$psi))
}

.meta4diag_interval <- function(x, row, conf_level) {
  alpha <- (1 - conf_level) / 2
  qcols <- grep("quant$", colnames(x), value = TRUE)
  probs <- suppressWarnings(as.numeric(sub("quant$", "", qcols)))
  if (length(qcols) < 2L || anyNA(probs))
    stop("meta4diag did not return the requested posterior quantiles.", call. = FALSE)
  c(est = x[row, "mean"],
    lwr = x[row, qcols[which.min(abs(probs - alpha))]],
    upr = x[row, qcols[which.min(abs(probs - (1 - alpha)))]] )
}

.fit_sroc_meta4diag <- function(data, study_col, sroc_type, conf_level, n_grid,
  posterior_samples, seed) {
  if (!requireNamespace("meta4diag", quietly = TRUE) || !requireNamespace("INLA", quietly = TRUE))
    stop("backend = 'meta4diag' requires both meta4diag and INLA.", call. = FALSE)
  attached_inla <- "package:INLA" %in% search()
  if (!attached_inla) {
    base::library("INLA", character.only = TRUE)
    on.exit(detach("package:INLA", unload = FALSE, character.only = TRUE), add = TRUE)
  }
  d <- data
  d$studynames <- as.character(data[[study_col]])
  alpha <- (1 - conf_level) / 2
  model <- meta4diag::meta4diag(d, model.type = 1, link = "logit",
    quantiles = c(alpha, .5, 1 - alpha), nsample = posterior_samples,
    seed = seed, verbose = FALSE)
  sf <- model$summary.fixed
  sh <- model$summary.hyperpar
  mu <- c(sf["mu", "mean"], -sf["nu", "mean"])
  psi <- matrix(c(sh["var_phi", "mean"],
    -sh["cor", "mean"] * sqrt(sh["var_phi", "mean"] * sh["var_psi", "mean"]),
    -sh["cor", "mean"] * sqrt(sh["var_phi", "mean"] * sh["var_psi", "mean"]),
    sh["var_psi", "mean"]), 2, 2)
  slope <- .sroc_slope(psi, sroc_type)
  if (slope <= 0)
    warning("Selected formula is flat or decreasing versus FPR; AUC is omitted. Consider type 5 and inspect model suitability.", call. = FALSE)
  curve <- function(x) .sroc_values(x, mu, slope)
  auc_raw <- meta4diag::AUC(model, sroc.type = sroc_type, est.type = "mean")
  qauc <- grep("quant$", names(auc_raw), value = TRUE)
  qprob <- suppressWarnings(as.numeric(sub("quant$", "", qauc)))
  auc <- c(est = if (slope > 0) unname(auc_raw["est"]) else NA_real_,
    lwr = if (length(qauc)) unname(auc_raw[qauc[which.min(abs(qprob - alpha))]]) else NA_real_,
    upr = if (length(qauc)) unname(auc_raw[qauc[which.min(abs(qprob - (1 - alpha)))]] ) else NA_real_)
  ss <- model$summary.expected.accuracy
  fixed <- stats::cov(t(rbind(model$samples.fixed["mu", ], -model$samples.fixed["nu", ])))
  ellipse_roc <- function(v) {
    e <- ellipse::ellipse(v, centre = mu, level = conf_level)
    data.frame(sp = 1 - stats::plogis(e[, 2]), se = stats::plogis(e[, 1]))
  }
  fpr <- data$FP / (data$FP + data$TN)
  grid <- seq(max(.001, min(fpr)), min(.999, max(fpr)), length.out = n_grid)
  list(model = model, backend = "bayes", model_type = "meta4diag::meta4diag",
    sroc_type = sroc_type, curve_parameters = list(mu = mu, slope = slope),
    input_data = data, interval_type = "95% CrI",
    metrics = list(sensitivity = .meta4diag_interval(ss, "mean(Se)", conf_level),
      specificity = .meta4diag_interval(ss, "mean(Sp)", conf_level), auc = auc,
      pauc = if (slope > 0 && diff(range(fpr)) > 0) stats::integrate(curve, min(fpr), max(fpr))$value else NA_real_),
    plot_data = list(confidence = ellipse_roc(fixed), prediction = ellipse_roc(fixed + psi),
      sroc = data.frame(sp = 1 - grid, se = curve(grid)),
      studies = data.frame(Study = as.character(data[[study_col]]),
        specificity = data$TN / (data$TN + data$FP), sensitivity = data$TP / (data$TP + data$FN))),
    random_effects = list(covariance = psi))
}

fit_sroc <- function(data, backend = c("frequency", "bayes", "legacy"),
  sroc_type = 5L, study_col = "study", year_col = "Year",
  conf_level = .95, n_grid = 1000, auc_boot = 0L, posterior_samples = 2000L, seed = 2026) {
  if (!is.data.frame(data)) stop("data must be a data.frame.", call. = FALSE)
  if (!is.numeric(sroc_type) || length(sroc_type) != 1L || is.na(sroc_type) ||
      !sroc_type %in% 1:5) stop("sroc_type must be an integer from 1 to 5.", call. = FALSE)
  backend <- match.arg(backend)
  if (!is.numeric(n_grid) || length(n_grid) != 1L || !is.finite(n_grid) ||
      n_grid < 100 || n_grid != floor(n_grid)) stop("n_grid must be an integer >= 100.", call. = FALSE)
  if (!is.numeric(conf_level) || length(conf_level) != 1L || conf_level <= 0 || conf_level >= 1)
    stop("conf_level must be strictly between 0 and 1.", call. = FALSE)
  if (!is.numeric(auc_boot) || length(auc_boot) != 1L || auc_boot < 0 || auc_boot != floor(auc_boot))
    stop("auc_boot must be a non-negative integer.", call. = FALSE)
  if (!is.numeric(posterior_samples) || length(posterior_samples) != 1L || posterior_samples < 20 || posterior_samples != floor(posterior_samples))
    stop("posterior_samples must be an integer >= 20.", call. = FALSE)
  if (!is.character(study_col) || length(study_col) != 1L || is.na(study_col))
    stop("study_col must be one column name.", call. = FALSE)
  if (!is.null(year_col) && (!is.character(year_col) || length(year_col) != 1L || is.na(year_col)))
    stop("year_col must be NULL or one column name.", call. = FALSE)
  .assert_dta_data(data, study_col)
  labels <- as.character(data[[study_col]])
  if (anyNA(labels) || any(!nzchar(trimws(labels)))) stop("Study labels must not be missing or empty.", call. = FALSE)
  if (!is.null(year_col) && year_col %in% names(data)) labels <- paste(labels, data[[year_col]])
  if (anyDuplicated(labels)) warning("Duplicate study/year labels: rows retained separately. Verify independent cohorts; repeated thresholds are not independent studies.", call. = FALSE)
  d <- data
  d[[study_col]] <- make.unique(labels, sep = " #")
  if (backend == "frequency") {
    out <- .fit_sroc_dtameta(d, study_col, as.integer(sroc_type), conf_level,
      n_grid, as.integer(auc_boot), seed)
    out$input_data <- data
    return(out)
  }
  if (backend == "bayes") {
    out <- .fit_sroc_meta4diag(d, study_col, as.integer(sroc_type), conf_level,
      n_grid, as.integer(posterior_samples), seed)
    out$input_data <- data
    return(out)
  }
  fit <- fit_bivariate_dta(d, study_col = study_col, n_grid = n_grid)
  Psi <- fit$model$Psi
  mu <- unname(stats::coef(fit$model)["(Intercept)", ])
  slope <- .sroc_slope(Psi, as.integer(sroc_type))
  rho <- Psi[1, 2] / sqrt(Psi[1, 1] * Psi[2, 2])
  if (abs(rho) > .999) warning("Between-study correlation is near its boundary; SROC formulas may coincide and estimates may be unstable.", call. = FALSE)
  if (slope <= 0) warning("Selected formula is flat or decreasing versus FPR; AUC is omitted. Consider type 5 and inspect model suitability.", call. = FALSE)
  fit$sroc_type <- as.integer(sroc_type)
  fit$curve_parameters <- list(mu = mu, slope = slope)
  fit$model_type <- paste("mada::reitsma SROC formula", sroc_type)
  fit$input_data <- data
  fit$plot_data$sroc$se <- .sroc_values(1 - fit$plot_data$sroc$sp, mu, slope)
  curve <- function(x) .sroc_values(x, mu, slope)
  observed <- range(data$FP / (data$FP + data$TN))
  fit$metrics$auc <- c(est = if (slope > 0) stats::integrate(curve, 0, 1, rel.tol = 1e-8)$value else NA_real_, lwr = NA_real_, upr = NA_real_)
  fit$metrics$pauc <- if (slope <= 0) NA_real_ else if (diff(observed) == 0) 0 else
    stats::integrate(curve, observed[1], observed[2], rel.tol = 1e-8)$value
  fit
}
