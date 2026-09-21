# Formula reference: local Stata midas.ado, SUMMARY ROC CURVE block.
# The fitting model is binomial-logit with correlated study random effects.
.fit_midas <- function(data, study_col, n_grid) {
  if (length(n_grid) != 1L || !is.finite(n_grid) || n_grid < 100 || n_grid != floor(n_grid)) stop("n_grid must be an integer >= 100.", call. = FALSE)
  k <- nrow(data)
  long <- data.frame(
    study = factor(rep(seq_len(k), 2)),
    endpoint = factor(rep(c("sens", "spec"), each = k), levels = c("sens", "spec")),
    event = c(data$TP, data$TN), failure = c(data$FN, data$FP)
  )
  model <- lme4::glmer(cbind(event, failure) ~ 0 + endpoint + (0 + endpoint | study),
    data = long, family = stats::binomial(), nAGQ = 1,
    control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000)))
  mu <- unname(lme4::fixef(model))
  V <- as.matrix(stats::vcov(model))
  Psi <- as.matrix(lme4::VarCorr(model)$study)
  if (any(!is.finite(c(mu, V, Psi)))) stop("Non-finite bivariate model estimates.", call. = FALSE)
  # Midas floors both random-effect variances at 0.001 in the curve formula.
  b <- (max(.001, Psi[2, 2]) / max(.001, Psi[1, 1]))^.25
  a <- mu[1] * b + mu[2] / b
  curve <- function(sp) stats::plogis((a - stats::qlogis(sp) / b) / b)
  sp <- seq(0, 1, length.out = n_grid)
  # Midas integrates over 500 points, independently of the plotting grid.
  auc_sp <- seq(0, 1, length.out = 500)
  auc_se <- curve(auc_sp)
  auc <- sum(diff(auc_sp) * (head(auc_se, -1) + tail(auc_se, -1)) / 2)
  observed_sp <- data$TN / (data$TN + data$FP)
  limits <- range(observed_sp)
  pauc <- if (diff(limits) == 0) 0 else stats::integrate(curve, limits[1], limits[2])$value
  # Midas uses sqrt(2 * F[2,k-2]) for confidence and prediction contours.
  contour <- function(covariance) {
    xy <- ellipse::ellipse(covariance, centre = mu, t = sqrt(2 * stats::qf(.95, 2, k - 2)), npoints = 500)
    data.frame(sp = stats::plogis(xy[, 2]), se = stats::plogis(xy[, 1]))
  }
  summary <- function(i) setNames(stats::plogis(mu[i] + c(0, -1, 1) * stats::qnorm(.975) * sqrt(V[i, i])), c("est", "lwr", "upr"))
  list(model = model, backend = "midas", model_type = "Midas-style binomial GLMM (Laplace)",
    metrics = list(sensitivity = summary(1), specificity = summary(2),
      auc = c(est = auc, lwr = NA_real_, upr = NA_real_), pauc = pauc),
    plot_data = list(confidence = contour(V), prediction = contour(V + Psi),
      sroc = data.frame(sp = sp, se = curve(sp)),
      studies = data.frame(Study = as.character(data[[study_col]]), specificity = observed_sp, sensitivity = data$TP / (data$TP + data$FN))),
    random_effects = list(covariance = Psi, alpha = a, beta = b),
    diagnostics = list(singular = lme4::isSingular(model), convergence = model@optinfo$conv$lme4$messages))
}
