# Rscript tools/validate-hsroc.R path/to/HFphf_meta_analysis.qmd
devtools::load_all(".", quiet = TRUE)
suppressPackageStartupMessages(library(dplyr))
qmd <- paste(readLines(commandArgs(TRUE)[1], encoding = "UTF-8", warn = FALSE), collapse = "\n")
for (name in c("df_h2in", "df_h2out", "df_hfain", "df_hfaout")) {
  block <- regmatches(qmd, regexpr(paste0("(?s)", name, " <- data.frame\\(.*?\\n\\)"), qmd, perl = TRUE))
  env <- new.env()
  eval(parse(text = block)[[1]], envir = env)
  d <- env[[name]]
  fit <- fit_bivariate_dta(d)
  ref <- mada::reitsma(d, method = "reml", correction = .5, correction.control = "single")
  fpr <- 1 - fit$plot_data$sroc$sp
  reference_curve <- mada::sroc(ref, fpr = fpr, type = "ruttergatsonis")
  stopifnot(isTRUE(all.equal(fit$plot_data$sroc$se, unname(reference_curve[,2]), tolerance = 1e-10)),
    all(diff(fit$plot_data$sroc$se) >= 0))
  mu <- as.numeric(stats::coef(ref))
  slope <- sqrt(ref$Psi[1,1] / ref$Psi[2,2])
  x <- seq(0, 1, length.out = 100001)
  y <- plogis(mu[1] + slope * (qlogis(x) - mu[2]))
  area <- sum(diff(x) * (head(y,-1) + tail(y,-1)) / 2)
  stopifnot(abs(area - fit$metrics$auc[["est"]]) < 1e-5)
  cat(sprintf("%s: Se=%.6f Sp=%.6f AUC=%.6f; curve and integration checks passed\n", name, fit$metrics$sensitivity[1], fit$metrics$specificity[1], fit$metrics$auc[1]))
}
