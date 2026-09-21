# Usage: Rscript tools/validate-qmd.R path/to/HFphf_meta_analysis.qmd
devtools::load_all(".", quiet = TRUE)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
qmd <- paste(readLines(commandArgs(TRUE)[1], encoding = "UTF-8", warn = FALSE), collapse = "\n")
start <- regexpr("fit_metandi <- function", qmd, fixed = TRUE)[1]
end <- regexpr("calc_ppv_npv_from_metandi <- function", qmd, fixed = TRUE)[1]
original <- new.env(parent = globalenv())
# Only the two selected function definitions are evaluated, never the document.
eval(parse(text = substr(qmd, start, end - 1)), envir = original)
for (name in c("df_h2in", "df_h2out", "df_hfain", "df_hfaout")) {
  block <- regmatches(qmd, regexpr(paste0("(?s)", name, " <- data.frame\\(.*?\\n\\)"), qmd, perl = TRUE))
  env <- new.env()
  eval(parse(text = block)[[1]], envir = env)
  d <- env[[name]]
  ref <- original$fit_metandi(d)
  got <- fit_bivariate_dta(d)
  stopifnot(isTRUE(all.equal(got$metrics[c("se", "sp", "auc")], ref$metrics, tolerance = 1e-10)),
    isTRUE(all.equal(got$plot_data, ref$plot_data, tolerance = 1e-10)))
  p <- plot_sroc(got, show_legend = TRUE)
  stopifnot(inherits(p, "ggplot"))
  cat(name, "matches original QMD; AUC =", got$metrics$auc[1], "\n")
  if (name == "df_h2in") plot_sroc(got, output_file = file.path(tempdir(), "qmd-h2in-sroc.png"), width = 6, height = 6)
}
cat("Preview:", file.path(tempdir(), "qmd-h2in-sroc.png"), "\n")
