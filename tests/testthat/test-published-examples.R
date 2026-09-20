# Published-data checks for the diagnostic-meta helpers.
# Sources accessed 2026-09-20:
# - Tada et al. 2021, PMID 34364907: H2FPEF low-score rule-out sensitivity 97%
#   in 194 HFpEF and 178 non-HFpEF participants.
# - Forsyth et al. 2021, doi:10.1002/ehf2.13612, Table 3: HFA-PEFF >5 was
#   present in 58.2% (39/67) of HFpEF and 44.4% (16/36) of non-HFpEF patients.

# Tada et al.: score 0-1 is the rule-out category.  The QMD row's 188/194
# sensitivity reproduces the article's rounded 97% result.
testthat::test_that("published examples reproduce reported values", {
  tada_h2_ruleout <- data.frame(study = "Tada 2021", TP = 188, FP = 100, FN = 6, TN = 78)
  tada_result <- dtaMetaBridge:::.forest_data(tada_h2_ruleout)$studies
  testthat::expect_equal(round(tada_result$sens, 2), 0.97)
  testthat::expect_equal(tada_result$TP + tada_result$FN, 194)
  testthat::expect_equal(tada_result$FP + tada_result$TN, 178)

# Forsyth et al.: HFA-PEFF >5 is a rule-in result. These counts reproduce both
# article percentages without any model fitting.
  forsyth_hfa_rulein <- data.frame(study = "Forsyth 2021", TP = 39, FP = 16, FN = 28, TN = 20)
  forsyth_result <- dtaMetaBridge:::.forest_data(forsyth_hfa_rulein)$studies
  testthat::expect_equal(round(forsyth_result$sens, 3), 0.582)
  testthat::expect_equal(round(forsyth_result$spec, 3), 0.556)

# Reddy et al. 2022, PMCID PMC9280610, reports 2 false negatives among score-0
# and 15 among score-1 patients. This independently confirms the QMD rule-out
# row's 17 false-negative events, but not its complete 2x2 table.
  reddy_h2_ruleout <- data.frame(study = "Reddy 2022", TP = 366, FP = 50, FN = 17, TN = 52)
  testthat::expect_equal(reddy_h2_ruleout$FN, 2 + 15)
})
