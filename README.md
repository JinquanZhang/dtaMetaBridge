# dtaMetaBridge

`dtaMetaBridge` connects paired `meta::metaprop()` analyses with standard
diagnostic-test-accuracy meta-analysis.

It provides paired forest plots, a Reitsma bivariate random-effects model,
Rutter-Gatsonis HSROC curves, and prevalence-dependent post-test probabilities.

## Installation

```r
remotes::install_github("YOUR_GITHUB_USER/dtaMetaBridge")
```

## Typical workflow

```r
library(meta)
library(dtaMetaBridge)

meta_sens <- metaprop(TP, TP + FN, studlab = study, data = dta,
                      method = "GLMM", method.tau = "ML")
meta_spec <- metaprop(TN, TN + FP, studlab = study, data = dta,
                      method = "GLMM", method.tau = "ML")

# Forest-plot diamond uses the random-effects estimates retained by meta.
plot_sensspec_forest_meta(meta_sens, meta_spec, output_file = "forest.png")

# Reitsma / Rutter-Gatsonis HSROC model.
fit <- fit_bivariate_meta(meta_sens, meta_spec)
plot_sroc(fit)
posttest_probability(fit, prevalence = c(0.1, 0.3, 0.5))
```

## Interpretation

Forest-plot diamonds are independent univariate random-effects summaries. The
SROC summary point is from a joint Reitsma model, so they can differ slightly.
The HSROC AUC is a cross-study summary measure and is not interchangeable with
the ROC AUC of an individual study's continuous score.
