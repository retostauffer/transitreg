
## Version 0.1-9011-90XX

* Further (major) changes as the package is still under development.
* Revamped how 'discrete' and 'pseudo-discrete' (continuous) models
  are handled, as well as how to depict this in `Transition` disbributions.
  For short, if `breaks` are integers, the model/distribution is assumed to 
  be a discrete count data distribution. Else a pseudo-discrete continuous
  distribution.
* Convinced `topmodels` to properly handle both (discrete and pseudo-discrete)
  distributions using a custom `proresiduals` S3 method.
* Refined plotting (to be in line with the `distributions3` plotting functions).
* Updated man pages where needed (tbc).
* `predict()` now supports `type = "mean"` to calculate the expectation (weighted mean).
* The `transitreg()` formula interface should not properly handle 
  transformed conversion objects (i.e., `sqrt(y)`, `I(y^2)`, ...).
* Implemented GitHub action for running tinytests on push; tests need to be
  extended to cover all features/use cases.
* Re-implemented censoring (uncensored, left-censored, right-censored, or both).
  Censoring points are pseudo-bins with a width of 0, handled as point mass probabilities.
* Refractured how breaks are calculated for different scenarios.


## Version 0.1-9010

* Migrated back to main branch, repository renamed to final
  package name `transitreg`.
* Major internal changes, refined several functions and methods,
  updated examples and documentation (tbc).
* Renamed `pmax` to `mode`.

## Version 0.1-9004

* Major changes to the C code (removing the
  dependency to know the original (numeric) response), fixed a few
  smaller bugs.
* Added proper dpqr support in addition to distributions3 methods:
  dtransit, ptransit, qtransit, and rtransit.
* Added proper handling for dpq functions if outside support.
* Migrated to roxygen2 documentation.
* Updated the prediction methods and functions.
* Added (full) topmodels support.

## Version 0.1-9003

* Renamed the package from it's working title `TransitionModels`
  to `transitreg` (Transitional Model Regression).
* Various smaller and larger changes to both increase functionality
  and performance (speed and memory footprint).
* Added distributions3 support for the new model class via the
  `Transition()` distribution object.

## Version 0.1-9002

* Removed original base-R implementation which have been replaced
    by vectorised versions and C functions.
* Added an experimental \code{engine = "glmnet"}.

## Version 0.1-9001

* Added proper `NA` handling to `tm_calc_cdf`, `tm_calc_pdf`,
    and `tm_calc_pmax` in C.
* Made the default development mode less verbose.
* Fixed R-version of `predict(...)` to properly deal with missing values
    when `type = "pmax"` or `type = "quantile"` (may be obsolete in the near future).
* Implemented C-function for quantile prediction/evaluation.

## Version 0.1-9000

* Evaluation/calculation of PDF, CDF, and 'pmax' implemented
  in C for higher performance.
* Added S3 method to extract model coefficients.
* Several smaller speed improvements and code/coding style touch-ups.
