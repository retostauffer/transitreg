
## Version 0.1-9004

* Still in development, major changes to the C code (removing the
  dependency to know the original (numeric) response), fixed a few
  smaller bugs.
* Migrated to roxygen2 documentation.
* Updated the prediction methods and functions.

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
