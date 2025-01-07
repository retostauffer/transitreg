
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
