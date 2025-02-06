

#' Estimat Transition Probabilities via `glmnet`
#'
#' Experimental engine used to estimate the transition probabilities
#' in a [transitreg()] model.
#'
#' @param formula An object of class `formula`.
#' @param data A data.frame containing the required data to fit a binomial
#'        regression model (given `formula`) using [glmnet::cv.glmnet()].
#' @param nfolds Integer, defaults to `10L`. Number of cross-folds for the
#'        glmnet model.
#' @param \dots forwarded to [glmnet::cv.glmnet()].
#'
#' @importFrom mgcv interpret.gam smoothCon
#'
#' @rdname transitreg_glmnet
#' @export
transitreg_glmnet <- function(formula, data, nfolds = 10, ...) {
    ## TODO(R): This is very experimental! Known issues:
    ## - Nfolds = 10 splits the data into 10 folds, ignoring that
    ##   each observation shows up 1-N times in 'data'.
    ## - Users have no option to control glmnet.
    ## - This is really just a test!
    warning("TODO(R): Experimental implementation!")
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package 'glmnet' is required but not installed.")
    }

    # Extract smooth terms
    smooth_terms <- interpret.gam(as.formula(formula))$smooth

    # Build basis functions using smoothCon()
    fn1 <- function(term, x) {
        smoothCon(eval(term), data = x, knots = NULL)[[1]]
    }
    basis_list <- lapply(smooth_terms, fn1, x = data)

    # Combine basis functions into a design matrix
    fn2 <- function(basis) {
        X <- basis$X
        colnames(X) <- sprintf("%s_%d", basis$label, seq_len(ncol(X)))
        return(X)
    }
    X <- do.call(cbind, lapply(basis_list, fn2)) # Design matrix

    # Convert to matrix for glmnet
    y <- data[, response_name(formula)]

    # Fit the model using glmnet
    mod <- glmnet::cv.glmnet(X, y, nfolds = nfolds, family = "binomial", ...)

    # We need to store the model matrix and the response to be used
    # when calling predict later on.
    attr(mod, "X") <- X
    attr(mod, "y") <- y

    class(mod) <- c("transitreg_glmnet", class(mod))
    print(class(mod))
    return(mod)
}

#' @param object an object of class `transitreg_glmnet`.
#' @param type character, defaults to `"response"`.
#' @param s value(s) of the penalty parameter, defaults to `"lambda.min"`.
#'        See [glmnet::predict.glmnet()] for details.
#' @param \dots currently unused.
#'
#' @exportS3Method predict transitreg_glmnet
#' @rdname transitreg_glmnet
predict.transitreg_glmnet <- function(object, type = "response", s = "lambda.min", ...) {
    # Removing custom class and call predict.glmnet below
    class(object) <- class(object)[!class(object) == "transitreg_glmnet"]

    return(predict(object, newx = attr(object, "X"), s = s, type = type))
}


