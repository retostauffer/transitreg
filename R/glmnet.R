

transitreg_glmnet <- function(formula, data, nfolds = 10, ...) {
    stopifnot(requireNamespace("glmnet"))

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

predict.transietreg_glmnet <- function(object, type = "response", s = "lambda.min", ...) {
    # Removing custom class and call predict.glmnet below
    class(object) <- class(object)[!class(object) == "transitreg_glmnet"]

    return(predict(object, newx = attr(object, "X"), s = s, type = type))
}


