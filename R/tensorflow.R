transitreg_tensorflow <- function(formula, data, nlayers = 2, units = 20,
                                  epochs = 1000, batch_size = 16, activation = "relu", dropout = 0.1,
                                  validation_split = 0.2, patience = 10, trace = 0,
                                  ncores = 1, ...)
{
  ## batchsize: 16-64, 32-128, 128-512+

  if (!requireNamespace("keras3", quietly = TRUE)) {
    stop("Package 'keras3' is required but not installed.")
  }
  if (!requireNamespace("tensorflow", quietly = TRUE)) {
    stop("Package 'tensorflow' is required but not installed.")
  }

  # Set number of threads (if tensorflow is available)
  if (!is.null(ncores)) {
    available_cores <- parallel::detectCores(logical = TRUE)
    ncores <- as.integer(max(1L, min(ncores, available_cores - 1L)))
    tf <- tensorflow::tf
    tf$config$threading$set_intra_op_parallelism_threads(ncores)
    tf$config$threading$set_inter_op_parallelism_threads(ncores)
  }

  # Response
  y <- as.integer(data[[response_name(formula)]])
  y <- tensorflow::tf$convert_to_tensor(y, dtype = "float32")

  # Regressors
  X <- model.matrix(formula, data = data)
  X <- tensorflow::tf$convert_to_tensor(X, dtype = "float32")

  # Ensure units vector has correct length
  units <- as.integer(rep(units, length.out = nlayers))

  # Create model
  mod <- keras3::keras_model_sequential()

  # Add first layer with input_shape
  mod <- mod |>
    keras3::layer_dense(units = units[1], activation = activation, input_shape = c(as.integer(ncol(X))))
  if (dropout > 0) {
    mod <- mod |> keras3::layer_dropout(rate = dropout)
  }

  # Add remaining hidden layers (if any)
  if (nlayers > 1) {
    for (i in 2:nlayers) {
      mod <- mod |>
        keras3::layer_dense(units = units[i], activation = activation)
      if (dropout > 0) {
        mod <- mod |> keras3::layer_dropout(rate = dropout)
      }
    }
  }

  # Output layer
  mod <- mod |>
    keras3::layer_dense(units = 1, activation = "sigmoid")

  # Compile model
  mod <- mod |>
    keras3::compile(
      optimizer = "adam",
      loss = "binary_crossentropy",
      metrics = "accuracy"
    )

  # Callback
  callbacks <- list(
    keras3::callback_early_stopping(patience = patience)
  )

  # Fit model
  history <- mod |>
    keras3::fit(
      x = X, y = y,
      epochs = epochs,
      batch_size = batch_size,
      validation_split = validation_split,
      verbose = trace,
      callbacks = callbacks
    )

  # Store formula for prediction
  attr(mod, "formula") <- formula
  class(mod) <- c("transitreg_tensorflow", class(mod))

  return(mod)
}

predict.transitreg_tensorflow <- function(object, newdata, ...) {
  # Remove custom class
  class(object) <- setdiff(class(object), "transitreg_tensorflow")

  X <- model.matrix(attr(object, "formula"), data = newdata)
  X <- tensorflow::tf$convert_to_tensor(X, dtype = "float32")

  probs <- object |> predict(X)

  return(probs)
}

#if(FALSE) {
#  data("WeatherGermany", package = "WeatherGermany")

#  MUC <- subset(WeatherGermany, id == 1262)

#  MUC$year <- as.POSIXlt(MUC$date)$year + 1900
#  MUC$yday <- as.POSIXlt(MUC$date)$yday

#  b <- transitreg(Tmax ~ theta + yday + year, data = MUC, breaks = 200,
#    engine = "tensorflow", ncores = 4, nlayers = 2, units = 200, dropout = 0.01, batch_size = 32,
#    trace = 2)

#  wormplot(b)

#  p <- predict(b, type = "mean")
#  plot(p, MUC$Tmax)
#  abline(0, 1, col = 4, lwd = 2)
#}

