% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tensorflow.R
\name{transitreg_tensorflow}
\alias{transitreg_tensorflow}
\alias{predict.transitreg_tensorflow}
\title{Estimat Transition Probabilities via \code{tensorflow}/\code{keras3}}
\usage{
transitreg_tensorflow(
  formula,
  data,
  nlayers = 2,
  units = 20,
  epochs = 1000,
  batch_size = 16,
  activation = "relu",
  dropout = 0.1,
  validation_split = 0.2,
  patience = 10,
  trace = 0,
  ncores = 1,
  ...
)

\method{predict}{transitreg_tensorflow}(object, newdata, ...)
}
\arguments{
\item{formula}{An object of class \code{formula}.}

\item{data}{A data.frame containing the required data to fit a tensorflow
model (given \code{formula}) using \code{\link[keras3:keras_model_sequential]{keras3::keras_model_sequential()}}.}

\item{nlayers}{Number of layers.}

\item{units}{... (TODO)}

\item{epochs}{... (TODO)}

\item{batch_size}{... (TODO)}

\item{activation}{... (TODO)}

\item{dropout}{... (TODO)}

\item{validation_split}{... (TODO)}

\item{patience}{... (TODO)}

\item{trace}{... (TODO)}

\item{ncores}{Number of cores to be used to train the model.}

\item{\dots}{currently unused.}

\item{object}{an object of class \code{transitreg_tensorflow}.}

\item{newdata}{Data frame with the variables required for the prediction.}
}
\description{
Experimental engine used to estimate the transition probabilities
in a \code{\link[=transitreg]{transitreg()}} model.
}
\author{
Niki
}
