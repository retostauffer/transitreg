\name{transitreg-package}
\alias{transitreg-package}

\docType{package}

\title{
  Flexible Estimation and Visualization of Transition Models
}

\description{
  The \pkg{transitreg} package provides tools for estimating and visualizing 
  transition models, a flexible framework for modeling count or continuous responses 
  without assuming fixed distributional forms. These models allow for dynamic 
  probabilistic learning using binary GLM-type approaches and can handle scenarios 
  such as zero-inflated or overdispersed data.
}

\details{
  \pkg{transitreg} supports the following key features:
  \itemize{
    \item \code{\link{transitreg}}: Fit transition models to count or continuous response data, 
      leveraging generalized additive model (GAM) techniques for estimation.
    \item \code{\link{transitreg_data}}: Transform raw data into the format required for 
      transition model estimation.
    \item \code{\link{transitreg_dist}}: Visualize probability density functions (PDFs) and 
      raw count data modeled with transition models.
  }

  By allowing the response distribution to be data-driven and nonparametric, 
  the package overcomes limitations of fixed distributional assumptions, such as those 
  in Poisson or negative binomial models.

  Additional features include:
  \itemize{
    \item Support for zero-inflated and hurdle-type models.
    \item Integration with \pkg{mgcv} for efficient GAM-based estimation.
    \item Support for both count and continuous response variables, with flexible 
      options for discretizing continuous responses.
  }
}

\author{
  Developed by Nikolaus Umlauf and Reto Stauffer.

  Maintainer: Nikolaus Umlauf (\email{Nikolaus.Umlauf@uibk.ac.at})
}

\references{
  Berger, M., & Tutz, G. (2021). Transition Models for Count Data: A 
  Flexible Alternative to Fixed Distribution Models. Statistical 
  Methods & Applications, 30, 1259--1283. \doi{10.1007/s10260-021-00558-6}.
}

\seealso{
  \code{\link{transitreg}}, \code{\link{transitreg_data}}, \code{\link{transitreg_dist}}
}

\keyword{package}

