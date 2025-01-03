\name{tmdist}
\alias{tmdist}
\alias{random.tmdist}
\alias{pdf.tmdist}
\alias{log_pdf.tmdist}
\alias{cdf.tmdist}
\alias{quantile.tmdist}
\title{Create an tmdist Distribution}
\usage{
tmdist(probs, bins)

\method{as.data.frame}{tmdist}(x, ...)
\method{random}{tmdist}(x, n = 1L, drop = TRUE, ...)
\method{pdf}{tmdist}(d, x, drop = TRUE, elementwise = NULL, ncores = NULL, ...)
\method{log_pdf}{tmdist}(d, x, drop = TRUE, elementwise = NULL, ncores = NULL, ...)
\method{cdf}{tmdist}(d, x, drop = TRUE, elementwise = NULL, ncores = NULL, ...)
\method{quantile}{tmdist}(x, probs, drop = TRUE, elementwise = NULL, ncores = NULL, ...)
}
\arguments{
    \item{x}{object of class \code{tmdist} (if main argument), else a numeric
        vector at which the CDF, PDF, or quantile should be evaluated.}
    \item{probs}{Numeric vector of transition probabilities}
    \item{bins}{Center of the (pseudo-)bins corresponding to the transition
        probabilites. Numeric, same length as \code{probs}.}
    \item{...}{Currently not used.}
    \item{n}{The number of samples to draw. Defaults to `1L`.}
    \item{drop}{logical. Should the result be simplified to a vector if possible?}
    \item{d}{An `tmdist` object created by a call to [tmdist()].}
    \item{elementwise}{logical. Should each distribution in \code{x} be evaluated
        at all elements of \code{probs} (\code{elementwise = FALSE}, yielding a matrix)?
        Or, if \code{x} and \code{probs} have the same length, should the evaluation be
        done element by element (\code{elementwise = TRUE}, yielding a vector)? The
        default of \code{NULL} means that \code{elementwise = TRUE} is used if the
        lengths match and otherwise \code{elementwise = FALSE} is used.}
    \item{ncores}{number of cores to use with OpenMP. If \code{NULL} (default)
        all but two available cores will be used if OpenMP is supported.}
}
\value{
A `tmdist` object.

TODO(R): Missing.
}
\description{
TODO(R): Missing
}
\details{
TODO(R): Missing
}
\examples{
data <- structure(list(bin = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
    12, 13, 14, 15, 16, 17, 18, 19), tp = c(0.998987937051858, 0.997723761984778,
    0.994959741987477, 0.989144986733817, 0.977718308296039, 0.958030163744388,
    0.929870939128467, 0.895691370872762, 0.856730189019159, 0.812772393792261,
    0.767294535298433, 0.727450916579927, 0.695208938905971, 0.662724982750024,
    0.619413076608789, 0.561464014434395, 0.493907965967658, 0.424410082007259,
    0.357570963995397, 0.295847728009519), binmid = c(-3.16187734522192,
    -2.83533036950331, -2.5087833937847, -2.18223641806609, -1.85568944234748,
    -1.52914246662887, -1.20259549091026, -0.876048515191646, -0.549501539473035,
    -0.222954563754424, 0.103592411964187, 0.430139387682798, 0.75668636340141,
    1.08323333912002, 1.40978031483863, 1.73632729055724, 2.06287426627585,
    2.38942124199446, 2.71596821771308, 3.04251519343169)), row.names = c(NA,
    20L), class = "data.frame")

X <- tmdist(data$tp, data$binmid)
X

mean(X)

random(X, 10)

pdf(X, 2)
log_pdf(X, 2)

cdf(X, 4)
quantile(X, 0.7)
}
\concept{tmdist distribution}
\concept{Normal distribution}
\concept{distributions}
