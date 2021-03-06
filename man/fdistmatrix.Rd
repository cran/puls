% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fdistmatrix.R
\name{fdistmatrix}
\alias{fdistmatrix}
\title{Distance Between Functional Objects}
\usage{
fdistmatrix(fd, subrange, distmethod)
}
\arguments{
\item{fd}{A functional data object \code{fd} of \code{fda} package.}

\item{subrange}{A vector of two values indicating the value range of
functional object to calculate on.}

\item{distmethod}{The method for calculating the distance matrix. Choose
between \code{"usc"} and \code{"manual"}. \code{"usc"} uses \code{\link[fda.usc:metric.lp]{fda.usc::metric.lp()}}
function while \code{"manual"} uses squared distance between functions. See
Details.}
}
\value{
A distance matrix with diagonal value and the upper half.
}
\description{
Calculate the distance between functional objects over the defined range.
}
\details{
If choosing \code{distmethod = "manual"}, the L2 distance between all pairs of
functions \eqn{y_i(t)} and \eqn{y_j(t)} is given by:
\deqn{d_R(y_i, y_j) = \sqrt{\int_{a_r}^{b_r} [y_i(t) - y_j(t)]^2 dt}.}
}
\examples{
library(fda)
# Examples taken from fda::Data2fd()
data(gait)
# Function only works on two dimensional data
gait <- gait[, 1:5, 1]
gaitbasis3 <- create.fourier.basis(nbasis = 5)
gaitfd3 <- Data2fd(gait, basisobj = gaitbasis3)

fdistmatrix(gaitfd3, c(0.2, 0.4), "usc")
}
