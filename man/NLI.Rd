% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NLI.R
\name{NLI}
\alias{NLI}
\title{Estimate NL Index}
\usage{
NLI(xip3, gp3)
}
\arguments{
\item{xip3}{numeric, value of the csi parameter}

\item{gp3}{numeric, vector of the 'latent' categories assigned to each rating point}
}
\value{
The NL index value
}
\description{
Estimate NL Index
}
\details{
(Details here)
}
\examples{
N <- 1000
pai.sim <- 0.8
xi.sim <- 0.3
g.sim <- c(1,1,2,4,2)
cats <- 5
set.seed(1234567)
}
\references{
M. Manisera and P. Zuccolotto (2014) Modeling rating data with Nonlinear CUB models. Computational Statistics and Data Analysis, 78, pp. 100–118

M. Manisera and P. Zuccolotto (2014) Nonlinear CUB models: The R code. Statistical Software - Statistica & Applicazioni, Vol. XII, n. 2, pp. 205-223
}
\author{
Paola Zuccolotto, Marica Manisera, Sandri Marco
}
