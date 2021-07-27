#' Simulate ratings according to a NLCUB model
#'
#' @author Paola Zuccolotto, Marica Manisera, Sandri Marco
#' @param Nsim  integer, sample size
#' @param paisim numeric, value of the \code{pai} parameter to be used for simulation
#' @param xisim	 numeric, value of the \code{csi} parameter to be used for simulation
#' @param gsim numeric vector, the 'latent' categories assigned to each rating point
#' @details It is necessary that the name of the team is contained in the column corresponding to the description
#' @return  A data frame with N ratings simulated according to a NLCUB model with parameters \code{paisim}, \code{csisim}, \code{gsim}
#' @references M. Manisera and P. Zuccolotto (2014) Modeling rating data with Nonlinear CUB models. Computational Statistics and Data Analysis, 78, pp. 100–118
#' @references M. Manisera and P. Zuccolotto (2014) Nonlinear CUB models: The R code. Statistical Software - Statistica & Applicazioni, Vol. XII, n. 2, pp. 205-223
#' @examples
#' N <- 1000
#' pai.sim <- 0.8
#' xi.sim <- 0.3
#' g.sim <- c(1,1,2,4,2)
#' cats <- 5
#' set.seed(1234567)
#' dataNLCUB <- simNLCUB(N, pai.sim, xi.sim, g.sim)
#' @export
#' @importFrom  stats rbinom
#' @importFrom  stats runif

simNLCUB <- function(Nsim, paisim, xisim, gsim) {
  m <- length(gsim)
  br.sim <- c(0, cumsum(gsim))
  Wn <- rbinom(Nsim, sum(gsim)-1, 1-xisim)+1
  Wn <- cut(Wn, breaks=br.sim, labels=FALSE)
  Un <- sample.int(m, size = Nsim, replace = T)
  rand <- runif(Nsim, min=0, max=1)
  r <- Wn
  for (i in 1:Nsim) {
    if (rand[i] > paisim) r[i] <- Un[i]
  }
  return(r)
}
