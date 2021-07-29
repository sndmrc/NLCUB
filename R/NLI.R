#' Estimate NL Index
#'
#' @author Paola Zuccolotto, Marica Manisera, Sandri Marco
#' @param xip3  numeric, value of the csi parameter
#' @param gp3 numeric, vector of the 'latent' categories assigned to each rating point
#' @details (Details here)
#' @return  The NL index value
#' @references M. Manisera and P. Zuccolotto (2014) Modeling rating data with Nonlinear CUB models. Computational Statistics and Data Analysis, 78, pp. 100–118
#' @references M. Manisera and P. Zuccolotto (2014) Nonlinear CUB models: The R code. Statistical Software - Statistica & Applicazioni, Vol. XII, n. 2, pp. 205-223
#' @examples
#' N <- 1000
#' pai.sim <- 0.8
#' xi.sim <- 0.3
#' g.sim <- c(1,1,2,4,2)
#' cats <- 5
#' set.seed(1234567)
#' @export
NLI <- function(xip3, gp3) {
  m <- length(gp3)
  pt.matrix <- probtrans(xipr=xip3, gpr=gp3, print.matrix=TRUE)$transition_probability_matrix
  
  card <- sum(is.finite(pt.matrix))
  if (m %% 2 != 0) {max.NLI <- sqrt(1/4)} #odd 
  if (m %% 2 == 0) {max.NLI <- sqrt( 1/4 - 1/(4* card^2))} #even
  
  sigma <- sd(pt.matrix,na.rm=TRUE)*(card-1)/card
  NLindex <- sigma/max.NLI
  return(NLindex)
}