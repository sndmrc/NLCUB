#' Plots a NLCUB model from a 'NLCUB' object
#'
#' @author Paola Zuccolotto, Marica Manisera, Sandri Marco
#' @param x an object of class \code{NLCUB}
#' @param log.scale logical, if TRUE, the logarithmic scale is used; if FALSE, the linear scale is used
#' @param freq.table.p1 logical, if TRUE, the data in r is the vector of the m observed frequencies (frequency table)
#' @param ... other graphical parameters.
#' @details (Details here)
#' @return  A ggplot object
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
#' datitab <- as.matrix(table(dataNLCUB))
#' est <- NLCUB(datitab, g=g.sim, freq.table=TRUE)
#' plot(est)
#' @method plot NLCUB
#' @export
#' @importFrom  gridExtra grid.arrange

plot.NLCUB <- function(x, log.scale=TRUE, freq.table.p1=TRUE, ...) {
  
  p1 <- NLCUBplot(rp1=x$tabr, paip1=x$pai, xip1=x$csi, gp1=x$g, freq.table.p1=freq.table.p1)
  p2 <- Transplot(xip2=x$csi, gp2=x$g, log.scale=log.scale)
  p <- grid.arrange(p1, p2, nrow=1)
  invisible(p)
}