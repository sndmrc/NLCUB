#' Plot observed relative frequencies nad fitted probabilities for different ratings
#'
#' @author Paola Zuccolotto, Marica Manisera, Sandri Marco
#' @param rp1  numeric observed ratings (either the vector of microdata or the vector of the m observed frequencies (frequency table)
#' @param paip1	numeric, value of the pai parameter
#' @param xip1	numeric, value of the csi parameter
#' @param gp1		numeric vector, 'latent' categories assigned to each rating point
#' @param freq.table.p1 logical, if TRUE, the data in r is the vector of the m observed frequencies (frequency table)
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
#' @export
#' @importFrom  ggplot2 ggplot
#' @importFrom  ggplot2 aes
#' @importFrom  ggplot2 geom_point
#' @importFrom  ggplot2 geom_line
#' @importFrom  ggplot2 geom_hline
#' @importFrom  ggplot2 geom_vline
#' @importFrom  ggplot2 labs
#' @importFrom  ggplot2 ylim
#' @importFrom  ggplot2 theme_bw

NLCUBplot <- function(rp1, paip1, xip1, gp1, freq.table.p1=TRUE) {
  ratings <- NULL;  obs_rel_freq <- NULL; fit_probs <- NULL;
  m <- length(gp1)
  if (!freq.table.p1) {
    tabr < -as.matrix(table(rp1))
    if(length(tabr) < m) {
      rr <- matrix(0,m,1)
      rr[as.integer(row.names(tabr)),] <- tabr
      tabr <- rr
    }
  } else {
    tabr <- as.matrix(rp1)
  }
  df <- data.frame(ratings = 1:m, 
                   obs_rel_freq = tabr/sum(tabr),
                   fit_probs = thfr(pait=paip1,xit=xip1,gt=gp1)$Fit)
  
  p <- ggplot(data=df, aes(x=ratings, y=obs_rel_freq)) +
    geom_point(size=2) +
    geom_point(aes(y = fit_probs), shape=1, size=3) +
    geom_line(aes(y = fit_probs), linetype=2, size=1) +
    geom_hline(yintercept=0) +
    ylim(c(0,1)) + 
    labs(x = 'Ratings', y = 'Observed relative frequencies (dots) and\nfitted probabilities (circles)') +
    theme_bw()
  
  return(p)
}