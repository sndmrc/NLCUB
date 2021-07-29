#' Transition plot
#'
#' @author Paola Zuccolotto, Marica Manisera, Sandri Marco
#' @param xip2  numeric, value of the csi parameter
#' @param gp2	numeric vector, 'latent' categories assigned to each rating point
#' @param log.scale logical, if TRUE, the logarithmic scale is used; if FALSE, the linear scale is used
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
#' @importFrom  ggplot2 geom_text
#' @importFrom  ggplot2  xlim
#' @importFrom  ggplot2  scale_y_continuous 
#' @importFrom  ggplot2  theme 
#' @importFrom  ggplot2  element_blank 
#' @importFrom  ggplot2  coord_cartesian 
#' @importFrom  ggplot2  %+replace% 

Transplot <- function(xip2, gp2, log.scale=TRUE) {
  x <- NULL; y <- NULL; labs <- NULL; grp=NULL;
  m <- length(gp2)
  prob.trans <- probtrans(xipr=xip2, gpr=gp2)$transition_probabilities
  
  if (log.scale){
    prob.trans <- -log(prob.trans)
  }
  a <- rbind(0, as.matrix(prob.trans))
  a <- cumsum(a)
  a <- a/max(a)
  
  df <- data.frame(x = 1:m, y = a, 
                  labs=paste0("NLCUB",c(1:m)))
  df_vline <- data.frame(x=rep(1:m,2), y=c(rep(-0.05,m),df$y), grp=factor(rep(1:m, 2)))
  df_hline <- data.frame(x=c(rep(0.5,m),df$x), y=rep(df$y,2), grp=factor(rep(1:m, 2)))  

  p <- ggplot(data=df, aes(x=x, y=y)) +
    geom_point(size=2) +
    geom_line(size=1) +
    geom_text(aes(label=labs), x=0.5, hjust=-0.2, vjust=-0.2) +
    geom_line(data=df_hline, aes(x=x, y=y, group=grp), linetype=2, inherit.aes=F) +
    geom_line(data=df_vline, aes(x=x, y=y, group=grp), linetype=2, inherit.aes=F) +    
    scale_y_continuous(breaks = round(df$y,3), minor_breaks=NULL, limits=c(-.05,1.05)) +
    coord_cartesian(expand=FALSE) +
    xlim(c(0.5,m+0.05)) +
    labs(x='Ratings', y='Perceived ratings') +
    theme_bw() %+replace% theme(panel.grid = element_blank())
  
  return(p)
}