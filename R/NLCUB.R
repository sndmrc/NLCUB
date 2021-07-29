#' Fitting NLCUB models
#'
#' @author Paola Zuccolotto, Marica Manisera, Sandri Marco
#' @param r numeric, observed ratings (either the vector of microdata or the vector of the m observed frequencies (frequency table)
#' @param g	numeric vector, 'latent' categories assigned to each rating point; if \code{g} is declared, pai and xi are estimated for fixed g;if g is not declared, model selection is performed (see method) in order to determine optimal g
#' @param m integer, number of categories of the response scale	(active only when \code{g} is not declared)
#' @param maxT numeric, maximum value for T (must be maxT > m-1, default = 2m-1); active only when \code{g} is not declared
#' @param param0 	numeric, starting values for pai and xi (default: c(0.5,0.5))
#' @param freq.table	logical, if TRUE, the data in r is the vector of the m observed frequencies (frequency table) (default=TRUE)
#' @param method  character, "NM" (likelihood based - Melder-Mead maximization) - "EM" (likelihood based - EM algorithm)
#' @param dk  numeric, proportion of 'don't know' responses; if declared, in addition to the estimate of pai, the estimated of pai adjusted for the presence of dk responses is provided
#' @return  A list with the following estimates:
#' @return * parameter estimates (pai, xi, and g)
#' @return * fitted frequencies
#' @return * diss index
#' @return * average transition probabilities
#' @return * transition probabilities matrix
#' @return * unconditional transition probability
#' @return * expected number of one-rating-point increments
#' @return The command can also display two graphs: observed vs fitted frequencies + transition plot
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
#' @export
#' @importFrom  maxLik maxLik
#' @importFrom  stats weighted.mean
#' @importFrom  stats dbinom
#' @importFrom  stats uniroot
#' @importFrom  stats sd
#' @importFrom  graphics box
#' @importFrom  graphics axis
#' @importFrom  graphics points
#' @importFrom  graphics abline
#' @importFrom  graphics text
#' @importFrom  graphics lines
#' @importFrom  graphics par

NLCUB <- function(r, g=NULL, m=NULL, maxT=NULL, param0=c(0.5,0.5), 
                  freq.table=TRUE, method="EM", dk=NULL) {

  if (!is.null(g)) {
    if (!is.null(m)) {
      warning("Parameter m: unused argument")
    }
    if (!is.null(maxT)) {
      warning("Parameter maxT: unused argument")
    }
    g.est <- g
  }


  if (is.null(m)) {
    if (is.null(g)) {
      stop("Please, specify m")
    } else {
      m <- length(g)
    }
  }

  if (is.null(maxT)) {
    maxT <- 2*m-1
  }


  if (maxT<=(m-1)) {
    stop("Please, specify maxT > m-1")
  }

  if (!freq.table) {
    tabr<-as.matrix(table(r))
    if(length(tabr) < m) {
      rr <- matrix(0,m,1)
      rr[as.integer(row.names(tabr)),] <- tabr
      tabr <- rr
    }
  } else {
    tabr <- as.matrix(r)
  }


  ################################# METHOD NM
  if (method=="NM" & is.null(g)){

    c <- maxT-m+2
    comb<- expand.grid(rep(list(1:c),m))
    selcomb <- which(rowSums(comb)<=(maxT+1))
    comb <- comb[selcomb,]
    ncomb <- nrow(comb)
    max.all <- array(NA,ncomb)
    code.all <- array(NA,ncomb)
    for(i in 1:ncomb) {
      g.i <- as.matrix(comb[i,])
      est.i <- MaxLikNLCUB.fix.g(tabrf=tabr,gf=g.i,param0f=param0)
      max.all[i] <- est.i$maximum
      code.all[i] <- est.i$code
    }
    sel.code1 <- which(code.all!=0)
    max.all[sel.code1] <- NA
    sel <- which(max.all==max(max.all,na.rm=T))[1]
    g.est <- as.matrix(comb[sel,])


  }#end method NM

  ################################# METHOD EM
  if (method=="EM" & is.null(g)) {

    c <- maxT-m+2
    comb<- expand.grid(rep(list(1:c),m))
    selcomb <- which(rowSums(comb)<=(maxT+1))
    comb <- comb[selcomb,]
    ncomb <- nrow(comb)
    max.all <- array(NA,ncomb)
    code.all <- array(NA,ncomb)
    for(i in 1:ncomb){
      g.i <- as.matrix(comb[i,])
      est.i <- EM.NLCUB.fix.g(tabrEM=tabr, gEM=g.i, param0EM=param0, tolerEM=0.000001, maxiterEM=500)
      max.all[i] <- est.i$maximum
    }
    sel <- which(max.all==max(max.all))[1]
    g.est <- as.matrix(comb[sel,])

  }#end method EM

  ################################# ESTIMATION
  if (method=="NM") {
    est <- MaxLikNLCUB.fix.g(tabrf=tabr, gf=g.est, param0f=param0)
  } else if (method=="EM") {
    est <- EM.NLCUB.fix.g(tabrEM=tabr, gEM=g.est, param0EM=param0, tolerEM=0.000001, maxiterEM=500)
  }

  frt <- thfr(pait=est$estimate[1], xit=est$estimate[2], gt=g.est)$Fit
  dissind <- dissNLCUB(rd=tabr, paid=est$estimate[1], xid=est$estimate[2], gd=g.est)
  probs <- probtrans(xipr=est$estimate[2], gpr=g.est, print.matrix=TRUE)
  NL <- NLI(xip3=est$estimate[2], gp3=g.est)

  pai.adj <- est$estimate[1]*(1-dk)

  if (is.null(dk)) {pai.adj <- "Parameter pai has not been adjusted for dk responses"}

  if (method=="NM") {
    out <- list(est,"pai"=est$estimate[1],"csi"=est$estimate[2],"Varmat"=solve(-est$hessian),
                "Infmat"=(-est$hessian)/sum(tabr),"g"=as.vector(g.est),"Fit"=frt,"diss"=dissind,
                "transprob"=probs$transition_probabilities,
                "transprob_mat"=probs$transition_probability_matrix,
                "uncondtransprob"=probs$unconditioned_transition_probability,
                "mu"=probs$expected_number_one.rating.point_increments,"NL_index"=NL,
                "pai_adjusted_for_dk"=pai.adj, tabr=tabr)
  } else if (method=="EM") {
    out <- list(est,"pai"=est$estimate[1],"csi"=est$estimate[2],"g"=as.vector(g.est),"Fit"=frt,"diss"=dissind,
                "transprob"=probs$transition_probabilities, 
                "transprob_mat"=probs$transition_probability_matrix,
                "uncondtransprob"=probs$unconditioned_transition_probability,
                "mu"=probs$expected_number_one.rating.point_increments,"NL_index"=NL,
                "pai_adjusted_for_dk"=pai.adj, tabr=tabr)
  }
  class(out) <- append("NLCUB", class(out))
  return(out)
}


#' @noRd
###################################################################
# LIKELIHOOD FUNCTION
# compute the likelihood function for a NLCUB model
loglik <- function(param,tabrll,gll) {
  pai <- param[1]
  xi <- param[2]
  mll <- length(gll)
  sg <- c(0,cumsum(gll))
  w <- matrix(0,mll,sum(gll))
  for(i in 1:mll){
    w[i,((sg[i]+1):sg[i+1])] <- 1
  }
  R <- c(1:sum(gll))
  SBinom <- w%*%((choose(sum(gll)-1,R-1)*((1-xi)^(R-1))*(xi^(sum(gll)-R))))
  ll <- sum(tabrll*(log(pai*(SBinom)+(1-pai)*(1/mll))))
  ll
}

#' @noRd
#######################################################################
#MAXIMUM LIKELIHOOD ESTIMATON
#######################################################################
#######################################################################
# NELDER-MEAD MAXIMIZATION
# with fixed g
# inputs:
# tabrf 	= vector of the m observed frequencies (frequency table)
# gf		= vector of the 'latent' categories assigned to each rating point
# param0f 	= starting values for pai and xi
# output: Maximum Likelihood estimates of pai and csi, in this order
MaxLikNLCUB.fix.g <- function(tabrf,gf,param0f) {
  N <- sum(tabrf)
  A <- rbind(diag(2),-diag(2))
  B <- rbind(matrix(0,2,1),matrix(1,2,1))
  est <- maxLik(loglik, start=param0f,tabrll=tabrf,gll=gf,
                constraints=list(ineqA=A,ineqB=B)
  )
  est
}

#' @noRd
#######################################################################
# EM ALGORITHM
# with fixed g
# inputs:
# tabrEM 	= vector of the m observed frequencies (frequency table)
# gEM		= vector of the 'latent' categories assigned to each rating point
# param0EM 	= starting values for pai and xi
# tolerEM 	= tolerance
# output: Maximum Likelihood estimates of pai and csi, in this order, obtained with EM algorithm
EM.NLCUB.fix.g <- function(tabrEM,gEM,param0EM,tolerEM,maxiterEM) {
  N <- sum(tabrEM)
  G <- c(0,cumsum(gEM))
  m <- length(gEM)
  T <- sum(gEM) - 1
  pai <- param0EM[1]
  xi <- param0EM[2]
  k <- 0
  diffpai <- 1
  diffxi <- 1
  while (diffpai > tolerEM | diffxi > tolerEM){
    k <- k+1
    Br <- thfr(pai[k],xi[k],gEM)$Br
    pt <- thfr(pai[k],xi[k],gEM)$Fit
    # E-step
    tau <- (pai[k]*Br)/pt
    # M-step
    pai <- c(pai,weighted.mean(tau,tabrEM))
    f <- function(xi){
      weighted.mean(tau/thfr(pai[k],xi,gEM)$Br *(dbinom(G[2:(m+1)]-1, T-1, 1-xi)-dbinom(G[1:m]-1, T-1, 1-xi)) ,tabrEM)
    }
    xi <- c(xi,uniroot(f,interval=c(10^-15,1-10^-15))$root)
    diffpai <- abs(pai[k+1]-pai[k])
    diffxi <- abs(xi[k+1]-xi[k])
    check <- "normal convergence"
    if (k==maxiterEM){
      diffpai <- 0
      diffxi <- 0
      check <- "iterations stopped (maxiter)"
    }
  } #end while

  pailist<-pai
  xilist<-xi
  pai <- pailist[k+1]
  xi <- xilist[k+1]

  Br <- thfr(pai,xi,gEM)$Br

  Ipai <- 1/N * sum (tabrEM*(((Br-(1/m))^2)/(thfr(pai,xi,gEM)$Fit^2)))

  A1 <- (pai*T*(dbinom(G[2:(m+1)]-1, T-1, 1-xi)-dbinom(G[1:m]-1, T-1, 1-xi))^2)/(thfr(pai,xi,gEM)$Fit^2)
  A2 <- ((T-1) * (dbinom(G[2:(m+1)]-1, T-2, 1-xi)- dbinom(G[2:(m+1)]-2, T-2, 1-xi)- dbinom(G[1:m]-1, T-2, 1-xi) + dbinom(G[1:m]-2, T-2, 1-xi)))/(thfr(pai,xi,gEM)$Fit)

  Ixi <- T/N *pai*sum(tabrEM*(A1-A2))

  B1 <- sum ( (pai * tabrEM * T * (dbinom(G[2:(m+1)]-1, T-1, 1-xi)-dbinom(G[1:m]-1, T-1, 1-xi))* (Br-(1/m)) )/ (thfr(pai,xi,gEM)$Fit^2)  )
  B2 <- sum ( (tabrEM * T * (dbinom(G[2:(m+1)]-1, T-1, 1-xi)-dbinom(G[1:m]-1, T-1, 1-xi)))/ (thfr(pai,xi,gEM)$Fit)  )

  Ipaixi <- (-T)/(m*N)*sum(tabrEM*((dbinom(G[2:(m+1)]-1, T-1, 1-xi)-dbinom(G[1:m]-1, T-1, 1-xi)))/(thfr(pai,xi,gEM)$Fit^2))

  I<-matrix(c(Ipai,Ipaixi,Ipaixi,Ixi),2,2)

  V <- 1/N*solve(I)

  ll <- loglik(c(pai,xi),tabrEM,gEM)

  list("pailist"=pailist,"xilist"=xilist,"estimate"=c(pai,xi),"maximum"=ll,"k"=k,"check"=check,"Varmat"=V,"Infmat"=I)
}



#' @noRd
#####################################################################
# FITTED FREQUENCIES
# inputs:
# pait 		= value of the pai parameter
# xit		= value of the csi parameter
# gt		= vector of the 'latent' categories assigned to each rating point
# output: m fitted frequencies according to the NLCUB models with parameters pai, csi, g
thfr <- function(pait,xit,gt) {
  m <- length(gt)
  sg <- c(0,cumsum(gt))
  R <- c(1:sum(gt))
  w <- matrix(0,m,sum(gt))
  for(i in 1:m) {
    w[i,((sg[R[i]]+1):sg[R[i]+1])] <- 1
  }
  SBinom <- w%*%((choose(sum(gt)-1,R-1)*((1-xit)^(R-1))*(xit^(sum(gt)-R))))
  tf <- pait*(SBinom)+(1-pait)*(1/m)
  list("Fit"=tf, "Br"=SBinom)
}


#' @noRd
#####################################################################
# DISS INDEX
# inputs:
# rd 		     = observed ratings (either the vector of microdata or the vector of the m observed frequencies (frequency table)
# paid 		     = value of the pai parameter
# xid		     = value of the csi parameter
# gd		     = vector of the 'latent' categories assigned to each rating point
# freq.table.diss    = logical flag: if TRUE, the data in r is the vector of the m observed frequencies (frequency table)
# output: value of the dissimilarity index
dissNLCUB <- function(rd,paid,xid,gd,freq.table.diss=TRUE) {
  m <- length(gd)
  if (freq.table.diss==FALSE){
    tabr<-as.matrix(table(rd))
    if(length(tabr) < m){
      rr<-matrix(0,m,1)
      rr[as.integer(row.names(tabr)),]<-tabr
      tabr<-rr
    }
  }else{tabr<-as.matrix(rd)}

  fr <- tabr/sum(tabr)
  diss <- 1/2*sum(abs(fr-thfr(pait=paid,xit=xid,gt=gd)$Fit))
  #cat("Dissimilarity index=",diss)
  diss
}

#' @noRd
#####################################################################
# OBSERVED VS. FITTED FREQUENCIES PLOT
# inputs:
# rp1 		     = observed ratings (either the vector of microdata or the vector of the m observed frequencies (frequency table)
# paip1		     = value of the pai parameter
# xip1		     = value of the csi parameter
# gp1		     = vector of the 'latent' categories assigned to each rating point
# freq.table.p1      = logical flag: if TRUE, the data in r is the vector of the m observed frequencies (frequency table)
# output:
# the observed vs. fitted frequencies plot
NLCUBplot2 <- function(rp1,paip1,xip1,gp1,freq.table.p1=TRUE){
  m <- length(gp1)
  if (!freq.table.p1) {
    tabr < -as.matrix(table(rp1))
    if(length(tabr) < m) {
      rr <- matrix(0,m,1)
      rr[as.integer(row.names(tabr)),]<-tabr
      tabr <- rr
    }
  } else {
    tabr <- as.matrix(rp1)
  }

  par(mar=c(5.1, 5.1, 4.1, 1.1))
  plot(tabr/sum(tabr),axes=FALSE,main=paste("Diss = ",round(dissNLCUB(rd=tabr,paid=paip1,xid=xip1,gd=gp1),4),sep=""),
       xlab='Ratings',pch=16,cex.lab=0.8,
       ylab='Observed relative frequencies (dots) and\nfitted probabilities (circles)',ylim=c(0,1))
  box()
  axis(1,1:m)
  axis(2)
  points(thfr(pait=paip1,xit=xip1,gt=gp1)$Fit, type='b', lty='dashed', cex=1.4)
  abline(c(0,0),c(0,0))
}

#' @noRd
#####################################################################
# TRANSITION PROBABILITIES, UNCONDITIONED TRANSITION PROBABILITY, EXPECTED NO. OF ONE-RATING-POINT INCREMENTS
# inputs:
# xipr		 = value of the csi parameter
# gpr		 = vector of the 'latent' categories assigned to each rating point
# print.matrix	 = logical flag: if TRUE, the matrix of the transition probabilities for all t,s is given as an output
# output: a list containing
#  $transition_probability_matrix (if print.matrix=TRUE)
#  $transition_probabilities
#  $unconditioned_transition_probability
#  $expected_number_one.rating.point_increments
probtrans <- function(xipr,gpr,print.matrix=FALSE) {
  m <- length(gpr)
  n <- sum(gpr)-1
  br <- c(0,cumsum(gpr))-1

  probtrans.i <- matrix(0,n-1,m-1)
  probcong.i <- matrix(0,n-1,m-1)

  ######## probtrans from t=1 to t=last but one step
  ######## compute the probability of increase in the next step
  for (i in 1:(n-1)) {
    for (r in 2:m) {
      num <- dbinom(br[r],i,1-xipr)
      den <- sum(dbinom((br[r-1]+1):br[r],i,1-xipr))
      probtrans.i[i,r-1] <- (1-xipr)*num/den
      probcong.i[i,r-1] <- (1-xipr)*num
    }
  }

  ######## adding Step 0
  probtrans.i <- rbind(matrix(NaN,1,m-1),probtrans.i)
  probtrans.i[1,1] <- 0
  if(br[2]==0) probtrans.i[1,1] <- 1-xipr

  prob.trans <- colMeans(probtrans.i, na.rm=T)

  probcong.i <- rbind(matrix(0,1,m-1),probcong.i)
  if(br[2]==0) probcong.i[1,1] <- 1-xipr

  prob.trans.unc <- mean(rowSums(probcong.i))

  mi <- prob.trans.unc*n #expected number of one-rating-point increments

  if (print.matrix) {
    return(list(transition_probability_matrix=probtrans.i,
                transition_probabilities=prob.trans,
                unconditioned_transition_probability=prob.trans.unc,
                expected_number_one.rating.point_increments=mi))
  } else {
    return(list(transition_probabilities=prob.trans,
                unconditioned_transition_probability=prob.trans.unc,
                expected_number_one.rating.point_increments=mi))
  }
}


#' @noRd
#####################################################################
# TRANSITION PLOT
# inputs:
# xip2		= value of the csi parameter
# gp2		= vector of the 'latent' categories assigned to each rating point
# log.scale	= logical flag: if TRUE, the logarithmic scale is used (if FALSE, the linear scale is used)
# output:
# the transition plot
transplot2 <- function(xip2,gp2,log.scale=TRUE) {
  m <- length(gp2)
  prob.trans <- probtrans(xipr=xip2,gpr=gp2)$transition_probabilities

  if (log.scale==TRUE){
    prob.trans <- -log(prob.trans)}

  a <- rbind(0,as.matrix(prob.trans))
  a <- cumsum(a)
  a <- a/max(a)

  plot(a,type='b',
       ylim=c(-0.1,max(a)+0.1),
       lty=1,lwd=4,
       xlim=c(0.7,m+0.2),xlab='Ratings',
       xaxt='n',yaxt='n',ylab='Perceived ratings',
       cex.lab=0.8,bty='l')

  text(c(1:m),array(-0.1,m), labels=c(1:m))
  tx <- a
  text(array(0.8,m),tx, labels=paste("NLCUB",c(1:m),sep=""),pos=3,cex=0.6)

  for(i in 1:m){
    lines(c(i,i),c(-0.1,tx[i]),lty='dotted')
  }

  for(i in 1:m){
    lines(c(0.7,i),c(tx[i],tx[i]),lty='dotted')
  }
}




