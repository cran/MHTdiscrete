#' The adjusted p-values for Modified Hochberg step-up FWER controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values and all attaianble p-values.
#'
#' @usage
#' MHoch.p.adjust(p,p.set)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param p.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis..
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \eqn{p}).
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' MHoch.p.adjust(p,p.set)
#' @export


MHoch.p.adjust <- function(p,p.set){
  o <- order(p); ro <- order(o); m <- length(p)
  sort.p <- p[o]; sort.p.set <- p.set[o]
  adjP <- numeric(m); pCDF <- matrix(NA,m,m)
  for(j in 1:m){
    for(i in m:1){
      pCDF[i,j] <- max(sort.p.set[[j]][sort.p.set[[j]] <= sort.p[i]],0)
      c <- sum(pCDF[i,i:m])
      adjP[i] <- ifelse(i==m,c,min(adjP[i+1],c))
    }
  }
  return(adjP[ro])
}

#' The number of rejected hypotheses for Roth's step-up FWER controlling procedure.
#'
#' The function for calculating the number of rejected hypotheses (rejection region) based on original available p-values, all attaianble p-values and the given significant level.
#'
#' @usage
#' Roth.rej(p,p.set,alpha)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param p.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis..
#' @param alpha the given significant level for Roth's procedure, the default value is 0.05.
#' @return
#' An integer value of the number of rejected hypotheses.
#' @references
#' Roth, A. J. (1999).
#' Multiple comparison procedures for discrete test statistics.
#' \emph{Journal of statistical planning and inference}, \strong{82}: 101-117.
#'
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' Roth.rej(p,p.set,0.05)
#' @export



Roth.rej <- function(p,p.set,alpha=0.05){
  m <- length(p); minP <- sapply(p.set,min)
  p_R1 <- p[minP < alpha]
  Q <-Inf
  for (j in seq_along(p_R1)){
    if(sort(p_R1,decreasing = T)[j] < alpha/j){
      Q <- j; break
    }
    else{next}
  }
  M <- c()
  for (j in seq_along(p)){
    M[j] <- sum(minP<alpha/j)
    if (M[j]<=j){
      K <- j; break
    }
    else {next}
  }
  p_RK <- if(M[K]<K) c(sort(p[minP < alpha/K],decreasing = T),rep(0,K-M[K])) else      sort(p[minP < alpha/K],decreasing = T)
  W <- Inf
  for (j in seq_along(p_RK)){
     if(max(p_RK[j],p[minP<alpha/j & minP>=alpha/K]) < alpha/j){
      W <- j; break
    }
    else{next}
  }
   return(min(Q,W))
}

#' The adjusted p-values for Roth's step-up FWER controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values and all attaianble p-values.
#'
#' @usage
#' Roth.p.adjust(p,p.set,digits)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param p.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis..
#' @param digits minimal number of significant digits for the adjusted p-values, the default value is 4, see \code{\link[base]{print.default}}.
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \eqn{p}).
#' @references
#' Roth, A. J. (1999).
#' Multiple comparison procedures for discrete test statistics.
#' \emph{Journal of statistical planning and inference}, \strong{82}: 101-117.
#'
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' Roth.p.adjust(p,p.set,digits=5)
#' @export

Roth.p.adjust <- function(p,p.set,digits=4){
  m <- length(p)
  adjP <- c()
  for (i in 1:m){
    ##### setting the accuracy as 4 digitals
    init.alpha <- numeric(); eps <- 10^(-1*c(1:digits))

    init.alpha[1] <- 1
    while(p[i] <= init.alpha[1]/Roth.rej(p,p.set,init.alpha[1]) & init.alpha[1] <= 1){
      init.alpha[1] <- init.alpha[1]-eps[1]
    }

    for (j in 2:digits){
      init.alpha[j] <- init.alpha[j-1]+eps[j-1]
      while(p[i] <=init.alpha[j]/Roth.rej(p,p.set,init.alpha[j]) ){
        init.alpha[j] <- init.alpha[j]-eps[j]
      }
    }
    adjP[i] <- init.alpha[digits]+eps[digits]
  }
  return(adjP)
}

