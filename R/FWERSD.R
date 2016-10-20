#' The adjusted p-values for Modified Holm step-down FWER controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values and all attaianble p-values.
#'
#' @usage
#' MHolm.p.adjust(p,p.set)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param p.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis.
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \eqn{p}).
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' MHolm.p.adjust(p,p.set)
#' @export

MHolm.p.adjust <- function(p,p.set){
  o <- order(p); ro <- order(o); m <- length(p)
  sort.p <- p[o]; sort.p.set <- p.set[o]
  adjP <- numeric(m); pCDF <- matrix(NA,m,m)
  for(j in 1:m){
    for(i in 1:m){
      pCDF[i,j] <- max(sort.p.set[[j]][sort.p.set[[j]] <= sort.p[i]],0)
      c <- min(1,sum(pCDF[i,i:m]))
      adjP[i] <- ifelse(i==1,c,max(adjP[i-1],c))
    }
  }
  return(adjP[ro])
}

#' The adjusted p-values for Tarone-Holm step-down FWER controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values and all attaianble p-values.
#'
#' @usage
#' TH.p.adjust(p,p.set)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param p.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis.
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \eqn{p}).
#' @references
#' Hommel, G., & Krummenauer, F. (1998).
#' Improvements and modifications of Tarone's multiple test procedure for discrete data.
#' \emph{Biometrics}, \strong{54}: 673-681.
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' TH.p.adjust(p,p.set)
#' @export


TH.p.adjust <- function(p,p.set){
  o <- order(p); ro <- order(o); m <- length(p)
  sort.p <- p[o]; sort.p.set <- p.set[o]; adjP <- numeric(m)
  j <- 1
  while (j <= m){
    minP <- sort(sapply(sort.p.set[j:m],min))
    for (i in 1:(m-j+1)){
      if (sort.p[j]>=max(minP)){q <- m-j+1}
      else if (sort.p[j]>=minP[i] & sort.p[j]<minP[i+1]){q <- i}
    }
    c <- min(1,q*sort.p[j])
    adjP[j] <- ifelse(j==1,c,max(adjP[j-1],c))
    j <- j+1
  }
  return(adjP[ro])
}


