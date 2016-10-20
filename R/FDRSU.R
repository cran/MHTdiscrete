#' The adjusted p-values for Modified Benjamini-Hochberg (BH) step-up FDR controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values and all attaianble p-values. 
#'
#' @usage
#' MBH.p.adjust(p,p.set)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param p.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis.
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \eqn{p}).
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' MBH.p.adjust(p,p.set)
#' @export

MBH.p.adjust <- function(p,p.set){
      o <- order(p); ro <- order(o); m <- length(p)
      sort.p <- p[o]; sort.p.set <- p.set[o]
      adjP <- numeric();pCDF <- matrix(NA,m,m)
      for(j in 1:m){
        for(i in m:1){
          pCDF[i,j] <- max(sort.p.set[[j]][sort.p.set[[j]] <= sort.p[i]],0)
          c <- sum(pCDF[i,1:m])/i
          adjP[i] <- ifelse(i==m,c,min(adjP[i+1],c))
        }
      }
      return(adjP[ro])
    }
    
#' The adjusted p-values for Modified Benjamini-Yekutieli (BY) step-up FDR controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values and all attaianble p-values.
#'
#' @usage
#' MBY.p.adjust(p,p.set)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param p.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis.
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \eqn{p}).
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' MBY.p.adjust(p,p.set)
#' @export


MBY.p.adjust <- function(p,p.set){
       o <- order(p); ro <- order(o); m <- length(p);C <- sum(1/c(1:m))
      sort.p <- p[o]; sort.p.set <- p.set[o]
      adjP <- numeric();pCDF <- matrix(NA,m,m)
      for(j in 1:m){
        for(i in m:1){
          pCDF[i,j] <- max(sort.p.set[[j]][sort.p.set[[j]] <= sort.p[i]],0)
          c <- min(1,sum(pCDF[i,1:m])*C/i)
          adjP[i] <- ifelse(i==m,c,min(adjP[i+1],c))
        }
      }
      return(adjP[ro])
    }




