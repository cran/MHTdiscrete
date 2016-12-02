#' The adjusted p-values for Modified Benjamini-Hochberg (BH) step-up FDR controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values and all attaianble p-values.
#'
#' @usage
#' MBH.p.adjust(p,p.set)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param p.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis.
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}).
#' @seealso \code{\link{MBY.p.adjust}},  \code{\link{MBL.p.adjust}}
#' @author Yalin Zhu
#' @references
#' Benjamini, Y., and Hochberg, Y. (1995).
#'  Controlling the false discovery rate: a practical and powerful approach to multiple testing.
#'  \emph{Journal of the Royal Statistical Society Series B}, \strong{57}: 289-300.
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' MBH.p.adjust(p,p.set)
#' @export

MBH.p.adjust <- function(p,p.set){
  o <- order(p); ro <- order(o); m <- length(p)
  sort.p <- p[o]; sort.p.set <- p.set[o]
  adjP <- numeric();pCDF <- matrix(NA,m,m)
  for(i in m:1){
    for(j in 1:m){
      pCDF[i,j] <- max(sort.p.set[[j]][sort.p.set[[j]] <= sort.p[i]],0)
    }
    c <- sum(pCDF[i,1:m])/i
    adjP[i] <- ifelse(i==m,c,min(adjP[i+1],c))
  }
  return(adjP[ro])
}

#' The adjusted p-values for Gilbert-Tarone-BH step-up FDR controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values and all attaianble p-values.
#'
#' @usage
#' GTBH.p.adjust(p,p.set)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param p.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis.
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}).
#' @seealso \code{\link{GTBY.p.adjust}},  \code{\link{MBH.p.adjust}},  \code{\link{MBY.p.adjust}}
#' @author Yalin Zhu
#' @references
#' Gilbert, P. B. (2005).
#' A modified false discovery rate multiple-comparisons procedure for discrete data, applied to human immunodeficiency virus genetics.
#' \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, \strong{54}: 143-158.
#'
#' Benjamini, Y., and Hochberg, Y. (1995).
#'  Controlling the false discovery rate: a practical and powerful approach to multiple testing.
#'  \emph{Journal of the Royal Statistical Society Series B}, \strong{57}: 289-300.
#'
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' GTBH.p.adjust(p,p.set)
#' @export

GTBH.p.adjust <- function(p,p.set){
  o <- order(p); ro <- order(o); m <- length(p); adjP <- numeric(m)
  sort.p <- p[o]; sort.p.set <- p.set[o]
  minP <- sort(sapply(sort.p.set,min))
  for (j in m:1){
    for (i in 1:m ){
      if(sort.p[j]>=max(minP)){q <- m}
      else if (sort.p[j]>=minP[i] & sort.p[j]<minP[i+1]){q <- i}
    }
    adjP[j] <-  ifelse(j==q,sort.p[j],min(adjP[j+1],q*sort.p[j]/j))
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
#' A numeric vector of the adjusted p-values (of the same length as \code{p}).
#' @seealso  \code{\link{MBH.p.adjust}},  \code{\link{MBL.p.adjust}}
#' @author Yalin Zhu
#' @references
#' Benjamini, Y., and Yekutieli, D. (2001).
#' The control of the false discovery rate in multiple testing under dependency.
#' \emph{Annals of Statistics}, \strong{29}: 1165-1188.
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' MBY.p.adjust(p,p.set)
#' @export


MBY.p.adjust <- function(p,p.set){
  o <- order(p); ro <- order(o); m <- length(p); C <- sum(1/c(1:m))
  sort.p <- p[o]; sort.p.set <- p.set[o]
  adjP <- numeric();pCDF <- matrix(NA,m,m)
  for(i in m:1){
    for(j in 1:m){
      pCDF[i,j] <- max(sort.p.set[[j]][sort.p.set[[j]] <= sort.p[i]],0)
    }
    c <- min(1,sum(pCDF[i,1:m])*C/i)
    adjP[i] <- ifelse(i==m,c,min(adjP[i+1],c))
  }
  return(adjP[ro])
}


#' The adjusted p-values for Gilbert-Tarone-BY step-up FDR controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values and all attaianble p-values.
#'
#' @usage
#' GTBY.p.adjust(p,p.set)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param p.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis.
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}).
#' @seealso \code{\link{GTBH.p.adjust}},  \code{\link{MBH.p.adjust}},  \code{\link{MBY.p.adjust}}
#' @author Yalin Zhu
#' @references
#' Gilbert, P. B. (2005).
#' A modified false discovery rate multiple-comparisons procedure for discrete data, applied to human immunodeficiency virus genetics.
#' \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, \strong{54}: 143-158.
#'
#' Benjamini, Y., and Yekutieli, D. (2001).
#' The control of the false discovery rate in multiple testing under dependency.
#' \emph{Annals of Statistics}, \strong{29}: 1165-1188.
#'
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' GTBY.p.adjust(p,p.set)
#' @export

GTBY.p.adjust <-  function(p,p.set){
  o <- order(p); ro <- order(o); m <- length(p); adjP <- numeric(m)
  sort.p <- p[o]; sort.p.set <- p.set[o]
  minP <- sort(sapply(sort.p.set,min))
  for (j in m:1){
    for (i in 1:m ){
      if(sort.p[j]>=max(minP)){q <- m}
      else if (sort.p[j]>=minP[i] & sort.p[j]<minP[i+1]){q <- i}
    }
    C <- sum(1/c(1:q))
    adjP[j] <-  ifelse(j==q,sort.p[j],min(adjP[j+1],q*C*sort.p[j]/j))
  }
  return(adjP[ro])
}
