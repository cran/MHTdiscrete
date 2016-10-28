#' The adjusted p-values for Modified Benjamini-Liu (BL) step-down FDR controlling procedure.
#'
#' The function for calculating the adjusted p-values based on original available p-values and all attaianble p-values.
#'
#' @usage
#' MBL.p.adjust(p,p.set)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param p.set a list of numeric vectors, where each vector is the vector of all attainable p-values containing the available p-value for the corresponding hypothesis.
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \eqn{p}).
#' @author Yalin Zhu
#' @references
#' Benjamini, Y., and Liu, W. (1999).
#' A step-down multiple hypotheses testing procedure that controls the false discovery rate under independence.
#'  \emph{Journal of Statistical Planning and Inference}, \strong{82}: 163-170.
#' @note The MBL procedure for discrete data controls FDR under the specific dependence assumption where the joint distribution of statistics from true nulls are independent of the joint distribution of statistics from false nulls.
#' @examples
#' p <- c(pbinom(1,8,0.5),pbinom(1,5,0.75),pbinom(1,6,0.6))
#' p.set <-list(pbinom(0:8,8,0.5),pbinom(0:5,5,0.75),pbinom(0:6,6,0.6))
#' MBL.p.adjust(p,p.set)
#' @export

MBL.p.adjust <- function(p,p.set){
  o <- order(p); ro <- order(o); m <- length(p)
  sort.p <- p[o]; sort.p.set <- p.set[o]
  adjP <- numeric(m); pCDF <- matrix(NA,m,m)
  for(i in 1:m){
    for(j in i:m){
      pCDF[i,j] <- max(sort.p.set[[j]][sort.p.set[[j]] <= sort.p[i]],0)
    }
    c <- (m-i+1)/m*(1-prod(1-pCDF[i,i:m]))
    adjP[i] <- ifelse(i==1,c,max(adjP[i-1],c))
  }
  return(adjP[ro])
}
