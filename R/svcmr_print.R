#' Print function for svc objects.
#' @description Print function for svc objects.
#' @export
#' @param x An object of type svc.
#' @param ... Not used.
print.svc <- function(x, ...) {
  dimR <- dim(x$R)
  cat("Covariance model:\n")
  print(x$form)
  cat("\nDimension of R:\n")
  cat(dimR[1], "x", dimR[2], "\n")
  cat("\nDistinct elements in off-diagonal:")
  if(methods::is(x$R, "ddiMatrix")) {
    cat("\nR is diagonal")
  } else {
    sumR <- Matrix::summary(x$R)
    sumROff <- table(sumR[sumR$i < sumR$j, "x"])
    print(sumROff)
  }
}

#' Print function for mc objects.
#' @description Print function for mc objects.
#' @export
#' @param x An object of type mc.
#' @param ... Not used
print.mc <- function(x, ...) {
  dimX <- dim(x$X)
  cat("Mean model:\n")
  print(x$form)
  cat("\nDimension of X:\n")
  cat(dimX[1], "x", dimX[2], "\n")
}

#' Print function for summary.fitm
#' @description Print function for summary.fitm objects.
#' @export
#' @param x An object of type summary.svc.
#' @param ... Not used
print.summary.fitm <- function(x, ...) {
  print(format(c("Log likelihood" = x$logl,
                 "Deviance" = x$dev,
                 "AIC" = x$AIC,
                 "BIC" = x$BIC)), quote = F)
  cat("\nFitted parameters:\n")
  stats::printCoefmat(x$est, digits = 3, signif.stars = FALSE, eps.Pvalue = 0.001)
}
