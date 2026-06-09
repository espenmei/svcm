#' Print function for svc objects.
#' @description Print function for svc objects.
#' @export
#' @param x An object of type svc.
#' @param ... Not used.
print.svc <- function(x, ...) {

  cat("Covariance model:\n")
  print(x$form)

  if(is.null(x$R)) {
    return(invisible(x))
  }
  dimR <- dim(x$R)
  cat("\nDimension of R:\n")
  cat(dimR[1], "x", dimR[2], "\n")
  cat("\nDistinct elements in off-diagonal:")
  if(methods::is(x$R, "diagonalMatrix")) {
    cat("\nR is diagonal\n")
  } else {
    sumR <- Matrix::summary(x$R)
    off_vals <- sumR[sumR$i < sumR$j, "x"]
    if(length(off_vals) == 0) {
      cat("\nNo non-zero off-diagonal elements\n")
    } else {
      print(table(off_vals))
    }
  }
  invisible(x)
}

#' Print function for mc objects.
#' @description Print function for mc objects.
#' @export
#' @param x An object of type mc.
#' @param ... Not used
print.mc <- function(x, ...) {

  cat("Mean model:\n")
  print(x$form)

  if(is.null(x$X)) {
    return(invisible(x))
  }
  dimX <- dim(x$X)
  cat("\nDimension of X:\n")
  cat(dimX[1], "x", dimX[2], "\n")
  invisible(x)
}

#' Print function for summary.svcm
#' @description Print function for summary.svcm objects.
#' @export
#' @param x An object of type summary.svcm.
#' @param ... Not used
print.summary.svcm <- function(x, ...) {
  print(format(c("Log likelihood" = x$logl,
                 "Deviance" = x$dev,
                 "AIC" = x$AIC,
                 "BIC" = x$BIC)), quote = FALSE)
  cat("\nFitted parameters:\n")
  stats::printCoefmat(x$est, digits = 3, signif.stars = FALSE, eps.Pvalue = 0.001)
  invisible(x)
}
