#' Computes functions of parameters
#' @description Evaluates expressions containing functions of model parameters in objects of class \code{svcm}.
#' @export
#' @param mod an instance of an \code{svcm} model.
#' @param form An expression to be evaluated.
compute <- function(mod, form) {
  res <- eval(substitute(form), envir = mod$env_comp)
  return(res)
}

#' Summary function for svcm models.
#' @description Summary function for class \code{"svcm"}.
#' @export
#' @param object An object of type \code{svcm}.
#' @param ... Not used.
summary.svcm <- function(object, ...) {

  if(is.null(object$opt)) {
    stop("The model has not been fitted.")
  }
  theta <- object$opt$par
  ses <- rep(NA, length(theta))
  if(!is.null(object$H)) {
    ses <- sqrt(diag(vcov(object)))
  }
  ll <- logLik(object)
  z <- theta / ses
  est <- data.frame(Estimate = theta,
                    `Std. Error` = ses,
                    `z value` = z,
                    `Pr(>|z|)` =  pchisq(z^2, 1, lower.tail = FALSE),
                    row.names = names(theta), check.names = FALSE)

  ret <- structure(list(N = attr(ll, "nobs"),
                        K = length(theta),
                        logl = ll,
                        dev = object$opt$objective,
                        AIC = AIC(object),
                        BIC = BIC(object),
                        est = est),
                   class = "summary.svcm")
  return(ret)
}

# Now also AIC and BIC should work.
#' Computes log-likelihood for fitm objects.
#' @description Computes log-likelihood for fitm objects.
#' @export
#' @param object An object of type fitm.
#' @param ... Not used.
logLik.svcm <- function(object, ...) {
  ll <- -0.5 * object$opt$objective
  #attr(ll, "nobs") <- length(object$dat$y[object$dat$keepy])
  attr(ll, "nobs") <- length(object$dat$y)
  attr(ll, "df") <- length(object$opt$par)
  class(ll) <- "logLik"
  return(ll)
}

#' Returns fitted model parameters.
#' @description Returns fitted model parameters for \code{svcm} objects.
#' @export
#' @param object An object of type \code{svcm}.
#' @param ... Not used.
#' @return Vector of fitted parameters.
coef.svcm <- function(object, ...) {
  return(object$opt$par)
}

#' Returns covariance matrix of fitted model parameters.
#' @description Returns covariance matrix of fitted model parameters for \code{svcm} objects.
#' @export
#' @param object An object of type \code{svcm}.
#' @param ... Not used.
#' @return Covariance matrix of fitted parameters.
vcov.svcm <- function(object, ...) {
  if(is.null(object$H)) {
    stop("Parameter covariance matrix is not available. Try fitting model with se = TRUE.")
  } else {
    return(solve(0.5 * object$H))
  }
}

#' Deviance tables.
#' @description Produce deviance tables for \code{svcm} objects.
#' @export
#' @param object An object of type \code{svcm}.
#' @param ... list of \code{svcm} models.
#' @return Anova.
anova.svcm <- function(object, ...) {
  mNms <- sapply(sys.call()[-1], as.character)
  dots <- list(...)
  cFits <- sapply(dots, inherits, "svcm")
  if(any(cFits)) {
    fits <- c(list(object), dots[cFits])
  } else {
    fits <- list(object)
  }
# C
  logLiks <- lapply(fits, logLik)
  #nobs <- sapply(logLiks, attr, "nobs")
  logLiksu <- unlist(logLiks)
  dfs <- sapply(logLiks, attr, "df")

  ord <- order(dfs)
  logLiksu <- logLiksu[ord]
  dfs <- dfs[ord]
  fits <- fits[ord]
  mNms <- mNms[ord]

  chisq <- 2 * c(NA, diff(logLiksu))
  dfChisq <- c(NA, diff(dfs))
  pval <- stats::pchisq(chisq, dfChisq, lower.tail = F)

  vals <- data.frame(Df = dfs,
                     AIC = sapply(fits, AIC),
                     BIC = sapply(fits, BIC),
                     logLik = logLiksu,
                     deviance = -2 * logLiksu,
                     Chisq = chisq,
                     `Chi Df` = dfChisq,
                     `Pr(>Chisq)` = pval,
                     row.names = mNms, check.names = FALSE)
  res <- structure(vals, class = c("anova", class(vals)))
  return(res)
}
