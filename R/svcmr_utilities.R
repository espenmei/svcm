# #' Compute hessian
# #' @description Computes hessian for fitted models of type \code{fitm}.
# #' @export
# #' @param object A fitted model of type \code{fitm}.
#SEm <- function(fitm) {
#  H <- numDeriv::hessian(fitm$svcm$objective, fitm$fit$par)
#  dimnames(H) <- list(names(fitm$fit$par), names(fitm$fit$par))
#  fitm$hessian <- H
#  return(fitm)
#}

#' Computes functions of parameters
#' @description \code{compute} is a generic function for evaluating expressions containing functions of parameters in models.
#' @export
#' @param object Either a model of type \code{svcm} or a fitted model of type \code{fitm}.
#' @param ... Further function arguments.
compute <- function(object, ...) {
  UseMethod("compute")
}

#' Computes functions of parameters
#' @description Evaluates expressions containing functions of model parameters in objects of class \code{svcm}.
#' @export
#' @param object An object of class \code{svcm}.
#' @param form An expression to be evaluated.
#' @param ... Not used.
compute.svcm <- function(object, form, ...) {
  en <- lapply(object$mos, "[[", "values")
  names(en) <- lapply(object$mos, "[[", "name")
  res <- eval(substitute(form), envir = en)
  return(res)
}

#' Computes functions of parameters
#' @description Evaluates expressions containing functions of model parameters in objects of class \code{fitm}.
#' @export
#' @param object An object of class \code{fitm}.
#' @param form An expression to be evaluated.
#' @param ... Not used.
compute.fitm <- function(object, form, ...) {
  en <- lapply(object$svcm$mos, "[[", "values")
  names(en) <- lapply(object$svcm$mos, "[[", "name")
  res <- eval(substitute(form), envir = en)
  return(res)
}

#' Summary function for svcmr models.
#' @description Summary function for class \code{"fitm"}.
#' @export
#' @param object An object of type \code{fitm}.
#' @param ... Not used.
summary.fitm <- function(object, ...) {
  theta <- object$fit$par
  ses <- rep(NA, length(theta))
  if(!is.null(object$hessian)) {
    ses <- sqrt(diag(vcov(object)))
  }
  zVal <- theta / ses
  ll <- logLik(object)
  ret <- structure(list(N = attr(ll, "nobs"),
                        K = length(theta),
                        logl = ll,
                        dev = object$fit$objective,
                        AIC = AIC(object),
                        BIC = BIC(object),
                        est = data.frame(Estimate = theta,
                                         `Std. Error` = ses,
                                         `Z value` = zVal,
                                         `Pr(>|z|)` =  pchisq(zVal^2, 1, lower.tail = FALSE),
                                         row.names = names(theta), check.names = FALSE)),
                   class = "summary.fitm")
  return(ret)
}

# Now also AIC and BIC should work.
#' Computes log-likelihood for fitm objects.
#' @description Computes log-likelihood for fitm objects.
#' @export
#' @param object An object of type fitm.
#' @param ... Not used.
logLik.fitm <- function(object, ...) {
  ll <- -0.5 * object$fit$objective
  attr(ll, "nobs") <- length(object$y)
  attr(ll, "df") <- length(object$fit$par)
  class(ll) <- "logLik"
  return(ll)
}

#' Returns fitted model parameters.
#' @description Returns fitted model parameters for fitm objects.
#' @export
#' @param object An object of type fitm.
#' @param ... Not used.
#' @return Fitted parameters.
coef.fitm <- function(object, ...) {
  return(object$fit$par)
}

#' Returns fitted model parameters.
#' @description Returns detailed table of model parameters for objects of class \code{summary.fitm}.
#' @export
#' @param object An object of class summary.fitm.
#' @param ... Not used.
#' @return \code{data.frame} of model parameters with details.
coef.summary.fitm <- function(object, ...) {
  return(object$est)
}

#' Returns covariance matrix of fitted model parameters.
#' @description Returns covariance matrix of fitted model parameters for fitm objects.
#' @export
#' @param object An object of type fitm.
#' @param ... Not used.
#' @return Covariance matrix of fitted parameters.
vcov.fitm <- function(object, ...) {
  if(is.null(object$hessian)) {
    stop("Parameter covariance matrix is not available. Try fitting model with se = TRUE.")
  } else {
    return(solve(0.5 * object$hessian))
  }
}

#' Deviance tables.
#' @description Produce deviance tables for fitm objects.
#' @export
#' @param object An object of type fitm.
#' @param ... Not used.
#' @return Anova.
anova.fitm <- function(object, ...) {
  mNms <- sapply(sys.call()[-1], as.character)
  dots <- list(...)
  cFits <- sapply(dots, inherits, "fitm")
  if(any(cFits)) {
    fits <- c(list(object), dots[cFits])
  } else {
    fits <- list(object)
  }

  logLiks <- lapply(fits, logLik)
  nobs <- sapply(logLiks, attr, "nobs")
  logLiksu <- unlist(logLiks)
  #Dfs <- vapply(logLiks, attr, FUN.VALUE = numeric(1), "df")
  Dfs <- sapply(logLiks, attr, "df")
  chisq <- 2 * c(NA, logLiksu[1] - logLiksu[-1])
  dfChisq <- c(NA, Dfs[1] - Dfs[-1])
  pval <- stats::pchisq(chisq, dfChisq, lower.tail = F)
  vals <- data.frame(Df = Dfs,
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
