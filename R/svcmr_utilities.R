#' Computes functions of parameters.
#' @description Can be used to evaluate expressions containing functions of parameters.
#' @export
#' @param object Either a model og type svcm or a fitted model of type fitm.
#' @param ... Further function arguments.
compute <- function(object, ...) {
  UseMethod("compute")
}

#' Computes functions of parameters in model objects.
#' @description Can be used to evaluate expressions containing functions of model objects.
#' @export
#' @param object A model of type svcm.
#' @param form An expression to be evaluated.
#' @param ... Not used.
compute.svcm <- function(object, form, ...) {
  en <- lapply(object$mos, "[[", "values")
  names(en) <- lapply(object$mos, "[[", "name")
  res <- eval(substitute(form), envir = en)
  return(res)
}

#' Computes functions of parameters in fitted model objects.
#' @description Can be used to evaluate expressions containing functions of parameters.
#' @export
#' @param object A Fitted model of type fitm.
#' @param form An expression to be evaluated.
#' @param ... Not used.
compute.fitm <- function(object, form, ...) {
  en <- lapply(object$svcm$mos, "[[", "values")
  names(en) <- lapply(object$svcm$mos, "[[", "name")
  res <- eval(substitute(form), envir = en)
  return(res)
}

# Summary functions
# -----------------------------------------
#' Summary function for fitm objects.
#' @description Summary function for fitm objects.
#' @export
#' @param object An object of type fitm.
#' @param ... Not used.
summary.fitm <- function(object, ...) {
  theta <- object$fit$par
  ses <- rep(NA, length(theta))
  if(!is.null(object$hessian)) {
    ses <- sqrt(diag(vcov(object)))
  }
  zVal <- theta / ses
  ret <- structure(list(N = nobs(object),
                        K = length(theta),
                        logl = logLik(object),
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

# #' Generic logLik.
# #' @export
# #' @param object An object of type fitm.
# #' @param ... Not used.
# logLik <- function(object, ...) {
#   UseMethod("logLik")
# }
# Now also noba, AIC and BIC should work.
#' Computes log-likelihood for fitm objects.
#' @description Computes log-likelihood for fitm objects.
#' @export
#' @param object An object of type fitm.
#' @param ... Not used.
logLik.fitm <- function(object, ...) {
  ll <- -0.5 * object$fit$objective
  attr(ll, "nobs") <- length(object$svcm$y)
  attr(ll, "df") <- length(object$fit$par)
  class(ll) <- "logLik"
  return(ll)
}

# # Not needed?
# #' Generic nobs.
# #' @export
# #' @param object fitm object.
# #' @param ... Not used.
# nobs <- function(object, ...) {
#   UseMethod("nobs")
# }
# #' Returns number of observations used to fit model.
# #' @description Returns number of observations used to fit model.
# #' @export
# #' @param object An object of type fitm.
# #' @param ... Not used.
# #' @return Number of observations.
# nobs.fitm <- function(object, ...) {
#   n <- length(object$svcm$y)
#   return(n)
# }

# Shouldn't be nencessay to implement as loglik and nobs are imp?
# #' Returns AIC for fitted model.
# #' @description Returns AIC for fitted model.
# #' @export
# #' @param object An object of type fitm.
# #' @param k k.
# #' @return AIC.
#AIC.fitm <- function(object, k = 2) {
#  lls <- logLik(object)
#  val <- -2 * as.numeric(lls) + k * attr(lls, "df")
#  return(val)
#}

# #' Returns BIC for fitted model.
# #' @description Returns BIC for fitted model.
# #' @export
# #' @param object An object of type fitm.
# #' @return BIC.
#BIC.fitm <- function(object) {
#  lls <- logLik(object)
#  nos <- attr(lls, "nobs")
#  val <- -2 * as.numeric(lls) + log(nos) * attr(lls, "df")
#  return(val)
#}

#' Returns fitted model parameters.
#' @description Returns fitted model parameters for fitm objects.
#' @export
#' @param object An object of type fitm.
#' @param ... Not used.
#' @return Fitted parameters.
coef.fitm <- function(object, ...) {
  return(object$fit$par)
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

#' Anova.
#' @description Anova.
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
  nobs <- sapply(fits, nobs)
  logLiks <- lapply(fits, logLik)
  logLiksu <- unlist(logLiks)
  Dfs <- vapply(logLiks, attr, FUN.VALUE = numeric(1), "df")
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
