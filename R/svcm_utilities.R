#' Evaluate an Expression in an svcm Model Environment
#' @description S3 method for \code{with()} to evaluate expressions in the model's
#' computing environment.
#' @export
#' @method with svcm
#' @param data an instance of an \code{svcm} model.
#' @param expr An expression to be evaluated.
#' @param ... Not used.
with.svcm <- function(data, expr, ...) {
  if(!inherits(data, "svcm")) {
    stop("data must be an object of class svcm.")
  }
  if(is.null(data$env_comp)) {
    stop("Model has no computing environment (env_comp).")
  }
  eval(substitute(expr), envir = data$env_comp)
}

.assert_fitted_svcm <- function(object, fun_name) {
  if(!inherits(object, "svcm")) {
    stop(fun_name, " only supports objects of class svcm.")
  }
  if(is.null(object$opt)) {
    stop("The model has not been fitted.")
  }
}

#' Summary function for svcm models.
#' @description Summary function for class \code{"svcm"}.
#' @export
#' @param object An object of type \code{svcm}.
#' @param ... Not used.
summary.svcm <- function(object, ...) {
  .assert_fitted_svcm(object, "summary.svcm")
  theta <- object$opt$par
  std_err <- rep(NA_real_, length(theta))
  if(!is.null(object$H)) {
    std_err <- sqrt(diag(vcov(object)))
  }
  ll <- logLik(object)
  z_value <- theta / std_err
  est <- data.frame(Estimate = theta,
                    `Std. Error` = std_err,
                    `z value` = z_value,
                    `Pr(>|z|)` =  stats::pchisq(z_value^2, 1, lower.tail = FALSE),
                    row.names = names(theta), check.names = FALSE)

  structure(list(N = attr(ll, "nobs"),
                 K = length(theta),
                 logl = ll,
                 dev = object$opt$objective,
                 AIC = AIC(object),
                 BIC = BIC(object),
                 est = est),
            class = "summary.svcm")
}

# Now also AIC and BIC should work.
#' Computes log-likelihood for svcm objects.
#' @description Computes log-likelihood for svcm objects.
#' @export
#' @param object An object of type svcm.
#' @param ... Not used.
logLik.svcm <- function(object, ...) {
  .assert_fitted_svcm(object, "logLik.svcm")
  ll <- -0.5 * object$opt$objective
  attr(ll, "nobs") <- length(object$dat$y)
  attr(ll, "df") <- length(object$opt$par)
  class(ll) <- "logLik"
  ll
}

#' Returns fitted model parameters.
#' @description Returns fitted model parameters for \code{svcm} objects.
#' @export
#' @param object An object of type \code{svcm}.
#' @param ... Not used.
#' @return Vector of fitted parameters.
coef.svcm <- function(object, ...) {
  .assert_fitted_svcm(object, "coef.svcm")
  object$opt$par
}

#' Returns covariance matrix of fitted model parameters.
#' @description Returns covariance matrix of fitted model parameters for \code{svcm} objects.
#' @export
#' @param object An object of type \code{svcm}.
#' @param ... Not used.
#' @return Covariance matrix of fitted parameters.
vcov.svcm <- function(object, ...) {
  .assert_fitted_svcm(object, "vcov.svcm")
  if(is.null(object$H)) {
    stop("Parameter covariance matrix is not available. Try fitting model with se = TRUE.")
  }
  solve(0.5 * object$H)
}

#' Deviance tables.
#' @description Produce deviance tables for \code{svcm} objects.
#' @export
#' @param object An object of type \code{svcm}.
#' @param ... list of \code{svcm} models.
#' @return Anova.
anova.svcm <- function(object, ...) {
  .assert_fitted_svcm(object, "anova.svcm")
  dots <- list(...)
  if(length(dots) > 0 && !all(vapply(dots, inherits, logical(1), "svcm"))) {
    stop("All models in ... must be of class svcm.")
  }
  if(length(dots) > 0 && !all(vapply(dots, function(x) !is.null(x$opt), logical(1)))) {
    stop("All models in ... must be fitted before calling anova().")
  }

  fits <- c(list(object), dots)

  call_expr <- match.call(expand.dots = FALSE)
  dot_calls <- call_expr$...
  dot_names <- if(length(dot_calls) == 0) {
    character(0)
  } else {
    vapply(dot_calls, deparse1, character(1))
  }
  model_names <- c(deparse1(call_expr$object), dot_names)

  log_liks <- lapply(fits, logLik)
  loglik_values <- unlist(log_liks)
  dfs <- vapply(log_liks, attr, numeric(1), "df")
  n_obs <- vapply(log_liks, attr, numeric(1), "nobs")
  if(length(unique(n_obs)) > 1) {
    warning("Models were fit to different numbers of observations.")
  }

  order_idx <- order(dfs)
  fits <- fits[order_idx]
  model_names <- model_names[order_idx]
  loglik_values <- loglik_values[order_idx]
  dfs <- dfs[order_idx]

  chisq <- 2 * c(NA_real_, diff(loglik_values))
  df_chisq <- c(NA_real_, diff(dfs))
  p_values <- rep(NA_real_, length(chisq))
  keep <- !is.na(chisq) & !is.na(df_chisq) & df_chisq > 0
  p_values[keep] <- stats::pchisq(chisq[keep], df_chisq[keep], lower.tail = FALSE)

  anova_table <- data.frame(Df = dfs,
                     AIC = vapply(fits, AIC, numeric(1)),
                     BIC = vapply(fits, BIC, numeric(1)),
                     logLik = loglik_values,
                     deviance = -2 * loglik_values,
                     Chisq = chisq,
                     `Chi Df` = df_chisq,
                     `Pr(>Chisq)` = p_values,
                     row.names = model_names, check.names = FALSE)
  structure(anova_table, class = c("anova", class(anova_table)))
}
