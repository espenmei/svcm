# Relies on these
#library(Matrix)
#library(sparseMVN)

# Used during optimization for updating model objects in computing environment.
# Assigns values according to equal labels.
.UpdateValues <- function(mo, values, env_comp) {
  ind_mo <- which(mo$labels %in% names(values))
  ind_values <- charmatch(mo$labels, names(values), nomatch = 0)
  mo$values[ind_mo] <- values[ind_values[ind_values > 0]]
  assign(mo$name, mo$values, envir = env_comp)
}

.GetFreeValues <- function(mo) {
  return(mo$values[mo$free])
}

.GetFreeLabels <- function(mo) {
  return(mo$labels[mo$free])
}

#' Constructor for parameter matrix
#' @description Creates an object of type pm used to represent parameters in a model.
#' @export
#' @param nrow Number of matrix rows.
#' @param ncol Number of matrix columns.
#' @param labels Vector or matrix of characters giving labels to parameters. Equal characters can be used to define equality constraints.
#' @param values Vector or matrix with parameter values. Values for free parameters are used as starting values during optimization
#' and values for fixed parameters are constants during optimization.
#' @param free Vector or matrix of logical values defining whether parameters are free or not.
#' @param name Character giving name to model object.
#' @return An object of class mo.
#' @examples
#' library(svcmr)
#' S <- pm(nrow = 2, ncol = 2,
#'         labels = c("s11", "s12", "s12", "s22"),
#'         values = c(2, 1, 1, 2),
#'         free = TRUE,
#'         name = "S")
pm <- function(nrow = NA,
               ncol = NA,
               labels = character(),
               free = FALSE,
               values = NA,
               name = NA) {
  free_mat <- matrix(free, nrow, ncol)
  # If values is matrix then:
  # Warning message:
  #  In Matrix(values, nrow, ncol, sparse = FALSE, doDiag = FALSE) :
  #  'nrow', 'ncol', etc, are disregarded for matrix 'data'
  # Fix?
  values_mat <- Matrix::Matrix(values, nrow, ncol, sparse = FALSE, doDiag = FALSE)
  labels_mat <- matrix(labels, nrow, ncol)
  if (is.na(name)) {
    stop("All objects of type mo must have a name.")
  }
  ret <- structure(list(free = free_mat,
                        values = values_mat,
                        labels = labels_mat,
                        name = name),
                   class = "pm")
  return(ret)
}

#' Constructor for variance component object
#' @description Creates a variance component object as a function of model objects and a relationship matrix.
#' @export
#' @param form An expression for the covariance model.
#' @param R A relationship matrix. Should be a sparse matrix from the Matrix package.
#' @return An object of class svc.
#' @examples
#' library(svcmr)
#' R <- Matrix::Diagonal(100)
#' L <- pm(nrow = 6, ncol = 2,
#'         labels = paste0("l", 1:12),
#'         values = c(1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1),
#'         free = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE,
#'                  FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
#'         name = "L")
#' S <- pm(nrow = 2, ncol = 2,
#'         labels = c("s11", "s12", "s12", "s22"),
#'         values = c(2, 1, 1, 2),
#'         free = TRUE,
#'         name = "S")
#' vc <- svc(form = L %*% S %*% t(L), R = R)
svc <- function(form, R = NULL) {
  if (is.null(R)) {
    stop("A relationship matrix must be supplied.")
  }
  if(!Matrix::isSymmetric(R)) {
    stop("Relationship matrix must be symmetric.")
  }
  #expr <- enexpr(form) # Turn the covariance formula into an expression. quote?
  expr <- substitute(form)
  ret <- structure(list(form = expr,
                        R = R),
                   class = "svc")
  return(ret)
}

#' Constructor for mean component object
#' @description Creates a mean component object used to define (part of) the mean structure in the model.
#' @export
#' @param form An expression for the mean model.
#' @param X A design matrix.
#' @return An object of class mc.
#' @examples
#' library(svcmr)
#' X <- cbind(1, rnorm(100))
#' B <- pm(nrow = 6, ncol = 2,
#'         labels = paste0("b", 1:12),
#'         values = 0,
#'         free = TRUE,
#'         name = "B")
#' mc <- mc(form = B, X = X)
mc <- function(form, X = NULL) {
  if (is.null(X)) {
    stop("A design matrix must be supplied.")
  }
  if (anyNA(X)) {
    stop("Missing values in X is not allowed.")
  }
  if (!is.numeric(X)) {
    stop("X must be numeric.")
  }
  expr <- substitute(form)
  ret <- structure(list(form = expr,
                        X = X),
                   class = "mc")
  return(ret)
}

#' Generic compute function.
#' @description Internal functions used to evaluate expressions containing functions of model objects.
#' @export
#' @importFrom Matrix t
#' @param object An object of type mc/svc.
#' @param ... Further function arguments.
.computeC <- function(object, ...) {
  UseMethod(".computeC")
}

#' Compute function for mc object.
#' @description Internal function used to evaluate expressions containing functions of model objects.
#' @export
#' @importFrom Matrix t
#' @param object An object of type mc.
#' @param env Computing environment.
#' @param ... Not used.
.computeC.mc <- function(object, env, ...) {
  # Use Matrix::t(object$X)) ?
  mci <- eval(object$form, env)
  mu <- as.vector(Matrix::t(mci %*% t(object$X)))
  return(mu)
}

#' Compute function for svc object.
#' @description Internal function used to evaluate expressions containing functions of model objects.
#' @export
#' @importFrom Matrix t
#' @param object An object of type svc.
#' @param env Computing environment.
#' @param ... Not used.
.computeC.svc <- function(object, env, ...) {
  vci <- eval(object$form, env)
  sigma <- vci %x% object$R
  return(sigma)
}

#' Creates a model
#' @description Creates a new model from model, mean and variance components objects
#' @export
#' @param Y Data described by model.
#' @param ... All relevant model, mean and variance components objects used to define the model.
#' @return An object of type svcm.
#' @examples
#' library(svcmr)
#' R <- Matrix::Diagonal(100)
#' L <- pm(nrow = 6, ncol = 2,
#'         labels = paste0("l", 1:12),
#'         values = c(1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1),
#'         free = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE,
#'                  FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
#'         name = "L")
#' S <- pm(nrow = 2, ncol = 2,
#'         labels = c("s11", "s12", "s12", "s22"),
#'         values = c(2, 1, 1, 2),
#'         free = TRUE,
#'         name = "S")
#' vc <- svc(form = L %*% S %*% t(L), R = R)
#' X <- cbind(1, rnorm(100))
#' B <- pm(nrow = 6, ncol = 2,
#'         labels = paste0("b", 1:12),
#'         values = 0,
#'         free = TRUE,
#'         name = "B")
#' mc1 <- mc(form = B, X = X)
#' Y <- matrix(rnorm(100 * 4), 100, 4)
#' mod <- svcm(Y, L, S, vc, B, mc1)
svcm <- function(Y, ...) {
  if (anyNA(Y)) {
    stop("Missing values in Y is not supported.")
  }
  if (!is.numeric(Y)) {
    stop("Y must be numeric.")
  }
  # Stack Y - the order is always var1[1], var1[2], ..., var1[n], var2[1], var2[2], ..., var2[n]
  y = c(Y)
  # Extract only objects of type mo, svc or mc and ignore anything else.
  input <- list(...)
  mos <- input[sapply(input, inherits, "pm")]
  svcs <- input[sapply(input, inherits, "svc")]
  mcs <- input[sapply(input, inherits, "mc")]
  if (length(mos) == 0 || length(svcs) == 0 || length(mcs) == 0) {
    stop("At least one pm, svc and mc object must be supplied.")
  }
  # Objective function
  # Note that objective is now working on mvs and svcs as they were given to
  # model object and not fit object. Weird? No, should be good!
  objective <- function(env_comp) {
    # ?
    M <- Reduce("+", lapply(mcs, .computeC, env_comp))
    S <- Reduce("+", lapply(svcs, .computeC, env_comp))
    lS <- Matrix::Cholesky(S)
    ll <- sparseMVN::dmvn.sparse(y, M, CH = lS, prec = FALSE)
    return(-2 * ll)
  }
  ret <- structure(list(y = y,
                        mos = mos,
                        svcs = svcs,
                        mcs = mcs,
                        objective = objective),
                   class = "svcm")
  return(ret)
}

#' Fit a model
#' @description Fits a model returned from svcm.
#' @export
#' @param svcm An object of type svcm.
#' @param se Should standard errors be computed?
#' @param ... Arguments passed to nlminb.
#' @return An object of class fitm.
#' @examples
#' library(svcmr)
#' R <- Matrix::Diagonal(100)
#' L <- pm(nrow = 6, ncol = 2,
#'         labels = paste0("l", 1:12),
#'         values = c(1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1),
#'         free = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE,
#'                  FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
#'         name = "L")
#' S <- pm(nrow = 2, ncol = 2,
#'         labels = c("s11", "s12", "s12", "s22"),
#'         values = c(2, 1, 1, 2),
#'         free = TRUE,
#'         name = "S")
#' vc <- svc(form = L %*% S %*% t(L), R = R)
#' X <- cbind(1, rnorm(100))
#' B <- pm(nrow = 6, ncol = 2,
#'         labels = paste0("b", 1:12),
#'         values = 0,
#'         free = TRUE,
#'         name = "B")
#' mc1 <- mc(form = B, X = X)
#' Y <- matrix(rnorm(100 * 4), 100, 4)
#' mod <- svcm(Y, L, S, vc, B, mc1)
#' \dontrun{
#' fit <- fitm(mod, se = TRUE)
#' }
fitm <- function(svcm, se = FALSE, ...) {
  if (!inherits(svcm, "svcm")) {
    stop("Only objects of type svcm are accepted.")
  }
  # Set start values
  theta_start <- unlist(lapply(svcm$mos, .GetFreeValues))
  names(theta_start) <- unlist(lapply(svcm$mos, .GetFreeLabels))
  # Keep only one (first) when equal labels (equality constraints)
  theta_start_u <- theta_start[!duplicated(names(theta_start))]
  # Create computing environment for fitting
  env_comp <- new.env()
  # Fitting function
  fit_objective <- function(theta) {
    # Update computing environment
    lapply(svcm$mos, .UpdateValues, theta, env_comp)
    # Compute objective
    return(svcm$objective(env_comp))
  }
  # optimize model
  cat("\niteration: objective:", names(theta_start_u), "\n")
  time_start <- proc.time()
  fit <- nlminb(theta_start_u, fit_objective, ...)
  time_used <- proc.time() - time_start
  if(fit$convergence != 0) {
    warning("Optimization may not have converged.",
            " \nnlminb convergence code: ", fit$convergence,
            " \nnlminb message: ", fit$message)
  }

  # Hessian at minimum
  H <- NULL
  if (se) {
    message("Computing standard errors.")
    #H <- optimHess(res$par, fit_objective)
    H <- numDeriv::hessian(fit_objective, fit$par)
  }
  # Update svcm object with values from solution before return.
  for (i in seq_along(svcm$mos)) {
    svcm$mos[[i]]$values <- get(svcm$mos[[i]]$name, envir = env_comp)
  }
  ret <- structure(list(fit = fit,
                        hessian = H,
                        time = time_used,
                        svcm = svcm),
                   class = "fitm")
  return(ret)
}
