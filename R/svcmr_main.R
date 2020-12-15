#' Updates model objects
#' @description Used during optimization for updating model objects in computing environment.
#' @export
#' @param mo model object.
#' @param values free parameters.
#' @param env_comp computing environment.
.updateValues <- function(mo, values, env_comp) {
  pos_free_mo <- which(mo$free) # positions of free values in mo
  ind_values <- pmatch(mo$labels[pos_free_mo], names(values), duplicates.ok = TRUE) # Positions of pos_free_mo in values
  mo$values[pos_free_mo] <- values[ind_values] # Update mo
  assign(mo$name, mo$values, envir = env_comp) # update environment
}

.getFreeValues <- function(mo) {
  return(mo$values[mo$free])
}

.getFreeLabels <- function(mo) {
  return(mo$labels[mo$free])
}

#' Constructor for parameter matrix
#' @description Creates an object of type \code{pm} used to represent parameters in a model.
#' @export
#' @param nrow Number of matrix rows.
#' @param ncol Number of matrix columns.
#' @param labels Vector/matrix with parameter labels. Equal characters defines equality constraints.
#' @param values Vector/matrix with parameter values. Values for free parameters are used as starting values for optimization
#' and values for fixed parameters are constants during optimization.
#' @param free Vector or matrix of logical values defining whether parameters are free or not.
#' @param name Character giving name to model object.
#' @return An object of class \code{pm}.
#' @examples
#' library(svcmr)
#' S <- pm(nrow = 2, ncol = 2,
#'         labels = c("s11", "s12", "s12", "s22"),
#'         values = c(2, 1, 1, 2),
#'         free = TRUE,
#'         name = "S")
pm <- function(nrow = NA, ncol = NA, labels = character(),
               free = FALSE, values = NA, name = NA) {

  free_mat <- matrix(free, nrow, ncol)
  values_mat <- Matrix::Matrix(c(values), nrow, ncol, sparse = FALSE, doDiag = FALSE)
  labels_mat <- matrix(labels, nrow, ncol)
  if (is.na(name)) {
    stop("A name is required.")
  }
  ret <- structure(list(free = free_mat,
                        values = values_mat,
                        labels = labels_mat,
                        name = name),
                   class = "pm")
  return(ret)
}

#' Constructor for intermediate calculation step object
#' @description Creates an object of type \code{ic} used to represent intermediate calculation steps in a model.
#' @export
#' @param form An expression for the calculation.
#' @param name Character giving name to object.
#' @return An object of class \code{ic}.
ic <- function(form, name = character()) {
  if (length(name) == 0) {
    stop("A name is required.")
  }
  ret <- structure(list(form = substitute(form),
                        name = name),
                   class = c("free", "ic"))
  return(ret)
}

#' Constructor for variance component object
#' @description Creates a variance component object as a function of model objects and (optionally) a relationship matrix.
#' @export
#' @param form An expression for the covariance model.
#' @param R A relationship matrix. Should be a sparse matrix from the \code{Matrix} package.
#' @return An object of class \code{svc}.
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
  if(is.null(R)) {
    ret <- .new_freesvc(form = substitute(form))
  } else {
    stopifnot(inherits(R, "Matrix"), Matrix::isSymmetric(R))
    ret <- .new_fixedsvc(form = substitute(form), R = R)
  }
  return(ret)
}
.new_svc <- function(..., class = character()) {
  ret <- structure(list(...),
                   class = c(class, "svc"))
  return(ret)
}
.new_fixedsvc <- function(form, R) {
  .new_svc(form = form, R = R, class = "fixedsvc")
}
.new_freesvc <- function(form) {
  .new_svc(form = form, class = "free")
}

#' Constructor for mean component object
#' @description Creates a mean component object used to define (part of) the mean structure in the model.
#' @export
#' @param form An expression for the mean model.
#' @param X A design matrix.
#' @return An object of class \code{mc}.
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
  if(is.null(X)) {
    ret <- .new_freemc(form = substitute(form))
  } else {
    stopifnot(!anyNA(X), is.numeric(X))
    ret <- .new_fixedmc(form = substitute(form), X = X)
  }
  return(ret)
}
.new_mc <- function(..., class = character()) {
  ret <- structure(list(...),
                   class = c(class, "mc"))
  return(ret)
}
.new_fixedmc <- function(form, X) {
  .new_mc(form = form, X = X, Xt = t(X), class = "fixedmc")
}
.new_freemc <- function(form) {
  .new_mc(form = form, class = "free")
}

#' Generic compute function
#' @description Internal functions used to evaluate expressions containing functions of model objects.
#' @export
#' @importFrom Matrix t
#' @param object An object of type mc/svc.
#' @param ... Further function arguments.
.computeC <- function(object, ...) {
  UseMethod(".computeC")
}

#' Compute function for mc object
#' @description Internal function used to evaluate expressions containing objects of type \code{fixedmc}.
#' @export
#' @importFrom Matrix t
#' @param object An object of type mc.
#' @param env Computing environment.
#' @param ... Not used.
.computeC.fixedmc <- function(object, env, ...) {
  mci <- eval(object$form, env)
  mu <- as.vector(Matrix::t(mci %*% object$Xt)) # Vec(t(B %*% t(X)))
  return(mu)
}

#' Compute function for svc objects
#' @description Internal function used to evaluate expressions containing objects of type \code{fixedsvc}.
#' @export
#' @param object An object of type svc.
#' @param env Computing environment.
#' @param ... Not used.
.computeC.fixedsvc <- function(object, env, ...) {
  vci <- eval(object$form, env)
  sigma <- vci %x% object$R
  return(sigma)
}

#' Compute function for svc object.
#' @description Internal function used to evaluate expressions containing objects of type \code{free}.
#' @export
#' @param object An object of type \code{free}.
#' @param env Computing environment.
#' @param ... Not used.
.computeC.free <- function(object, env, ...) {
  res <- eval(object$form, env)
  return(res)
}

#' Creates a model
#' @description Creates a new model
#' @export
#' @param Y Matrix of data described by model.
#' @param ... All relevant model, mean and variance components objects used to define the model.
#' @return An object of type svcm.
svcm <- function(Y, ...) {
  # Extract only objects of type pm, svc, mc and ic and ignore anything else.
  dots <- list(...)
  pms <- dots[sapply(dots, inherits, "pm")]
  svcs <-dots[sapply(dots, inherits, "svc")]
  mcs <- dots[sapply(dots, inherits, "mc")]
  ics <- dots[sapply(dots, inherits, "ic")]
  if (length(pms) == 0 || length(svcs) == 0 || length(mcs) == 0) {
    stop("At least one pm, svc and mc object must be supplied.")
  }

  ret <- structure(list(dat = datm(Y),
                        pms = pms,
                        svcs = svcs,
                        mcs = mcs,
                        ics = ics,
                        env_comp = new.env(), # Computing environment
                        opt = NULL,
                        H = NULL),
                   class = "svcm")
  return(ret)
}

#' Create a object for storing data
#' @description creates a \code{datm} object.
#' @export
#' @param Y data described by model.
#' @return an object of class \code{datm}.
datm <- function(Y) {
  if (!is.numeric(Y)) {
    stop("Y must be numeric.")
  }
  # Stack Y - the order is always var1[1], var1[2], ..., var1[n], var2[1], var2[2], ..., var2[n]
  y <- c(Y)
  # Find positions of non-missing values
  keepy <- !is.na(y)

  structure(list(Y = Y,
                 y = y,
                 keepy = keepy),
            class = "datm")
}

#' Fit a model
#' @description fits a \code{svcm} model.
#' @export
#' @param svcm an object of class \code{svcm}.
#' @param se should standard errors be computed?
#' @param ... arguments passed to \code{nlminb}.
#' @return an object of class \code{svcm}.
fitm <- function(svcm, se = FALSE, ...) {

  if (!inherits(svcm, "svcm")) {
    stop("Only objects of type svcm are accepted.")
  }

  if (!is.null(svcm$opt)) {
    stop("This model has already been fitted.")
  }

  fit_objective <- function(theta) {
    lapply(svcm$pms, .updateValues, theta, svcm$env_comp)
    lapply(svcm$ics, function(x) assign(x$name, .computeC(x, svcm$env_comp), envir = svcm$env_comp))
    return(objective(svcm))
  }

  # Set start values
  theta_start <- unlist(lapply(svcm$pms, .getFreeValues))
  names(theta_start) <- unlist(lapply(svcm$pms, .getFreeLabels))
  # Keep only one (first) when equal labels (equality constraints)
  theta_start_u <- theta_start[!duplicated(names(theta_start))]
  # optimize model
  cat("\niteration: objective:", names(theta_start_u), "\n")
  time_start <- proc.time()
  fit <- nlminb(theta_start_u, fit_objective, ...)
  fit$time <- proc.time() - time_start
  if(fit$convergence != 0) {
    warning("Optimization may not have converged.",
            " \nnlminb convergence code: ", fit$convergence,
            " \nnlminb message: ", fit$message)
  }

  svcm$opt = fit

  # Update pm objects with values from solution before return. Consider updating ic objects as well.
  for (i in seq_along(svcm$pms)) {
    svcm$pms[[i]]$values <- get(svcm$pms[[i]]$name, envir = svcm$env_comp)
  }

  # Hessian at minimum
  if(se) {
    message("Computing standard errors.")
    svcm$H <-compHess(fit_objective, fit$par)
  }

  return(svcm)
}

#' Objective function
#' @description Objective function
#' @export
#' @param mod An object of type \code{fitm}
#' @return Twice negative log likelihood.
objective <- function(mod) {

  M <- expectedM(mod)
  S <- expectedS(mod)
  y <- mod$dat$y[mod$dat$keepy]

  ch <- Matrix::Cholesky(S)
  rm <- y - M
  r2 <- Matrix::solve(ch, rm)
  deter <- 2 * Matrix::determinant(ch)$modulus
  dev <- log(2*pi) * length(y) + deter + sum(rm * r2)
  return(dev)
}

#' Compute hessian
#' @description Computes hessian
#' @export
#' @param fit_objective Function returned from \code{fitm} defining the objective function.
#' @param par Parameter vector.
#' @param ... Arguments passed to \code{numDeriv::hessian}.
#' @return Hessian matrix.
compHess <- function(fit_objective, par, ...) {
  H <- numDeriv::hessian(fit_objective, par, ...)
  dimnames(H) <- list(names(par), names(par))
  return(H)
}

#' Model implied mean vector
#' @description Compute model implied mean vector.
#' @export
#' @param svcm an instance of an \code{svcm} model.
#' @param dropmiss drop expectation for missing values?
#' @return vector of means.
expectedM <- function(svcm, dropmiss = T) {
  M <- Reduce("+", lapply(svcm$mcs, .computeC, svcm$env_comp))
  if(dropmiss) {
    M <- M[svcm$dat$keepy]
  }
  return(M)
}

#' Model implied covariance matrix
#' @description Compute model implied covariance matrix.
#' @export
#' @param svcm an instance of an \code{svcm} model.
#' @param dropmiss drop expectation for missing values?
#' @return matrix of (co)variances.
expectedS <- function(svcm, dropmiss = T) {
  S <- Reduce("+", lapply(svcm$svcs, .computeC, svcm$env_comp))
  if(dropmiss) {
    S <- S[svcm$dat$keepy, svcm$dat$keepy]
  }
  return(S)
}
