#' Generic update function
#.update <- function(object, x, env, ...) {
#  UseMethod(".update")
#}

#' Update function for pm objects
#.update.pm <- function(object, x, env, ...) {
#  pos_free_mo <- which(object$free) # positions of free values in mo
#  ind_x <- pmatch(object$labels[pos_free_mo], names(x), duplicates.ok = TRUE) # Positions of pos_free_mo in values
#  object$values[pos_free_mo] <- x[ind_x] # Update mo
#  assign(object$name, object$values, envir = env)
#}

#' Update function for ic objects
#' Denne er ikke riktig. Feil signatur. trenger ikke theta (x)
#.update.ic <- function(object, x, env, ...) {
#  assign(object$name, .computeC(x, env), envir = env)
#}

#' Return free values in parameter matrix as a vector
#' @param pm parameter matrix object.
.get_free_values <- function(pm) {
  return(pm$values[pm$free])
}

#' Return labels of free values in parameter matrix as a vector
#' @param pm parameter matrix object.
.get_free_labels <- function(pm) {
  return(pm$labels[pm$free])
}

#' Updates values in parameter matrix
#' @description Used during optimization for updating free values in a parameter matrix in both computing environment and model object.
#' @param pm parameter matrix object.
#' @param theta vector of free model parameters.
#' @param env_comp computing environment where the expression is evaluated and stored.
.update_pm <- function(pm, theta, env_comp) {
  pos_free_pm <- which(pm$free) # positions of free values in pm
  ind_theta <- pmatch(pm$labels[pos_free_pm], names(theta), duplicates.ok = TRUE) # Positions of pos_free_pm in values
  pm$values[pos_free_pm] <- theta[ind_theta] # Update pm
  assign(pm$name, pm$values, envir = env_comp) # update environment
}

#' Updates all parameter matrix objects of a model
#' @description Used during optimization for updating parameter matrices in computing environment and model object.
#' @param mod an instance of an \code{svcm} model.
#' @param theta vector of free model parameters.
.update_pms <- function(mod, theta) {
  lapply(mod$pms, .update_pm, theta, mod$env_comp)
}

#' Updates intermediate computations
#' @description Used during optimization for computing expressions in \code{ic} objects.
#' @param pm parameter matrix object.
#' @param env_comp computing environment where the expression is evaluated and stored.
.update_ic <- function(ic, env_comp) {
  assign(ic$name, .compute(ic, env_comp), envir = env_comp)
}

#' Updates all intermediate computation objects of a model
#' @description Used during optimization for updating intermediate computations in computing environment.
#' @param mod an instance of an \code{svcm} model.
.update_ics <- function(mod) {
  lapply(mod$ics, .update_ic, mod$env_comp)
}

#' Generic compute function
#' @description Internal functions used to evaluate expressions containing functions of model objects.
#' @export
#' @importFrom Matrix t
#' @param object An object of type \code{mc} / \code{svc}.
#' @param ... Further function arguments.
.compute <- function(object, ...) {
  UseMethod(".compute")
}

#' Compute function for mc object
#' @description Internal function used to evaluate expressions containing objects of type \code{fixedmc}.
#' @export
#' @importFrom Matrix t
#' @param object An object of type \code{mc}.
#' @param env Computing environment.
#' @param ... Not used.
.compute.fixedmc <- function(object, env, ...) {
  mci <- eval(object$form, env)
  mu <- as.vector(Matrix::t(mci %*% object$Xt)) # Vec(t(B %*% t(X)))
  return(mu)
}

#' Compute function for svc objects
#' @description Internal function used to evaluate expressions containing objects of type \code{fixedsvc}.
#' @export
#' @param object An object of type \code{svc}.
#' @param env Computing environment.
#' @param ... Not used.
.compute.fixedsvc <- function(object, env, ...) {
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
.compute.free <- function(object, env, ...) {
  res <- eval(object$form, env)
  return(res)
}

#' Model implied mean vector
#' @description Compute model implied mean vector.
#' @export
#' @param mod an instance of an \code{svcm} model.
#' @param drop_miss drop expectation for missing values?
#' @return vector of means.
expected_mean <- function(mod, drop_miss = T) {
  M <- Reduce("+", lapply(mod$mcs, .compute, mod$env_comp))
  if(drop_miss) {
    M <- M[mod$dat$keepy]
  }
  return(M)
}

#' Model implied covariance matrix
#' @description Compute model implied covariance matrix.
#' @export
#' @param svcm an instance of an \code{svcm} model.
#' @param drop_miss drop expectation for missing values?
#' @return matrix of (co)variances.
expected_cov <- function(svcm, drop_miss = T) {
  S <- Reduce("+", lapply(svcm$svcs, .compute, svcm$env_comp))
  if(drop_miss) {
    S <- S[svcm$dat$keepy, svcm$dat$keepy]
  }
  return(S)
}

#' Constructor for parameter matrix
#' @description Creates an object of type \code{pm} used to represent parameters in a model.
#' @export
#' @param nrow Number of matrix rows.
#' @param ncol Number of matrix columns.
#' @param labels Vector/matrix with parameter labels. Equal characters defines equality constraints.
#' @param values Vector/matrix with parameter values. Values for free parameters are used as starting values for optimization
#' and values for fixed parameters are constant during optimization.
#' @param free Vector or matrix of logical values defining parameters free or not.
#' @param name Character giving name to model object. Used to reference the parameter matrix in computations.
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
  .new_mc(form = form, Xt = t(X), class = "fixedmc")
}
.new_freemc <- function(form) {
  .new_mc(form = form, class = "free")
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
                        pms = pms, # is it really any point of storing other than pm?
                        svcs = svcs,
                        mcs = mcs,
                        ics = ics,
                        env_comp = new.env(), # Computing environment Should probably inherit prom parent environment to get relationship matrices for fixed types
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


#' Returns parameters of a model
#' @description returns non-duplicated parameters of a \code{svcm} model.
#' @export
#' @param svcm an object of class \code{svcm}.
#' @return vector of parameters.
theta <- function(svcm) {
  theta_start <- unlist(lapply(svcm$pms, .get_free_values))
  names(theta_start) <- unlist(lapply(svcm$pms, .get_free_labels))
  theta_start[!duplicated(names(theta_start))] # Keep only one (first) when equal labels (equality constraints)
}

#' Update parameter matrices in environment from theta
#' @description update all free elements in parameter matrices.
#' @export
#' @param mod An object of type \code{fitm}.
#' @param theta A vector of parameter values.
update_model <- function(mod, theta) {
  .update_pms(mod, theta) # update parameter matrices
  .update_ics(mod) # update / evaluate intermediate computations
}

fit_objective <- function(theta, mod) {
  #.update_pms(svcm, theta) # update parameter matrices
  #.update_ics(svcm) # update / evaluate intermediate computations
  update_model(mod, theta)
  return(objective(mod)) # total mean and covariance are updated in objective.
}

#' Fit a model
#' @description fits a \code{svcm} model.
#' @export
#' @param mod an instance of an \code{svcm} model.
#' @param se should standard errors be computed?
#' @param ... arguments passed to \code{nlminb}.
#' @return an object of class \code{svcm}.
fit_svcm <- function(mod, se = FALSE, ...) {

  if(!inherits(mod, "svcm")) {
    stop("Only objects of type svcm are accepted.")
  }

  if(!is.null(mod$opt)) {
    warning("This model has already been fitted.")
  }



  # optimize model
  #ctrl = get0("control", where = list(...))
  #if(!is.null(ctrl$trace)) {
    #if(ctrl$trace > 0) {
      #cat("\niteration: objective:", names(theta(mod)), "\n")
    #}
  #}
  time_start <- proc.time()
  fit <- nlminb(theta(mod), fit_objective, mod = mod, ...)
  fit$time <- proc.time() - time_start
  if(fit$convergence != 0) {
    warning("Optimization may not have converged.",
            " \nnlminb convergence code: ", fit$convergence,
            " \nnlminb message: ", fit$message)
  }

  mod$opt = fit

  # Hessian at minimum
  if(se) {
    message("Computing standard errors.")
    mod$H <- fd_hess_svcm(mod)
  }
  # This is necessary because they are not updated during optimisation and theta() uses these values.
  ## They are the starting values for optimisation and may be used to continue optimization. No need to update during optim.
  for(i in seq_along(mod$pms)) {
    mod$pms[[i]]$values <- get(mod$pms[[i]]$name, envir = mod$env_comp)
  }

  return(mod)
}

#' Objective function
#' @description Objective function
#' @export
#' @param mod An object of type \code{fitm}
#' @return Twice negative log likelihood.
objective <- function(mod) {

  M <- expected_mean(mod)
  S <- expected_cov(mod)
  y <- mod$dat$y[mod$dat$keepy]

  ch <- Matrix::Cholesky(S)
  rm <- y - M
  r2 <- Matrix::solve(ch, rm)
  deter <- 2 * Matrix::determinant(ch)$modulus
  dev <- log(2 * pi) * length(y) + deter + sum(rm * r2)
  return(dev)
}
# You should make this so that it can be called outside a model call. Either return fit_objective or change this to work on a model
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

#' Compute hessian
#' @description Computes hessian of parameters with finite differences
#' @export
#' @param mod An object of type \code{fitm}
#' @param ... Arguments passed to \code{numDeriv::hessian}.
#' @return Hessian matrix.
fd_hess_svcm <- function(mod, ...) {

  H <- numDeriv::hessian(fit_objective, mod$opt$par, mod = mod)
  dimnames(H) <- list(names(mod$opt$par), names(mod$opt$par))

  # Finite diff have now changed values in pms so need to be reset to par before return
  update_model(mod, mod$opt$par)
  # Dont think this is needed here as this should only be called on fitted models and these should already be updated solution
  #for(i in seq_along(mod$pms)) {
#    mod$pms[[i]]$values <- get(mod$pms[[i]]$name, envir = mod$env_comp)
 # }
  return(H)
}
