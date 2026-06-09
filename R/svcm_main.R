
#' @importFrom Matrix t
#' @importFrom methods as
#' @importFrom stats AIC BIC logLik nlminb vcov
"_PACKAGE"

#' Free values in parameter matrix
#' @param pm an instance of an \code{pm} object.
#' @return vector of values of free elements in parameter matrix (column wise)
.get_free_values <- function(pm) {
  return(pm$values[pm$free])
}

#' Labels of free values in parameter matrix
#' @param pm an instance of an \code{pm} object.
#' @return vector of labels of free elements in parameter matrix (column wise).
.get_free_labels <- function(pm) {
  return(pm$labels[pm$free])
}

#' Updates values in parameter matrix
#' @description Used during optimization for updating free values in a parameter matrix in both computing environment and model object.
#' @param pm an instance of an \code{pm} object.
#' @param theta vector of free model parameters.
#' @param env_comp computing environment where the expression is evaluated and stored.
.update_pm <- function(pm, theta, env_comp) {
  pos_free_pm <- which(pm$free) # positions of free values in pm
  ind_theta <- match(pm$labels[pos_free_pm], names(theta)) # Positions of pos_free_pm in values
  pm$values[pos_free_pm] <- theta[ind_theta]
  assign(pm$name, pm$values, envir = env_comp) # update environment
}

#' Updates all \code{pm} objects of a model
#' @description Used during optimization for updating parameter matrices in computing environment and model object.
#' @param m an instance of an \code{svcm} model.
#' @param theta vector of free model parameters.
.update_pms <- function(m, theta) {
  lapply(m$pms, .update_pm, theta, m$env_comp)
}

#' Updates intermediate computations
#' @description Used during optimization for computing expressions in \code{ic} objects.
#' @param ic and instance of an \code{ic} object.
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
#' @param object An object of type \code{mc} or \code{svc}.
#' @param ... Further function arguments.
.compute <- function(object, ...) {
  UseMethod(".compute")
}

#' Compute function for mc object
#' @description Internal function used to evaluate expressions containing objects of type \code{fixedmc}.
#' @export
#' @param object An object of type \code{mc}.
#' @param env Computing environment.
#' @param ... Not used.
#' @return vector/matrix with means for an object of type \code{mc}.
.compute.fixedmc <- function(object, env, ...) {
  mci <- eval(object$form, env)
  mu <- object$X %*% Matrix::t(mci) # Column means
  return(mu)
}

#' Compute function for svc objects
#' @description Internal function used to evaluate expressions containing objects of type \code{fixedsvc}.
#' @export
#' @param object An object of type \code{svc}.
#' @param env Computing environment.
#' @param ... Not used.
#' @return Matrix with covariances for an object of type \code{svc}.
.compute.fixedsvc <- function(object, env, ...) {
  vci <- eval(object$form, env)
  if(is.null(object$tmpl)) {
    return(vci %x% object$R) # not prepared: fall back to a full Kronecker product
  }
  # The sparsity pattern of vci %x% R is fixed (R is constant, vci is dense), so
  # only the numeric values change between iterations. Refill them directly
  # instead of rebuilding the sparse Kronecker product from scratch.
  sigma <- object$tmpl
  sigma@x <- as.numeric(vci)[object$vidx] * object$rval
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
expected_mean <- function(mod, drop_miss = TRUE) {
  Mm <- Reduce("+", lapply(mod$mcs, .compute, mod$env_comp))
  M  <- as.vector(Mm) # vec(Mm) = vec(X %*% t(B))
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
expected_cov <- function(svcm, drop_miss = TRUE) {
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
#' library(svcm)
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

#' Constructor for constant object
#' @description Creates a constant object for use in model expressions. Constants
#' are values (e.g. sparse matrices) that are referenced in \code{svc}, \code{mc},
#' or \code{ic} expressions but are not model parameters.
#' @export
#' @param value The constant value.
#' @param name Character giving name to reference the constant in model expressions.
#' @return An object of class \code{const}.
const <- function(value, name) {
  if (missing(name) || !is.character(name) || length(name) != 1L) {
    stop("A name is required.")
  }
  structure(list(value = value, name = name), class = "const")
}

#' Constructor for variance component object
#' @description Creates a variance component object as a function of model objects and (optionally) a relationship matrix.
#' @export
#' @param form An expression for the covariance model.
#' @param R A relationship matrix. Should be a sparse matrix from the \code{Matrix} package.
#' @return An object of class \code{svc}.
#' @examples
#' library(svcm)
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

# Precompute the fixed Kronecker pattern of a fixedsvc term so that
# .compute.fixedsvc can refill only the numeric values during optimization.
# Stores on the object:
#   tmpl - a CsparseMatrix with the full (dense-block) pattern of vci %x% R,
#   rval - the value of R at each stored nonzero,
#   vidx - the column-major index into vec(vci) scaling each stored nonzero.
.prepare_svc <- function(object, env) {
  if(!inherits(object, "fixedsvc")) {
    return(object)
  }
  vci <- eval(object$form, env)
  p <- nrow(vci)
  Rg <- as(as(object$R, "CsparseMatrix"), "generalMatrix")
  ones <- matrix(1, p, p)
  tmpl <- ones %x% Rg          # full block pattern; @x holds the R values
  Vpos <- matrix(seq_len(p * p), p, p)
  Rones <- Rg
  Rones@x <- rep(1, length(Rg@x))
  object$tmpl <- tmpl
  object$rval <- tmpl@x
  object$vidx <- as.integer((Vpos %x% Rones)@x)
  return(object)
}

#' Constructor for mean component object
#' @description Creates a mean component object used to define (part of) the mean structure in the model.
#' @export
#' @param form An expression for the mean model.
#' @param X A design matrix.
#' @return An object of class \code{mc}.
#' @examples
#' library(svcm)
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
    stopifnot(!anyNA(X)) # is.numeric(X) fail for Matrix but not matrix
    if(!is.null(dim(X))) { # Is matrix, not vector
      rm <- Matrix::rankMatrix(X)
      if(rm < min(dim(X))) {
        warning("The matrix X is not of full rank.")
      }
    }
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
  .new_mc(form = form, X = X, class = "fixedmc")
}
.new_freemc <- function(form) {
  .new_mc(form = form, class = "free")
}

#' Creates a model
#' @description Creates a new model
#' @export
#' @param Y Matrix of data described by model.
#' @param ... All relevant model objects: \code{pm}, \code{svc}, \code{mc},
#'   \code{ic}, and \code{const}.
#' @return An object of type \code{svcm}.
#' @examples
#' \donttest{
#' # Genetic factor model: y_j = b + l * eta_j + delta_j
#' # eta_j = a_j + e_j,  Cov(a) = sigma_a^2 * A,  Cov(e) = sigma_e^2 * I
#' # Cov(y) = (l l' sigma_a^2) x A  +  (l l' sigma_e^2 + TH) x I
#'
#' # Small genetic relationship matrix (two pedigree structures)
#' str1 <- matrix(c(1, 1/4, 1/4, 1/4, 1, 1/2, 1/4, 1/2, 1), 3, 3)
#' str2 <- matrix(c(1, 1/2, 1/8, 1/8, 1/2, 1, 1/8, 1/8,
#'                  1/8, 1/8, 1, 1/2, 1/8, 1/8, 1/2, 1), 4, 4)
#' A <- Matrix::bdiag(Matrix::bdiag(replicate(5, str1, simplify = FALSE)),
#'                    Matrix::bdiag(replicate(5, str2, simplify = FALSE)))
#' N <- nrow(A)
#'
#' # Simulate data
#' set.seed(6)
#' Lc  <- chol(A)
#' a   <- sqrt(2) * t(Lc) %*% rnorm(N)
#' e   <- rnorm(N, 0, sqrt(2))
#' Y   <- rep(1, N) %*% t(c(2, 2, 4, 4)) +
#'        (a + e) %*% t(c(1, 0.5, 0.5, 0.8)) +
#'        matrix(rnorm(N * 4), N, 4)
#'
#' # Model 1: R= argument handles the Kronecker product implicitly
#' m1 <- svcm(Y,
#'   pm(4, 1, paste0("l", 1:4), c(FALSE, TRUE, TRUE, TRUE), c(1, .5, .5, .5), "l"),
#'   pm(1, 1, "Sa1", TRUE, 1, "Sa"),
#'   pm(1, 1, "Se1", TRUE, 1, "Se"),
#'   pm(4, 4, sapply(1:4, \(i) paste0("th", 1:4, i)), diag(TRUE, 4), diag(4), "TH"),
#'   pm(4, 1, paste0("b", 1:4), TRUE, 0, "b"),
#'   svc(l %*% Sa %*% t(l), R = A),
#'   svc(l %*% Se %*% t(l) + TH, R = Matrix::Diagonal(N)),
#'   mc(b, X = rep(1, N))
#' )
#'
#' # Model 2: equivalent, using const() to pass A and I as named constants
#' # and writing the Kronecker products explicitly in the svc() expression
#' m2 <- svcm(Y,
#'   pm(4, 1, paste0("l", 1:4), c(FALSE, TRUE, TRUE, TRUE), c(1, .5, .5, .5), "l"),
#'   pm(1, 1, "Sa1", TRUE, 1, "Sa"),
#'   pm(1, 1, "Se1", TRUE, 1, "Se"),
#'   pm(4, 4, sapply(1:4, \(i) paste0("th", 1:4, i)), diag(TRUE, 4), diag(4), "TH"),
#'   pm(4, 1, paste0("b", 1:4), TRUE, 0, "b"),
#'   const(A, "A"),
#'   const(Matrix::Diagonal(N), "I"),
#'   svc(kronecker(l %*% Sa %*% t(l), A) + kronecker(l %*% Se %*% t(l) + TH, I)),
#'   mc(b, X = rep(1, N))
#' )
#'
#' # Both specifications produce identical starting objectives
#' stopifnot(isTRUE(all.equal(objective(m1), objective(m2))))
#'
#' m1 <- fit_svcm(m1)
#' m2 <- fit_svcm(m2)
#' stopifnot(isTRUE(all.equal(logLik(m1), logLik(m2))))
#' }
svcm <- function(Y, ...) {
  dots   <- list(...)
  pms    <- dots[sapply(dots, inherits, "pm")]
  svcs   <- dots[sapply(dots, inherits, "svc")]
  mcs    <- dots[sapply(dots, inherits, "mc")]
  ics    <- dots[sapply(dots, inherits, "ic")]
  consts <- dots[sapply(dots, inherits, "const")]
  if (length(pms) == 0 || length(svcs) == 0 || length(mcs) == 0) {
    stop("At least one pm, svc and mc object must be supplied.")
  }

  # Use the package namespace as the parent of the computing environment so
  # that all imported functions (e.g. Matrix::t) are visible in expressions
  # evaluated there, without requiring the user to attach Matrix explicitly.
  # getNamespace() is used rather than environment() because environment()
  # inside a function returns the execution frame, not the package namespace.
  ret <- structure(list(dat    = dat_svcm(Y),
                        pms    = pms,
                        svcs   = svcs,
                        mcs    = mcs,
                        ics    = ics,
                        consts = consts,
                        env_comp = new.env(parent = getNamespace("svcm")),
                        opt    = NULL,
                        H      = NULL),
                   class = "svcm")
  # Seed constants into the computing environment before any evaluation.
  lapply(consts, \(c) assign(c$name, c$value, envir = ret$env_comp))
  update_model(ret, theta(ret))
  # Precompute the fixed Kronecker patterns so each optimization step only
  # refreshes the numeric values of the variance components.
  ret$svcs <- lapply(ret$svcs, .prepare_svc, env = ret$env_comp)
  return(ret)
}

#' Create a object for storing data
#' @description creates a \code{dat_svcm} object for the response variable.
#' @export
#' @param Y response data described by model.
#' @return an object of class \code{dat_svcm}.
dat_svcm <- function(Y) {
  # Stack Y column-major: var1[1], var1[2], ..., var1[n], var2[1], ...
  y <- c(Y)
  # Find positions of non-missing values
  keepy <- !is.na(y)
  structure(list(Y = Y,
                 y = y[keepy],
                 keepy = keepy),
            class = "dat_svcm")
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
#' @param mod An object of type \code{svcm}.
#' @param theta A vector of parameter values.
update_model <- function(mod, theta) {
  .update_pms(mod, theta) # update parameter matrices
  .update_ics(mod) # update / evaluate intermediate computations
}

fit_objective <- function(theta, mod) {
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
  dots <- list(...)
  th <- theta(mod)
  if("control" %in% names(dots)) {
    ctrl <- dots[["control"]]
    if(!is.null(ctrl$trace) && ctrl$trace > 0) {
      cat("\niter: objective:", names(th), "\n", sep = "\t")
    }
  }
  time_start <- proc.time()
  fit <- nlminb(th, fit_objective, mod = mod, ...)
  fit$time <- proc.time() - time_start
  if(fit$convergence != 0) {
    warning("Optimization may not have converged.",
            " \nnlminb convergence code: ", fit$convergence,
            " \nnlminb message: ", fit$message)
  }
  mod$opt <- fit

  # Hessian at minimum
  if(se) {
    message("Computing standard errors.")
    mod$H <- fd_hess_svcm(mod)
  }
  # pm$values are not updated during optimisation; sync them back from env_comp
  # so that theta() returns the fitted values and the model can be re-fitted.
  for(i in seq_along(mod$pms)) {
    mod$pms[[i]]$values <- get(mod$pms[[i]]$name, envir = mod$env_comp)
  }
  return(mod)
}

#' Objective function
#' @description Objective function
#' @export
#' @param mod An object of type \code{svcm}
#' @return Twice negative log likelihood.
objective <- function(mod) {
  M <- expected_mean(mod)
  S <- expected_cov(mod)

  cS <- Matrix::Cholesky(S)
  r <- mod$dat$y - M
  iSr <- Matrix::solve(cS, r) # inv(S) %*% r
  ld <- Matrix::determinant(cS, sqrt = FALSE)$modulus # logdet(S)
  dev <- log(2 * pi) * length(r) + ld + sum(r * iSr)
  return(dev)
}

#' Compute hessian
#' @description Computes hessian of parameters with finite differences.
#' @export
#' @param mod An object of type \code{svcm}
#' @param ... Arguments passed to \code{numDeriv::hessian}.
#' @return Hessian matrix.
fd_hess_svcm <- function(mod, ...) {
  # Work on a clone so finite differencing does not mutate the caller's model.
  mod <- .clone_model(mod)
  H <- numDeriv::hessian(fit_objective, mod$opt$par, mod = mod)
  pnms <- names(mod$opt$par)
  dimnames(H) <- list(pnms, pnms)
  return(H)
}

#' Compute jacobian of expression
#' @description Computes jacobian of expression containing parameter matrices with finite differences. The output should be a vector.
#' Note that the jacobian is computed with regard to all free parameters in the model.
#' @export
#' @param mod An object of type \code{svcm}
#' @param form An expression for the calculation.
#' @param ... Arguments passed to \code{numDeriv::jacobian}.
#' @return Jacobian matrix.
fd_jacobian_svcm <- function(mod, form, ...) {
  # Work on a clone so finite differencing does not mutate the caller's model.
  mod_tmp <- .clone_model(mod)
  prs <- substitute(form)
  test_form <- eval(prs, mod_tmp$env_comp)
  if(!is.vector(test_form)) {
    stop("Only implemented for vector-valued functions/expressions")
  }
  fu <- function(th) {
    update_model(mod_tmp, th)
    res <- eval(prs, envir = mod_tmp$env_comp)
  }
  J <- numDeriv::jacobian(fu, theta(mod_tmp))
  colnames(J) <- names(theta(mod_tmp))
  return(J)
}

# Return a shallow copy of a model whose computing environment is an independent
# copy. Functions that perturb the parameters (finite differencing) can then
# mutate the clone freely without affecting the caller's model.
.clone_model <- function(mod) {
  mod$env_comp <- list2env(as.list.environment(mod$env_comp, all.names = TRUE),
                           parent = parent.env(mod$env_comp))
  return(mod)
}
