
# Shared helpers -------------------------------------------------------------

# Simulate a random intercept model with known parameters.
# Y[j, i] = eta[j] + eps[j, i], eta ~ N(mean, var_eta), eps ~ N(0, var_eps)
sim_random_intercept <- function(J, I, var_eta = 1, var_eps = 1, mean = 0) {
  eta <- rnorm(J, mean, sqrt(var_eta))
  eps <- matrix(rnorm(J * I, 0, sqrt(var_eps)), J, I)
  Y   <- matrix(rep(eta, I), J, I) + eps
  list(Y = Y, R = Matrix::Diagonal(J), X = matrix(1, J, 1), I = I)
}

# Build the model using a Cholesky / sqrt parameterisation to guarantee
# positive-definite covariance components during unconstrained optimisation:
#
#   LP  is I×1 (all elements equal, labeled "lp")
#   P   = tcrossprod(LP) = LP %*% t(LP)   [I×I rank-1 matrix]
#
#   LTH is I×I diagonal (diagonal labeled "lth")
#   TH  = tcrossprod(LTH) = LTH %*% t(LTH) [I×I diagonal matrix]
#
# True variances are recovered as lp^2 = var_eta and lth^2 = var_eps.
build_random_intercept_model <- function(d) {
  I <- d$I
  svcm(
    d$Y,
    pm(I, 1, rep("lp",  I), TRUE,          1,          "LP"),
    pm(I, I, rep("lth", I), diag(TRUE, I), diag(1, I), "LTH"),
    pm(I, 1, "u",           TRUE,          0,           "U"),
    ic(tcrossprod(LP),  name = "P"),
    ic(tcrossprod(LTH), name = "TH"),
    svc(P + TH, R = d$R),
    mc(U, X = d$X)
  )
}


# fit_svcm() behaviour -------------------------------------------------------

test_that("fit_svcm() returns svcm with opt populated", {
  set.seed(3719)
  d   <- sim_random_intercept(J = 50, I = 4)
  mod <- build_random_intercept_model(d)
  fit <- fit_svcm(mod)

  expect_s3_class(fit, "svcm")
  expect_false(is.null(fit$opt))
  expect_equal(fit$opt$convergence, 0)
})

test_that("fit_svcm() errors on non-svcm input", {
  expect_snapshot(error = TRUE, fit_svcm(list()))
})

test_that("fit_svcm() warns when model is already fitted", {
  set.seed(3719)
  d   <- sim_random_intercept(J = 50, I = 4)
  mod <- build_random_intercept_model(d)
  fit <- fit_svcm(mod)
  # Suppress the return value so only the warning text is snapshotted
  expect_snapshot(invisible(fit_svcm(fit)))
})


# Numerical recovery ---------------------------------------------------------

test_that("fit_svcm() recovers unit variances in a random intercept model", {
  set.seed(3719)
  d   <- sim_random_intercept(J = 500, I = 4, var_eta = 1, var_eps = 1, mean = 0)
  fit <- fit_svcm(build_random_intercept_model(d))
  th  <- theta(fit)

  # lp^2 = var_eta, lth^2 = var_eps
  expect_equal(th[["lp"]]^2,  1, tolerance = 0.1)
  expect_equal(th[["lth"]]^2, 1, tolerance = 0.1)
  expect_equal(th[["u"]],     0, tolerance = 0.1)
})

test_that("fit_svcm() recovers non-unit variances (var_eta=2, var_eps=0.5, mean=3)", {
  set.seed(8204)
  d   <- sim_random_intercept(J = 500, I = 4, var_eta = 2, var_eps = 0.5, mean = 3)
  fit <- fit_svcm(build_random_intercept_model(d))
  th  <- theta(fit)

  expect_equal(th[["lp"]]^2,  2,   tolerance = 0.2)
  expect_equal(th[["lth"]]^2, 0.5, tolerance = 0.1)
  expect_equal(th[["u"]],     3,   tolerance = 0.15)
})


# logLik() -------------------------------------------------------------------

test_that("logLik() returns a finite negative value after fitting", {
  set.seed(3719)
  d   <- sim_random_intercept(J = 50, I = 4)
  fit <- fit_svcm(build_random_intercept_model(d))
  ll  <- logLik(fit)

  expect_s3_class(ll, "logLik")
  expect_true(is.finite(ll))
  expect_lt(as.numeric(ll), 0)
})

test_that("logLik() sets nobs and df attributes", {
  set.seed(3719)
  d   <- sim_random_intercept(J = 50, I = 4)
  fit <- fit_svcm(build_random_intercept_model(d))
  ll  <- logLik(fit)

  expect_equal(attr(ll, "nobs"), 50 * 4)
  expect_equal(attr(ll, "df"),   3)  # lp, lth, u
})

test_that("logLik() errors when model has not been fitted", {
  set.seed(3719)
  d   <- sim_random_intercept(J = 50, I = 4)
  mod <- build_random_intercept_model(d)
  expect_snapshot(error = TRUE, logLik(mod))
})


# ic() environment path ------------------------------------------------------

test_that("with() returns P = tcrossprod(LP) with all elements equal to lp^2", {
  set.seed(3719)
  d   <- sim_random_intercept(J = 500, I = 4)
  fit <- fit_svcm(build_random_intercept_model(d))

  lp_hat <- theta(fit)[["lp"]]
  P_hat  <- as.matrix(with(fit, P))

  expect_equal(dim(P_hat), c(4L, 4L))
  expect_equal(P_hat, matrix(lp_hat^2, 4, 4), tolerance = 1e-10)
})
