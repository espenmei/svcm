
# Delta method: Cholesky SEs == direct SEs ------------------------------------
#
# The same single-factor CFA is fitted twice:
#   direct:   P and TH are free parameters (p, th1..4)
#   Cholesky: LP and LTH are free parameters (lp, lth1..4), P = LP %*% t(LP),
#             TH = LTH %*% t(LTH)
#
# The two fits maximise the same likelihood, so they produce identical point
# estimates (in direct-parameterisation units) and logLikelihoods.
#
# For the SEs we apply the delta method to the Cholesky fit:
#   J  = fd_jacobian_svcm(fit_chol, <expression>)  -- d(expression)/d(theta)
#   C  = vcov(fit_chol)                            -- solve(0.5 * H)
#   SE = sqrt(diag(J %*% C %*% t(J)))
#
# Note: using H directly (J %*% H %*% t(J)) is wrong because H is the Hessian
# of -2*logL (precision-like), not the covariance. vcov = solve(0.5 * H).

set.seed(5129)
I <- 4L; J <- 500L

eta_sim  <- rnorm(J, 0, sqrt(1.2))
eps_sim  <- sweep(matrix(rnorm(J * I), J, I), 2, sqrt(c(0.5, 0.8, 0.6, 0.7)), `*`)
Y_delta  <- outer(eta_sim, c(1, 0.7, 0.8, 0.6)) +
            sweep(eps_sim, 2, c(0.3, 0.5, -0.2, 0.1), `+`)
R <- Matrix::Diagonal(J)
X <- matrix(1, J, 1)

# Direct parameterisation: estimate p and th1..4 directly
fit_dir <- svcm(
  Y_delta,
  pm(I, 1, paste0("l",  1:I), c(FALSE, TRUE, TRUE, TRUE), 1,         "L"),
  pm(1, 1, "p",               TRUE,                       1,         "P"),
  pm(I, I, paste0("th", 1:I), diag(TRUE, I), diag(1, I),              "TH"),
  pm(I, 1, paste0("u",  1:I), TRUE,          0,                       "U"),
  svc(L %*% P %*% t(L) + TH, R = R),
  mc(U, X = X)
) |> fit_svcm(se = TRUE)

# Cholesky parameterisation: estimate lp and lth1..4, with
#   P  = LP  %*% t(LP)  = tcrossprod(LP)   (1x1 scalar)
#   TH = LTH %*% t(LTH) = tcrossprod(LTH)  (diagonal 4x4)
fit_chol <- svcm(
  Y_delta,
  pm(I, 1, paste0("l",  1:I), c(FALSE, TRUE, TRUE, TRUE), 1,         "L"),
  pm(1, 1, "lp",              TRUE,                       1,         "LP"),
  pm(I, I, paste0("lth", 1:I), diag(TRUE, I), diag(1, I),             "LTH"),
  pm(I, 1, paste0("u",  1:I), TRUE,           0,                      "U"),
  ic(tcrossprod(LP),  name = "P"),
  ic(tcrossprod(LTH), name = "TH"),
  svc(L %*% P %*% t(L) + TH, R = R),
  mc(U, X = X)
) |> fit_svcm(se = TRUE)

C <- vcov(fit_chol)  # = solve(0.5 * fd_hess_svcm(fit_chol))


# Identical likelihoods -------------------------------------------------------

test_that("direct and Cholesky parameterisations give identical logLik", {
  expect_equal(
    as.numeric(logLik(fit_chol)),
    as.numeric(logLik(fit_dir)),
    tolerance = 1e-6
  )
})


# Equivalent point estimates --------------------------------------------------

test_that("Cholesky estimates recover direct estimates after squaring", {
  th_c <- theta(fit_chol)
  th_d <- theta(fit_dir)

  # Small differences (~1e-5) are expected: both optimisers find the same MLE
  # but converge along different paths due to the reparameterisation.
  expect_equal(th_c[["lp"]]^2,  th_d[["p"]],   tolerance = 1e-4)
  for (i in 1:I) {
    expect_equal(th_c[[paste0("lth", i)]]^2, th_d[[paste0("th", i)]], tolerance = 1e-4,
                 label = paste0("lth", i, "^2 == th", i))
  }
})


# Delta-method SEs ------------------------------------------------------------
#
# fd_jacobian_svcm(mod, expr) returns d(expr)/d(theta): a (k x n) matrix
# where k = length of the expression and n = number of free parameters.
#
# Delta method: Var(g(theta)) = J * vcov * J^T
#               SE(g(theta))  = sqrt(diag(J * vcov * J^T))

test_that("delta-method SE for factor variance matches direct SE", {
  # J: d(lp^2) / d(theta) -- scalar expression, so J is 1 x n_params
  J_p <- fd_jacobian_svcm(fit_chol, as.vector(LP %*% t(LP)))
  se_p_delta  <- sqrt(as.numeric(J_p %*% C %*% t(J_p)))
  se_p_direct <- sqrt(diag(vcov(fit_dir)))[["p"]]

  expect_equal(se_p_delta, se_p_direct, tolerance = 1e-4)
})

test_that("delta-method SEs for residual variances match direct SEs", {
  # J: d(diag(LTH %*% t(LTH))) / d(theta) -- 4-element expression, so J is 4 x n_params
  J_th <- fd_jacobian_svcm(fit_chol, diag(as.matrix(LTH %*% t(LTH))))
  se_th_delta  <- sqrt(diag(J_th %*% C %*% t(J_th)))
  se_th_direct <- sqrt(diag(vcov(fit_dir)))[paste0("th", 1:I)]

  for (i in seq_len(I)) {
    expect_equal(se_th_delta[[i]], se_th_direct[[i]], tolerance = 1e-4,
                 label = paste0("th", i, " SE"))
  }
})
