
# Profiled fixed effects fe() ------------------------------------------------
#
# A model whose linear fixed effects are profiled out with fe(X) must reproduce
# exactly the same maximised log-likelihood and parameter estimates as the
# equivalent model that specifies those fixed effects explicitly with a free
# pm() + mc(). fe() gives each response its own coefficient vector (no
# constraints across responses), so the explicit model uses a free B matrix with
# distinct labels per response.

# Simulate a multivariate outcome with response-specific intercepts and slopes.
#   Y[, i] = b0[i] + b1[i] * z + eta + eps[, i]
sim_fe_data <- function(J = 200, I = 3, var_eta = 1, var_eps = 1) {
  z   <- rnorm(J)
  b0  <- c(1, -2, 0.5)[seq_len(I)]
  b1  <- c(0.5, 1.5, -1)[seq_len(I)]
  eta <- rnorm(J, 0, sqrt(var_eta))
  eps <- matrix(rnorm(J * I, 0, sqrt(var_eps)), J, I)
  fixed <- sweep(outer(z, b1), 2, b0, "+")   # J x I: b0[i] + b1[i] z
  Y   <- fixed + matrix(rep(eta, I), J, I) + eps
  list(Y = Y, R = Matrix::Diagonal(J), X = cbind(1, z), I = I)
}

# Variance side shared by both parameterisations (Cholesky / sqrt param).
.var_part <- function(I, R) {
  list(
    pm(I, 1, rep("lp",  I), TRUE,          1,          "LP"),
    pm(I, I, rep("lth", I), diag(TRUE, I), diag(1, I), "LTH"),
    ic(tcrossprod(LP),  name = "P"),
    ic(tcrossprod(LTH), name = "TH"),
    svc(P + TH, R = R)
  )
}

# Explicit model: free B (I x 2), one intercept and one slope per response.
build_explicit <- function(d) {
  I <- d$I
  do.call(svcm, c(
    list(d$Y),
    .var_part(I, d$R),
    list(
      pm(I, 2, sapply(1:2, \(j) paste0("b", 1:I, j)), TRUE, 0, "B"),
      mc(B, X = d$X)
    )
  ))
}

# Profiled model: same fixed effects handled by fe(X).
build_profiled <- function(d) {
  I <- d$I
  do.call(svcm, c(
    list(d$Y),
    .var_part(I, d$R),
    list(fe(d$X, labels = c("intercept", "slope")))
  ))
}


test_that("fe() reproduces the explicit-mean log-likelihood", {
  set.seed(101)
  d <- sim_fe_data()
  fit_e <- fit_svcm(build_explicit(d))
  fit_p <- fit_svcm(build_profiled(d))

  expect_equal(as.numeric(logLik(fit_e)),
               as.numeric(logLik(fit_p)),
               tolerance = 1e-6)
})

test_that("fe() reproduces the variance parameter estimates", {
  set.seed(101)
  d <- sim_fe_data()
  fit_e <- fit_svcm(build_explicit(d))
  fit_p <- fit_svcm(build_profiled(d))

  expect_equal(theta(fit_e)[["lp"]]^2,  theta(fit_p)[["lp"]]^2,  tolerance = 1e-4)
  expect_equal(theta(fit_e)[["lth"]]^2, theta(fit_p)[["lth"]]^2, tolerance = 1e-4)
})

test_that("fe() reproduces the fixed-effect coefficient estimates", {
  set.seed(101)
  d <- sim_fe_data()
  fit_e <- fit_svcm(build_explicit(d))
  fit_p <- fit_svcm(build_profiled(d))

  # Explicit B is stored column-major as b<i><j>: intercepts (j=1) then slopes
  # (j=2). Profiled beta is ordered response-major (per response: intercept,
  # slope). Reorder the explicit estimates to the profiled layout.
  I <- d$I
  B_hat <- theta(fit_e)[paste0("b", rep(1:I, times = 2), rep(1:2, each = I))]
  B_prof_layout <- as.numeric(t(matrix(B_hat, I, 2)))   # response-major

  expect_equal(unname(fit_p$beta), B_prof_layout, tolerance = 1e-4)
})

test_that("logLik df counts profiled coefficients (AIC/BIC comparable)", {
  set.seed(101)
  d <- sim_fe_data()
  fit_e <- fit_svcm(build_explicit(d))
  fit_p <- fit_svcm(build_profiled(d))

  expect_equal(attr(logLik(fit_e), "df"), attr(logLik(fit_p), "df"))
  expect_equal(AIC(fit_e), AIC(fit_p), tolerance = 1e-5)
  expect_equal(BIC(fit_e), BIC(fit_p), tolerance = 1e-5)
})

test_that("expected_mean() includes the profiled fixed-effect contribution", {
  set.seed(101)
  d <- sim_fe_data()
  fit_e <- fit_svcm(build_explicit(d))
  fit_p <- fit_svcm(build_profiled(d))

  expect_equal(as.numeric(expected_mean(fit_p)),
               as.numeric(expected_mean(fit_e)),
               tolerance = 1e-4)
})

test_that("fe() standard errors match the explicit model", {
  set.seed(101)
  d <- sim_fe_data()
  fit_e <- fit_svcm(build_explicit(d), se = TRUE)
  fit_p <- fit_svcm(build_profiled(d), se = TRUE)

  se_e <- sqrt(diag(vcov(fit_e)))
  se_p <- sqrt(diag(vcov(fit_p)))

  # Compare the profiled-coefficient standard errors (reordered as above).
  I <- d$I
  se_B <- se_e[paste0("b", rep(1:I, times = 2), rep(1:2, each = I))]
  se_B_layout <- as.numeric(t(matrix(se_B, I, 2)))
  se_beta <- se_p[fit_p$beta_labels]

  expect_equal(unname(se_beta), se_B_layout, tolerance = 1e-3)
})
