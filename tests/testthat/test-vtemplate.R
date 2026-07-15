
# Precomputed marginal covariance pattern (V template) ------------------------

make_two_term_model <- function(N_blocks = 50, P = 4, seed = 42) {
  set.seed(seed)
  str1 <- matrix(c(1, 1/4, 1/4, 1/4, 1, 1/2, 1/4, 1/2, 1), 3, 3)
  A <- Matrix::bdiag(replicate(N_blocks, str1, simplify = FALSE))
  N <- nrow(A)
  I <- Matrix::Diagonal(N)
  Y <- matrix(rnorm(N * P), N, P)
  svcm(Y,
    pm(P, 1, paste0("l", 1:P), c(FALSE, rep(TRUE, P - 1)), c(1, .5, .5, .5), "l"),
    pm(1, 1, "Sa1", TRUE, 1, "Sa"),
    pm(1, 1, "Se1", TRUE, 1, "Se"),
    pm(P, P, sapply(1:P, \(i) paste0("th", 1:P, i)), diag(TRUE, P), diag(P), "TH"),
    pm(P, 1, paste0("b", 1:P), TRUE, 0, "b"),
    svc(l %*% Sa %*% Matrix::t(l), R = A),
    svc(l %*% Se %*% Matrix::t(l) + TH, R = I),
    mc(b, X = rep(1, N))
  )
}

test_that("svcm() precomputes the V template for all-fixedsvc models", {
  mod <- make_two_term_model(N_blocks = 5)
  expect_false(is.null(mod$Vtmpl))
  expect_s4_class(mod$Vtmpl, "dgCMatrix")
  expect_length(mod$vc_idx, length(mod$svcs))
})

test_that("svcm() skips the V template when a free svc is present", {
  set.seed(1)
  J <- 10
  Y <- matrix(rnorm(J * 2), J, 2)
  R <- Matrix::Diagonal(J)
  P <- pm(2, 2, c("p11", "p12", "p12", "p22"), TRUE, diag(1, 2), "P")
  U <- pm(2, 1, c("u1", "u2"), TRUE, 0, "U")
  mod <- svcm(Y, P, U,
              const(R, "A"),
              svc(P %x% A),
              mc(U, X = matrix(1, J, 1)))
  expect_null(mod$Vtmpl)
  expect_no_error(objective(mod))
})

test_that("expected_cov() with V template matches naive summation", {
  mod <- make_two_term_model(N_blocks = 10)
  S_fast <- expected_cov(mod)
  S_naive <- Reduce("+", lapply(mod$svcs, .compute, mod$env_comp))
  expect_equal(as.matrix(S_fast), as.matrix(S_naive))

  # also after a parameter update
  th <- theta(mod)
  th[] <- th + 0.1
  update_model(mod, th)
  S_fast <- expected_cov(mod)
  S_naive <- Reduce("+", lapply(mod$svcs, .compute, mod$env_comp))
  expect_equal(as.matrix(S_fast), as.matrix(S_naive))
})

test_that("expected_cov() with V template handles missing data", {
  set.seed(7)
  str1 <- matrix(c(1, 1/4, 1/4, 1/4, 1, 1/2, 1/4, 1/2, 1), 3, 3)
  A <- Matrix::bdiag(replicate(5, str1, simplify = FALSE))
  N <- nrow(A)
  Y <- matrix(rnorm(N * 2), N, 2)
  Y[c(2, 9), 1] <- NA
  P <- pm(2, 2, c("p11", "p12", "p12", "p22"), TRUE, diag(1, 2), "P")
  U <- pm(2, 1, c("u1", "u2"), TRUE, 0, "U")
  mod <- svcm(Y, P, U, svc(P, R = A), mc(U, X = matrix(1, N, 1)))
  expect_false(is.null(mod$Vtmpl))
  S <- expected_cov(mod)
  expect_equal(dim(S), rep(sum(mod$dat$keepy), 2))
  expect_no_error(objective(mod))
})

test_that("fitting with V template reproduces the naive fit", {
  mod <- make_two_term_model(N_blocks = 10)
  mod_naive <- mod
  mod_naive$Vtmpl <- NULL
  expect_equal(objective(mod), objective(mod_naive))
  fit <- suppressWarnings(fit_svcm(mod, control = list(iter.max = 500, eval.max = 500)))
  fit_naive <- suppressWarnings(fit_svcm(mod_naive, control = list(iter.max = 500, eval.max = 500)))
  expect_equal(fit$opt$objective, fit_naive$opt$objective, tolerance = 1e-8)
  expect_equal(coef(fit), coef(fit_naive), tolerance = 1e-6)
})

test_that("V template refill is not slower than naive summation", {
  skip_on_cran()
  mod <- make_two_term_model(N_blocks = 400)
  mod_naive <- mod
  mod_naive$Vtmpl <- NULL

  n_rep <- 30
  # warm up dispatch
  invisible(expected_cov(mod)); invisible(expected_cov(mod_naive))
  t_fast <- system.time(
    for (i in seq_len(n_rep)) invisible(expected_cov(mod))
  )["elapsed"]
  t_naive <- system.time(
    for (i in seq_len(n_rep)) invisible(expected_cov(mod_naive))
  )["elapsed"]
  # Generous margin to absorb timer noise; the point is to catch regressions
  # where the template path would be substantially slower than naive Reduce.
  expect_lt(t_fast, t_naive * 1.25 + 0.05)
})
