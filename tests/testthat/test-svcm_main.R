
# Shared test data -----------------------------------------------------------

make_simple_model <- function(J = 10, seed = 4817) {
  set.seed(seed)
  Y <- matrix(rnorm(J * 2), J, 2, dimnames = list(NULL, c("y1", "y2")))
  R <- Matrix::Diagonal(J)
  X <- matrix(1, J, 1)
  P <- pm(2, 2, c("p11", "p12", "p12", "p22"), TRUE, diag(1, 2), "P")
  U <- pm(2, 1, c("u1", "u2"), TRUE, 0, "U")
  list(Y = Y, R = R, X = X, P = P, U = U)
}


# pm() -----------------------------------------------------------------------

test_that("pm() returns pm class with correct dimensions", {
  p <- pm(2, 3, paste0("l", 1:6), TRUE, 0, "L")
  expect_s3_class(p, "pm")
  expect_equal(dim(p$values), c(2, 3))
  expect_equal(dim(p$free), c(2, 3))
  expect_equal(dim(p$labels), c(2, 3))
  expect_equal(p$name, "L")
})

test_that("pm() recycles scalar free and values", {
  p <- pm(2, 2, paste0("a", 1:4), free = TRUE, values = 0.5, name = "A")
  expect_true(all(p$free))
  expect_true(all(p$values == 0.5))
})

test_that("pm() stores fixed elements correctly", {
  free <- diag(TRUE, 2)
  p <- pm(2, 2, c("a", "b", "c", "d"), free, diag(2, 2), "A")
  expect_equal(p$free, free)
  expect_false(p$free[1, 2])
  expect_false(p$free[2, 1])
})

test_that("pm() errors without a name", {
  expect_snapshot(
    error = TRUE,
    pm(2, 2, c("a", "b", "c", "d"), TRUE, 1)
  )
})


# ic() -----------------------------------------------------------------------

test_that("ic() returns ic class with correct structure", {
  i <- ic(A %*% t(A), name = "B")
  expect_s3_class(i, "ic")
  expect_true(inherits(i, "free"))
  expect_equal(i$name, "B")
  expect_type(i$form, "language")
})

test_that("ic() errors without a name", {
  expect_snapshot(
    error = TRUE,
    ic(A %*% t(A))
  )
})


# const() --------------------------------------------------------------------

test_that("const() returns const class with correct structure", {
  R <- Matrix::Diagonal(10)
  k <- const(R, "R")
  expect_s3_class(k, "const")
  expect_equal(k$name, "R")
  expect_equal(k$value, R)
})

test_that("const() errors without a name", {
  expect_snapshot(error = TRUE, const(Matrix::Diagonal(5)))
})

test_that("const() is assigned into env_comp and visible in expressions", {
  d <- make_simple_model()
  I <- ncol(d$Y)
  Z <- Matrix::Diagonal(I)
  mod <- svcm(d$Y, d$P, d$U,
              const(Z, "Z"),
              svc(Z %*% P %*% t(Z), R = d$R),
              mc(U, X = d$X))
  expect_true(exists("Z", envir = mod$env_comp, inherits = FALSE))
  expect_equal(get("Z", envir = mod$env_comp), Z)
})

test_that("t() resolves to Matrix::t in model expressions", {
  # Verifies that env_comp's parent is the svcm namespace (not the calling
  # frame), so imported functions like Matrix::t are found automatically.
  d <- make_simple_model()
  I <- ncol(d$Y)
  LP <- pm(I, 1, rep("lp", I), TRUE, 1, "LP")
  expect_no_error(
    svcm(d$Y, LP, d$U,
         ic(LP %*% t(LP), name = "P"),
         svc(P, R = d$R),
         mc(U, X = d$X))
  )
})


# svc() ----------------------------------------------------------------------

test_that("svc() without R returns free svc", {
  s <- svc(A + B)
  expect_s3_class(s, "svc")
  expect_true(inherits(s, "free"))
  expect_false(inherits(s, "fixedsvc"))
})

test_that("svc() with R returns fixedsvc", {
  R <- Matrix::Diagonal(10)
  s <- svc(A, R = R)
  expect_s3_class(s, "svc")
  expect_true(inherits(s, "fixedsvc"))
  expect_equal(s$R, R)
})

test_that("svc() stores the expression unevaluated", {
  s <- svc(L %*% S %*% t(L))
  expect_type(s$form, "language")
})

test_that("svc() errors when R is not a Matrix", {
  expect_snapshot(
    error = TRUE,
    svc(A, R = matrix(c(1, 0, 0, 1), 2, 2))
  )
})

test_that("svc() errors when R is not symmetric", {
  R_asym <- Matrix::Matrix(c(1, 2, 0, 1), 2, 2, sparse = TRUE)
  expect_snapshot(
    error = TRUE,
    svc(A, R = R_asym)
  )
})


# mc() -----------------------------------------------------------------------

test_that("mc() without X returns free mc", {
  m <- mc(B)
  expect_s3_class(m, "mc")
  expect_true(inherits(m, "free"))
  expect_false(inherits(m, "fixedmc"))
})

test_that("mc() with X returns fixedmc", {
  X <- matrix(1, 10, 1)
  m <- mc(B, X = X)
  expect_s3_class(m, "mc")
  expect_true(inherits(m, "fixedmc"))
  expect_equal(m$X, X)
})

test_that("mc() stores the expression unevaluated", {
  m <- mc(B %*% t(X))
  expect_type(m$form, "language")
})

test_that("mc() warns when X is rank-deficient", {
  X_rankdef <- matrix(1, 5, 2)  # both columns identical
  expect_snapshot(mc(B, X = X_rankdef))
})


# dat_svcm() -----------------------------------------------------------------

test_that("dat_svcm() returns correct class and structure", {
  Y <- matrix(1:6, 3, 2)
  d <- dat_svcm(Y)
  expect_s3_class(d, "dat_svcm")
  expect_equal(d$Y, Y)
  expect_length(d$y, 6)
  expect_length(d$keepy, 6)
  expect_true(all(d$keepy))
})

test_that("dat_svcm() correctly handles missing values", {
  Y <- matrix(c(1, NA, 3, 4, NA, 6), 3, 2)
  d <- dat_svcm(Y)
  expect_equal(sum(d$keepy), 4)
  expect_length(d$y, 4)
  expect_false(anyNA(d$y))
})

test_that("dat_svcm() stacks columns in column-major order", {
  Y <- matrix(1:6, 3, 2)
  d <- dat_svcm(Y)
  expect_equal(d$y, 1:6)
})


# svcm() ---------------------------------------------------------------------

test_that("svcm() returns svcm class", {
  d <- make_simple_model()
  mod <- svcm(d$Y, d$P, d$U, svc(P, R = d$R), mc(U, X = d$X))
  expect_s3_class(mod, "svcm")
})

test_that("svcm() populates pms, svcs, mcs", {
  d <- make_simple_model()
  mod <- svcm(d$Y, d$P, d$U, svc(P, R = d$R), mc(U, X = d$X))
  expect_length(mod$pms, 2)
  expect_length(mod$svcs, 1)
  expect_length(mod$mcs, 1)
})

test_that("svcm() initialises the computing environment", {
  d <- make_simple_model()
  mod <- svcm(d$Y, d$P, d$U, svc(P, R = d$R), mc(U, X = d$X))
  expect_true(is.environment(mod$env_comp))
  expect_true(exists("P", envir = mod$env_comp))
  expect_true(exists("U", envir = mod$env_comp))
})

test_that("svcm() errors when pm is missing", {
  d <- make_simple_model()
  expect_snapshot(
    error = TRUE,
    svcm(d$Y, svc(P, R = d$R), mc(U, X = d$X))
  )
})

test_that("svcm() errors when svc is missing", {
  d <- make_simple_model()
  expect_snapshot(
    error = TRUE,
    svcm(d$Y, d$P, d$U, mc(U, X = d$X))
  )
})

test_that("svcm() errors when mc is missing", {
  d <- make_simple_model()
  expect_snapshot(
    error = TRUE,
    svcm(d$Y, d$P, d$U, svc(P, R = d$R))
  )
})


# theta() --------------------------------------------------------------------

test_that("theta() returns a named numeric vector", {
  d <- make_simple_model()
  mod <- svcm(d$Y, d$P, d$U, svc(P, R = d$R), mc(U, X = d$X))
  th <- theta(mod)
  expect_type(th, "double")
  expect_named(th)
})

test_that("theta() deduplicates equality-constrained parameters", {
  d <- make_simple_model()
  mod <- svcm(d$Y, d$P, d$U, svc(P, R = d$R), mc(U, X = d$X))
  th <- theta(mod)
  # P has labels c("p11", "p12", "p12", "p22") — p12 appears twice
  expect_equal(sum(names(th) == "p12"), 1)
})

test_that("theta() returns correct parameter count", {
  d <- make_simple_model()
  mod <- svcm(d$Y, d$P, d$U, svc(P, R = d$R), mc(U, X = d$X))
  th <- theta(mod)
  # P: p11, p12, p22 (3 unique); U: u1, u2 (2) => 5 total
  expect_length(th, 5)
})

test_that("theta() returns starting values from pm()", {
  d <- make_simple_model()
  mod <- svcm(d$Y, d$P, d$U, svc(P, R = d$R), mc(U, X = d$X))
  th <- theta(mod)
  expect_equal(th[["p11"]], 1)
  expect_equal(th[["p22"]], 1)
  expect_equal(th[["u1"]], 0)
})


# update_model() -------------------------------------------------------------

test_that("update_model() propagates new values into the environment", {
  d <- make_simple_model()
  mod <- svcm(d$Y, d$P, d$U, svc(P, R = d$R), mc(U, X = d$X))
  new_theta <- theta(mod)
  new_theta["p11"] <- 99
  update_model(mod, new_theta)
  P_in_env <- get("P", envir = mod$env_comp)
  expect_equal(P_in_env[1, 1], 99)
})
