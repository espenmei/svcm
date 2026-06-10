
# Equivalence with lavaan: single-factor CFA ----------------------------------
#
# Model: Y[j] = L * eta[j] + u + eps[j]
#   eta[j] ~ N(0, p),  eps[j,i] ~ N(0, th_i) independently
#
# The same model is fitted by both packages using full-information ML.
# Point estimates should agree to ~1e-5; SEs differ by at most ~1e-3
# because lavaan uses an analytical information matrix whereas svcm uses
# a numerical Hessian.

skip_if_not_installed("lavaan")

# Simulate once and share across tests ---------------------------------------

set.seed(5129)
I <- 4L
J <- 500L

eta_sim <- rnorm(J, 0, sqrt(1.2))
eps_sim <- sweep(matrix(rnorm(J * I), J, I), 2, sqrt(c(0.5, 0.8, 0.6, 0.7)), `*`)
Y_fac   <- outer(eta_sim, c(1, 0.7, 0.8, 0.6)) +
           sweep(eps_sim, 2, c(0.3, 0.5, -0.2, 0.1), `+`)
colnames(Y_fac) <- paste0("y", 1:I)

# lavaan: standard single-factor CFA with mean structure
lav_mod <- "eta =~ 1*y1 + l2*y2 + l3*y3 + l4*y4"
fit_lav <- lavaan::cfa(lav_mod, as.data.frame(Y_fac), meanstructure = TRUE)
pe      <- lavaan::parameterEstimates(fit_lav)  # estimates and SEs

# svcm: equivalent model using direct parameterisation
# (avoids Cholesky so SE units match lavaan directly)
fit_cfa <- svcm(
  Y_fac,
  pm(I, 1, paste0("l",  1:I), c(FALSE, TRUE, TRUE, TRUE), 1, "L"),
  pm(1, 1, "p",               TRUE,                      1, "P"),
  pm(I, I, paste0("th", 1:I), diag(TRUE, I), diag(1, I),   "TH"),
  pm(I, 1, paste0("u",  1:I), TRUE,          0,             "U"),
  svc(L %*% P %*% Matrix::t(L) + TH, R = Matrix::Diagonal(J)),
  mc(U, X = matrix(1, J, 1))
) |> fit_svcm(se = TRUE)

th_cfa <- theta(fit_cfa)
se_cfa <- sqrt(diag(vcov(fit_cfa)))

# Helpers for pulling lavaan estimates / SEs by parameter type
lav_est <- function(op, lhs, rhs) pe$est[pe$op == op & pe$lhs == lhs & pe$rhs == rhs]
lav_se  <- function(op, lhs, rhs) pe$se[ pe$op == op & pe$lhs == lhs & pe$rhs == rhs]


# Tests ----------------------------------------------------------------------

test_that("svcm and lavaan log-likelihoods agree for factor model", {
  expect_equal(
    as.numeric(logLik(fit_cfa)),
    as.numeric(lavaan::logLik(fit_lav)),
    tolerance = 1e-3
  )
})

test_that("svcm and lavaan factor loading estimates and SEs agree", {
  for (i in 2:I) {
    lname <- paste0("l", i)
    yname <- paste0("y", i)
    expect_equal(th_cfa[[lname]], lav_est("=~", "eta", yname), tolerance = 1e-4,
                 label = paste(lname, "estimate"))
    expect_equal(se_cfa[[lname]], lav_se("=~", "eta", yname),  tolerance = 0.03,
                 label = paste(lname, "SE"))
  }
})

test_that("svcm and lavaan factor variance estimate and SE agree", {
  expect_equal(th_cfa[["p"]], lav_est("~~", "eta", "eta"), tolerance = 1e-4)
  expect_equal(se_cfa[["p"]], lav_se( "~~", "eta", "eta"), tolerance = 0.03)
})

test_that("svcm and lavaan residual variance estimates and SEs agree", {
  for (i in 1:I) {
    tname <- paste0("th", i)
    yname <- paste0("y", i)
    expect_equal(th_cfa[[tname]], lav_est("~~", yname, yname), tolerance = 1e-4,
                 label = paste(tname, "estimate"))
    expect_equal(se_cfa[[tname]], lav_se( "~~", yname, yname), tolerance = 0.03,
                 label = paste(tname, "SE"))
  }
})

test_that("svcm and lavaan intercept estimates and SEs agree", {
  for (i in 1:I) {
    uname <- paste0("u", i)
    yname <- paste0("y", i)
    # lavaan stores intercepts with rhs == "" (not the variable name)
    expect_equal(th_cfa[[uname]], lav_est("~1", yname, ""), tolerance = 1e-4,
                 label = paste(uname, "estimate"))
    expect_equal(se_cfa[[uname]], lav_se( "~1", yname, ""), tolerance = 0.03,
                 label = paste(uname, "SE"))
  }
})


# Missing data (MCAR) --------------------------------------------------------
#
# ~10 % of cells set to NA at random. svcm marginalises over missing values
# via the observed submatrix of the covariance; lavaan uses FIML
# (missing = "ML"). Both maximise the same likelihood, so estimates and SEs
# should agree to the same precision as the complete-data case.
# Note: with missing = "ML" lavaan also uses a numerical information matrix,
# which explains why SEs agree even more tightly here (~1e-6) than in the
# complete-data tests where lavaan uses an analytical matrix (~1e-3).

set.seed(2847)
Y_miss           <- Y_fac
Y_miss[sample(length(Y_miss), size = floor(0.10 * length(Y_miss)))] <- NA

fit_lav_miss <- lavaan::cfa(lav_mod, as.data.frame(Y_miss),
                            meanstructure = TRUE, missing = "ML")
pe_miss      <- lavaan::parameterEstimates(fit_lav_miss)

fit_cfa_miss <- svcm(
  Y_miss,
  pm(I, 1, paste0("l",  1:I), c(FALSE, TRUE, TRUE, TRUE), 1, "L"),
  pm(1, 1, "p",               TRUE,                      1, "P"),
  pm(I, I, paste0("th", 1:I), diag(TRUE, I), diag(1, I),   "TH"),
  pm(I, 1, paste0("u",  1:I), TRUE,          0,             "U"),
  svc(L %*% P %*% Matrix::t(L) + TH, R = Matrix::Diagonal(J)),
  mc(U, X = matrix(1, J, 1))
) |> fit_svcm(se = TRUE)

th_cfa_miss <- theta(fit_cfa_miss)
se_cfa_miss <- sqrt(diag(vcov(fit_cfa_miss)))

lav_est_m <- function(op, lhs, rhs) pe_miss$est[pe_miss$op == op & pe_miss$lhs == lhs & pe_miss$rhs == rhs]
lav_se_m  <- function(op, lhs, rhs) pe_miss$se[ pe_miss$op == op & pe_miss$lhs == lhs & pe_miss$rhs == rhs]


test_that("svcm and lavaan log-likelihoods agree with missing data", {
  expect_equal(
    as.numeric(logLik(fit_cfa_miss)),
    as.numeric(lavaan::logLik(fit_lav_miss)),
    tolerance = 1e-3
  )
})

test_that("svcm and lavaan factor loading estimates and SEs agree with missing data", {
  for (i in 2:I) {
    lname <- paste0("l", i)
    yname <- paste0("y", i)
    expect_equal(th_cfa_miss[[lname]], lav_est_m("=~", "eta", yname), tolerance = 1e-4,
                 label = paste(lname, "estimate"))
    expect_equal(se_cfa_miss[[lname]], lav_se_m( "=~", "eta", yname), tolerance = 1e-4,
                 label = paste(lname, "SE"))
  }
})

test_that("svcm and lavaan factor variance estimate and SE agree with missing data", {
  expect_equal(th_cfa_miss[["p"]], lav_est_m("~~", "eta", "eta"), tolerance = 1e-4)
  expect_equal(se_cfa_miss[["p"]], lav_se_m( "~~", "eta", "eta"), tolerance = 1e-4)
})

test_that("svcm and lavaan residual variance estimates and SEs agree with missing data", {
  for (i in 1:I) {
    tname <- paste0("th", i)
    yname <- paste0("y", i)
    expect_equal(th_cfa_miss[[tname]], lav_est_m("~~", yname, yname), tolerance = 1e-4,
                 label = paste(tname, "estimate"))
    expect_equal(se_cfa_miss[[tname]], lav_se_m( "~~", yname, yname), tolerance = 1e-4,
                 label = paste(tname, "SE"))
  }
})

test_that("svcm and lavaan intercept estimates and SEs agree with missing data", {
  for (i in 1:I) {
    uname <- paste0("u", i)
    yname <- paste0("y", i)
    expect_equal(th_cfa_miss[[uname]], lav_est_m("~1", yname, ""), tolerance = 1e-4,
                 label = paste(uname, "estimate"))
    expect_equal(se_cfa_miss[[uname]], lav_se_m( "~1", yname, ""), tolerance = 1e-4,
                 label = paste(uname, "SE"))
  }
})
