
# Equivalence with lme4 -------------------------------------------------------
#
# The random intercept model can be expressed either as:
#   - lme4:  y ~ 1 + (1 | subject)          (long format, REML = FALSE)
#   - svcm:  svc(P + TH, R = I), mc(U, X)   (wide format, Cholesky param)
#
# Both maximise the same multivariate normal likelihood, so logLik values
# should agree to floating-point precision and parameter estimates should
# match closely.

skip_if_not_installed("lme4")

# Simulate once and share across tests in this file
set.seed(6251)
I <- 4
J <- 300

eta <- rnorm(J)
eps <- matrix(rnorm(J * I), J, I)
Y   <- matrix(rep(eta, I), J, I) + eps

# lme4: long format
dat_long <- data.frame(
  y       = c(Y),
  subject = rep(seq_len(J), times = I)
)
fit_lme4 <- lme4::lmer(y ~ 1 + (1 | subject), data = dat_long, REML = FALSE)

# svcm: wide format, Cholesky parameterisation
R <- Matrix::Diagonal(J)
X <- matrix(1, J, 1)
mod_svcm <- svcm(
  Y,
  pm(I, 1, rep("lp",  I), TRUE,          1,          "LP"),
  pm(I, I, rep("lth", I), diag(TRUE, I), diag(1, I), "LTH"),
  pm(I, 1, "u",           TRUE,          0,           "U"),
  ic(tcrossprod(LP),  name = "P"),
  ic(tcrossprod(LTH), name = "TH"),
  svc(P + TH, R = R),
  mc(U, X = X)
)
fit_svcm <- fit_svcm(mod_svcm)


test_that("svcm and lme4 log-likelihoods agree", {
  ll_svcm <- as.numeric(logLik(fit_svcm))
  ll_lme4 <- as.numeric(logLik(fit_lme4))
  expect_equal(ll_svcm, ll_lme4, tolerance = 1e-4)
})

test_that("svcm and lme4 random intercept variance estimates agree", {
  # svcm: lp^2 = var_eta; lme4: VarCorr(fit)$subject[1]
  var_eta_svcm <- theta(fit_svcm)[["lp"]]^2
  var_eta_lme4 <- as.numeric(lme4::VarCorr(fit_lme4)$subject)
  expect_equal(var_eta_svcm, var_eta_lme4, tolerance = 1e-4)
})

test_that("svcm and lme4 residual variance estimates agree", {
  # svcm: lth^2 = var_eps; lme4: sigma^2
  var_eps_svcm <- theta(fit_svcm)[["lth"]]^2
  var_eps_lme4 <- sigma(fit_lme4)^2
  expect_equal(var_eps_svcm, var_eps_lme4, tolerance = 1e-4)
})

test_that("svcm and lme4 fixed intercept estimates agree", {
  intercept_svcm <- theta(fit_svcm)[["u"]]
  intercept_lme4 <- as.numeric(lme4::fixef(fit_lme4))
  expect_equal(intercept_svcm, intercept_lme4, tolerance = 1e-4)
})
