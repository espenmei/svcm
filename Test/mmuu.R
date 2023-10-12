library(svcm)
library(lavaan)
library(lme4)


# Random intercept model
# --------------------------
I = 4
J = 1000
K = 2
Y = matrix(NA, J, I, dimnames = list(NULL, paste0("y", 1:I)))
for(j in 1:J) {
  eta_j = rnorm(1)
  eta_k = rnorm(2)
  eps = rnorm(4)
  Y[j, ] = rep(eta_j, 4) + rep(eta_k, each = 2) + eps
}
Y[1:10, 1] = NA
X = cbind(1, matrix(runif(J * 2), J, 2))
datl = data.frame(y = c(Y),
                  rbind(X, X, X, X),
                  vr = rep(1:4, each = J),
                  j = rep(1:J, times = I),
                  k = rep(1:K, each = K * J))

# lmer fit
fitlme <- lmer(y~0+(X1 + X2 + X3) : factor(vr) + (1|j)+(1|j:k), datl, REML = F)
fitlme4 = lmer(y~1+(1|j), datl, REML = F)
anova(fitlme4, fitlme)
summary(fitlme)

# svcmr
R = Matrix::Diagonal(J)
datw = data.frame(Y)
mod <- svcm(Y,
            pm(4, 4, "th", diag(T, 4), diag(1, 4), "TH"),
            pm(4, 4, "p", T, 1, "P"),
            pm(2, 2, "l", T, 1, "L"),
            pm(4, 3, paste0("u", 1:12), T, 1:12, "U"),
            ic(diag(2) %x%  L, "LL"),
            svc((P + LL + TH), R),
            mc(U, X = X))
.compute(mod$mcs[[1]], mod$env_comp)
svcm:::.compute(mod$mcs[[1]], mod$env_comp)
svcm:::.compute.fixedmc(mod$mcs[[1]], mod$env_comp)
expected_mean(mod)
bench::bench_time({mod = fit_svcm(mod, se = F, control = list(trace = 1))})
summary(mod)

theta(mod)
mod$opt$par
compute(mod, P)
bench::bench_time(mod$dat$y[mod$dat$keepy])
bench::bench_time(expected_mean(mod))
bench::bench_time(expected_cov(mod))

-0.5 * objective(mod)
logLik(mod)
-0.5 * mod$opt$objective
logLik(fitlme)
compute(mod, LL + P + TH)
