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
datl = data.frame(y = c(Y), j = rep(1:J, times = I), k = rep(1:K, each = K * J))

# lmer fit
fitlme <- lmer(y~1+(1|j)+(1|j:k), datl, REML = F)
fitlme4 = lmer(y~1+(1|j), datl, REML = F)
anova(fitlme4, fitlme)
summary(fitlme)

# svcmr
R = Matrix::Diagonal(J)
X = matrix(1, J, 1)
datw = data.frame(Y)
mod <- svcm(Y,
            pm(4, 4, "th", diag(T, 4), diag(1, 4), "TH"),
            pm(4, 4, "p", T, 1, "P"),
            pm(2, 2, "l", T, 1, "L"),
            pm(4, 1, "u", T, 0, "U"),
            ic(diag(2) %x%  L, "LL"),
            svc((P + LL + TH), R),
            mc(U, X = X))
mod <- fit_svcm(mod, se = T, control = list(trace = 1))
summary(mod)

mod4 <- svcm(Y,
            pm(4, 4, "th", diag(T, 4), diag(1, 4), "TH"),
            pm(4, 4, "p", T, 1, "P"),
            pm(4, 1, "u", T, 0, "U"),
            svc((P + TH), R),
            mc(U, X = X))
mod4 <- fit_svcm(mod4, se = F, control = list(trace = 1))

mod6 <- svcm(Y,
             pm(4, 4, "th", diag(T, 4), diag(1, 4), "TH"),
             pm(4, 1, "u", T, 0, "U"),
             svc(TH, R),
             mc(U, X = X))
mod6 <- fit_svcm(mod6, se = F, control = list(trace = 1))

anova(mod, mod4)
anova(mod4, mod, mod6)
theta(mod)
mod$opt$par
compute(mod, P)

-0.5 * objective(mod)
logLik(mod)
-0.5 * mod$opt$objective
logLik(fitlme)
compute(mod, LL + P + TH)

trans = function(x) {
  x[3] / 4
}

Jz = jacobian(trans, theta(mod))
VC2 = Jz %*% vcov(mod) %*% t(Jz)
round(vcov(mod), 5)
sqrt(diag(VC2))


# MIMIC Model
# --------------------------
W = matrix(NA, J, 2, dimnames = list(NULL, paste0("w", 1:2)))
Y = matrix(NA, J, I, dimnames = list(NULL, paste0("y", 1:I)))
for(j in 1:J) {
  w = rnorm(2)
  eta = sum(c(2, 1.5) * w) + rnorm(1)
  eps = rnorm(4)
  Y[j, ] = c(1, 0.5, 0.5, 0.7) * eta + eps
  W[j, ] = w
}
datw = data.frame(Y, W)

# Lavaan
modl = "
eta =~ 1*y1 + l2*y2 + l3*y3 + l4*y4
eta ~ b1*w1 + b2*w2
"
fitlav = cfa(modl, datw, meanstructure = T)

# svcmr
mod2 = svcm(Y,
            pm(4, 1, paste0("l", 1:4), c(F, T, T, T), 1, "L"),
            pm(1, 1, "p", T, 1, "P"),
            pm(4, 4, paste0("th", 1:4), diag(T, 4), diag(1, 4), "TH"),
            pm(4, 1, paste0("u", 1:4), T, 0, "U"),
            pm(1, 2, paste0("b", 1:2), T, 0, "B"),
            svc(L %*% P %*% t(L) + TH, R = R),
            mc(U, X = X),
            mc(L %*% B, X = W))
fit2 = fitm(mod2, se = T, control = list(trace = 1))
summary(fit2)

logLik(fit2)
logLik(fitlav)

compute(fit2, B)
trans2 = function(x) {
  pos_b = which(grepl("b", names(theta(mod2))))
  b = x[pos_b]
  bmat = matrix(b, 1, 2)
  bmat %*% (4 * diag(2))
}

trans2(fit2$opt$par)
jz2 = jacobian(trans2, fit2$opt$par)
jz2 %*% vcov(fit2) %*% t(jz2)
# Missing
# --------------------------
Y[1:10, 1] = NA
datwmiss = data.frame(Y, W)

# Lavaan
fitlavmiss = cfa(modl, datwmiss, meanstructure = T, missing = "ML")

# svcmr
mod2miss = svcm(Y,
                pm(4, 1, paste0("l", 1:4), c(F, T, T, T), 1, "L"),
                pm(1, 1, "p", T, 1, "P"),
                pm(4, 4, paste0("th", 1:4), diag(T, 4), diag(1, 4), "TH"),
                pm(4, 1, paste0("u", 1:4), T, 0, "U"),
                pm(1, 2, paste0("b", 1:2), T, 0, "B"),
                svc(L %*% P %*% t(L) + TH, R = R),
                mc(U, X = X),
                mc(L %*% B, X = W))
mod2miss = fitm(mod2miss, se = T, control = list(trace = 1))

logLik(mod2miss)
logLik(fitlavmiss)

anova(mod, mod2miss)

# Delta
# --------------------------
mod2missL = svcm(Y,
                pm(4, 1, paste0("l", 1:4), c(F, T, T, T), 1, "L"),
                pm(1, 1, "lp", T, 1, "LP"),
                ic(LP %*% t(LP), name = "P"),
                pm(4, 4, paste0("th", 1:4), diag(T, 4), diag(1, 4), "TH"),
                pm(4, 1, paste0("u", 1:4), T, 0, "U"),
                pm(1, 2, paste0("b", 1:2), T, 0, "B"),
                svc(L %*% P %*% t(L) + TH, R = R),
                mc(U, X = X),
                mc(L %*% B, X = W))
mod2missL = fitm(mod2missL, se = T)


tra = function(theta) {
  #theta[4] = theta[4]^2
  #theta
  #l = theta[1:3]
  #c = l %*% t(l)
  #c[lower.tri(c)]
  theta[4]^2
}

J = numDeriv::jacobian(tra, mod2missL$opt$par)
C = solve(0.5 * mod2missL$H)
VC = t(J) %*% C %*% J
sqrt(diag(VC))
sqrt(J %*% C %*% t(J))
