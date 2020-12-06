library(svcmr)
library(lavaan)
library(lme4)

# Random intercept model
# --------------------------
I = 4
J = 10000
Y = matrix(NA, J, I, dimnames = list(NULL, paste0("y", 1:I)))
for(j in 1:J) {
  eta = rnorm(1)
  eps = rnorm(4)
  Y[j, ] = eta + eps
}
R = Matrix::Diagonal(J)
X= matrix(1, J, 1)
datw = data.frame(Y)
datl = data.frame(y = c(Y), j = rep(1:J, times = I))

# lmer fit
fitlme = lmer(y~1+(1|j), datl, REML = F)
summary(fitlme)

# svcmr
mod <- svcm(pm(4, 4, "th", diag(T, 4), diag(1, 4), "TH"),
            pm(4, 4, "p", T, 1, "P"),
            pm(4, 1, "u", T, 0, "U"),
            svc(P + TH, R = R),
            mc(U, X = X))
fit = fitm(Y, mod, se = F, control = list(trace = 1))
summary(fit)

logLik(fit)
logLik(fitlme)

# MIMIC Model
# --------------------------
W = matrix(NA, J, 2, dimnames = list(NULL, paste0("w", 1:2)))
Y = matrix(NA, J, I, dimnames = list(NULL, paste0("y", 1:I)))
for(j in 1:J) {
  w = rnorm(2)
  eta = sum(c(2, 1.5) * w) + rnorm(1)
  eps = rnorm(4)
  Y[j, ] = c(1,0.5, 0.5, 0.7) *  eta + eps
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
mod2 = svcm(pm(4, 1, paste0("l", 1:4), c(F, T, T, T), 1, "L"),
             pm(1, 1, "p", T, 1, "P"),
             pm(4, 4, paste0("th", 1:4), diag(T, 4), diag(1, 4), "TH"),
             pm(4, 1, paste0("u", 1:4), T, 0, "U"),
             pm(1, 2, paste0("b", 1:2), T, 0, "B"),
             svc(L %*% P %*% t(L) + TH, R = R),
             mc(U, X = X),
             mc(L %*% B, X = W))
fit2 = fitm(Y, mod2, se = F, control = list(trace = 1))
summary(fit2)

logLik(fit2)
logLik(fitlav)

# Missing
# --------------------------
Y[1:10, 1] = NA
datw = data.frame(Y, W)

# Lavaan
fitlavmiss = cfa(modl, datw, meanstructure = T, missing = "ML")

# svcmr
fit2miss = fitm(Y, mod2, se = F, control = list(trace = 1))

logLik(fit2miss)
logLik(fitlavmiss)
