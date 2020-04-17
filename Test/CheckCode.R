library(svcmr)
library(lavaan)

# Data
# --------------------------
N <- 100
W <- matrix(rnorm(N * 2), N, 2)
eta <- W %*% c(2, 3) + rnorm(N, 0, sqrt(2))
l <- matrix(c(1, 0.5, 0.5, 0.8), 4, 1)
b = c(2, 2, 4, 4)
Y <- matrix(NA, N, 4, dimnames = list(NULL, paste0("y", 1:4)))
for(i in 1:nrow(Y)) {
  Y[i, ] <- b + l %*% eta[i] + rnorm(4)
}
Renv <- Matrix::Diagonal(N)
X <- matrix(1, N, 1)
Y[1:10, 1] = NA

# svcmr
# --------------------------
# Failer hvis noen av de ikke-diagonale også heter "th". Hmm?
lab_th = matrix(NA, 4, 4)
#lab_th = matrix("th", 4, 4)
diag(lab_th) = c("th", "th", "th3", "th4")
mod <- svcm(pm(4, 1, paste0("l", 1:4), c(F, T, T, T), diag(1, 4), "L"),
            pm(4, 4, lab_th, diag(T, 4), diag(1, 4), "TH"),
            pm(1, 1, paste0("p", 1), T, 1, "P"),
            pm(4, 1, paste0("u", 1:4), T, 0, "U"),
            pm(1, 2, paste0("w", 1:2), T, 0, "G"),
            svc(L %*% P %*% t(L) + TH, R = Renv),
            mc(U, X = X),
            mc(L %*% G, X = W))
fit <- fitm(Y, mod, se = F, control = list(trace = 6))
summary(fit)
#fit$hessian <- compHess(fit$fit_objective, fit$fit$par)


lab_free = c(paste0("l", 1:4)[c(F, T, T, T)],
                lab_th[diag(T, 4)],
                paste0("p", 1)[T],
                paste0("u", 1:4)[T],
                paste0("w", 1:2)[T])
lab_values = lab_free[!duplicated(lab_free)]
values = runif(length(lab_values))
names(values) = lab_values
mo = pm(4, 4, lab_th, diag(T, 4), diag(1, 4), "TH")
mo = pm(4, 1, paste0("l", 1:4), c(F, T, T, T), diag(1, 4), "L")

# Lavaan
# --------------------------
dat = data.frame(Y, W)
modl = "
eta =~ 1*y1 + l2*y2 + l3*y3 + l4*y4
eta ~ w1*X1 + w2*X2
y1~~th*y1
y2~~th*y2
y3~~th3*y3
y4~~th4*y4
"
fitl = cfa(modl, dat, meanstructure = T,  missing = "ML")
summary(fitl)

logLik(fit)
logLik(fitl)

# Free
# --------------------------
mod2 <- svcm(pm(4, 1, paste0("l", 1:4), c(F, T, T, T), diag(1, 4), "L"),
            pm(4, 4, lab_th, diag(T, 4), diag(1, 4), "TH"),
            pm(1, 1, paste0("p", 1), T, 1, "P"),
            pm(4, 1, paste0("u", 1:4), T, 0, "U"),
            pm(1, 2, paste0("w", 1:2), T, 0, "G"),
            svc((L %*% P %*% t(L) + TH) %x% Renv),
            mc(U, X = X),
            mc(L %*% G, X = W))
fit2 <- fitm(Y, mod2, se = F, control = list(trace = 6))
summary(fit2)
