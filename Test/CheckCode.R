library(svcmr)
library(lavaan)

# Data
# --------------------------
N <- 100
# Common factor
W <- matrix(rnorm(N * 2), N, 2)
eta <- W %*% c(2, 3) + rnorm(N, 0, sqrt(2))
# Loadings
l <- matrix(c(1, 0.5, 0.5, 0.8), 4, 1)
# Intercepts
b = c(2, 2, 4, 4)
# Data
Y <- matrix(NA, N, 4, dimnames = list(NULL, paste0("y", 1:4)))
for(i in 1:nrow(Y)) {
  Y[i, ] <- b + l %*% eta[i] + rnorm(4)
}
# Relationship matrix for unique environmental deviations
Renv <- Matrix::Diagonal(N)
# Covariates
X <- matrix(1, N, 1)

Y[1:10, 1] = NA
# svcmr
# --------------------------
mod <- svcm(# Parameters
  pm(nrow = 4, ncol = 1, labels = paste0("l", 1:4), free = c(F, T, T, T), values = diag(1, 4), name = "L"),
  pm(nrow = 4, ncol = 4, labels = paste0("th", 1:16), free = diag(T, 4), values = diag(1, 4), name = "TH"),
  pm(nrow = 1, ncol = 1, labels = paste0("p", 1), free = T, values = 1, name = "P"),
  pm(nrow = 4, ncol = 1, labels = paste0("u", 1:4), free = T, values = 0, name = "U"),
  pm(nrow = 1, ncol = 2, labels = paste0("w", 1:2), free = T, values = 0, name = "G"),
  # Variance components
  svc(L %*% P %*% t(L) + TH, R = Renv),
  # Mean components
  mc(U, X = X),
  mc(L %*% G, X = W))
fit <- fitm(Y, mod, se = T, control = list(trace = 6))
summary(fit)

# Lavaan
# --------------------------
dat = data.frame(Y, W)
modl = "
eta =~ 1*y1 + l2*y2 + l3*y3 + l4*y4
eta ~ w1*X1 + w2*X2
"
fitl = cfa(modl, dat, meanstructure = T,  missing = "ML")
summary(fitl)

logLik(fit)
logLik(fitl)
