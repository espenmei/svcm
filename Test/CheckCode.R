library(svcmr)
library(lavaan)

N = 100
P = 4
Y = matrix(rnorm(N * P), N, P)
X = matrix(1, N, 1)
W = matrix(rnorm(N * 2), N, 2)
R = Diagonal(N)

# Common factors
l = pm(4, 1, paste0("l", 1:4), c(F, T, T, T), 1, "l")
c = pm(1, 1, "c", T, 1, "c")
# Specific factors
s = pm(4, 4, paste0("es", 1:16), free = diag(T, 4), values = c(diag(2, 4)), name = "s")
# VCs
svc1 = svc(l %*% (c %*% t(c)) %*% t(l) + s %*% t(s), R = R)
# Means
B = pm(4, 1, labels = paste0("b", 1:4), free = T, values = colMeans(Y), name = "B")
mc_X = mc(B, X = X)
A = pm(1, 2, labels = paste0("a", 1:2), free = T, values = 0, name = "A")
mc_W = mc(l %*% A, X = W)
# Combine
mod1 = svcm(Y, l, c, s, svc1, B, mc_X, A, mc_W)
fit1 = fitm(mod1, se = T, control = list(trace = 6))
summary(fit)

# Compute
# --------------------------
compute(mod1, l %*% (c %*% t(c)) %*% t(l))
compute(fit1, l %*% (c %*% t(c)) %*% t(l))

# Comparison lavaan
# --------------------------
colnames(Y) = c("v1", "v2", "v3", "v4")
colnames(W) = c("w1", "w2")
dat_lav = data.frame(Y, W)
mod_lav = "
c =~ v1 + v2 + v3 + v4
c~ a1*w1 + a2*w2
"
fit_lav = cfa(mod_lav, dat = dat_lav, meanstructure = T)
summary(fit_lav)

