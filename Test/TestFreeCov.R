library(svcmr)
library(OpenMx)

# Genetic relationships
str1 <- matrix(c(1, 1/4, 1/4, 1/4, 1, 1/2, 1/4, 1/2, 1), 3, 3)
str2 <- matrix(c(1, 1/2, 1/8, 1/8, 1/2, 1, 1/8, 1/8, 1/8, 1/8, 1, 1/2, 1/8, 1/8, 1/2, 1), 4, 4)
# Genetic relationship matrix
Rgen <- Matrix::bdiag(Matrix::bdiag(replicate(100, str1, simplify = FALSE)),
                      Matrix::bdiag(replicate(100, str2, simplify = FALSE)))

# Number of observations
N <- dim(Rgen)[1]
# Sample genotypic values
a <- MASS::mvrnorm(1, rep(0, N * 2), matrix(c(3, 1, 1, 2), 2, 2) %x% Rgen)
a1 <- a[1:N]
a2 <- a[(N + 1):(2 * N)]
# Sample environmental deviations
e1 <- rnorm(N, 0, sqrt(2))
e2 <- rnorm(N, 0, sqrt(3))
# Intercepts
b = c(2, 4)
# Data
Y <- matrix(NA, N, 2, dimnames = list(NULL, paste0("y", 1:2)))
for(i in 1:nrow(Y)) {
  Y[i, 1] <- b[1] + a1[i] + e1[i]
  Y[i, 2] <- b[2] + a2[i] + e2[i]
}

# Relationship matrix for unique environmental deviations
Renv <- Matrix::Diagonal(dim(Rgen)[1])

# Covariates
X <- matrix(1, N, 1)

mod <- svcm(# Parameters
  pm(nrow = 2, ncol = 2, labels = c("va11", "va21", "va21", "va22"), free = T, values = c(2, 1, 1, 2), name = "va"),
  pm(nrow = 2, ncol = 2, labels = c("ve11", "ve21", "ve21", "ve22"), free = T, values = c(2, 1, 1, 2), name = "ve"),
  pm(nrow = 2, ncol = 1, labels = paste0("u", 1:2), free = T, values = 0, name = "U"),
  # Variance components
  svc(va, R = Rgen),
  svc(ve, R = Renv),
  # Mean components
  mc(U, X = X))
fit <- fitm(Y, mod, se = TRUE, control = list(trace = 6))
summary(fit)

st = proc.time()
mod2 <- svcm(pm(2, 2, c("va11", "va21", "va21", "va22"), T, c(2, 1, 1, 2), "va"),
             pm(2, 2, c("ve11", "ve21", "ve21", "ve22"), T, c(2, 1, 1, 2), "ve"),
             pm(2, 1, paste0("u", 1:2), T, 0, "U"),
             svc(va %x% Rgen + ve %x% Renv),
             mc(U, X))
fit2 <- fitm(Y, mod2, se = TRUE, control = list(trace = 6))
en = proc.time() - st
summary(fit2)

mod3 <- svcm(# Parameters
  pm(nrow = 2, ncol = 2, labels = c("va11", "va21", "va21", "va22"), free = T, values = c(2, 1, 1, 2), name = "va"),
  pm(nrow = 2, ncol = 2, labels = c("ve11", "ve21", "ve21", "ve22"), free = T, values = c(2, 1, 1, 2), name = "ve"),
  pm(nrow = 2, ncol = 1, labels = paste0("u", 1:2), free = T, values = 0, name = "U"),
  # Variance components
  svc(rbind(cbind(va[1, 1] * Rgen, va[1, 2] * Rgen),
            cbind(va[2, 1] * Rgen, va[2, 2] * Rgen))),
  svc(rbind(cbind(ve[1, 1] * Renv, ve[1, 2] * Renv),
            cbind(ve[2, 1] * Renv, ve[2, 2] * Renv))),
  # Mean components
  mc(U, X = X))
fit3 <- fitm(c(Y), mod3, se = TRUE, control = list(trace = 6))
summary(fit3)

anova(fit, fit2, fit3)

# OpenMx
# -------------------------------------------
Xd = diag(2) %x% X
datOpenMx = data.frame(t(c(Y)))
mod2OpenMx <- mxModel(mxMatrix("Fu", 2, 2, T, c(2, 1, 1, 2), c("va11", "va21", "va21", "va22"), name = "va"),
                      mxMatrix("Fu", 2, 2, T, c(2, 1, 1, 2), c("ve11", "ve21", "ve21", "ve22"), name = "ve"),
                      mxMatrix("Fu", 1, 2 * N, T, 0, c(rep("u1", N), rep("u2", N)), name = "U"),
                      mxMatrix("Fu", N, N, F, as.matrix(Rgen), name = "Rgen"),
                      mxMatrix("Fu", N, N, F, as.matrix(Renv), name = "Renv"),
                      mxMatrix("Fu", 2, 2 * N, F, t(Xd), name = "Xt"),
                      mxMatrix("Fu", 1, 2, T, 0, c("u1", "u2"), name = "U"),
                      mxAlgebra(va %x% Rgen + ve %x% Renv, name = "C"),
                      mxAlgebra(U %*% Xt, name = "M"),
                      mxExpectationNormal("C", "M", colnames(datOpenMx)),
                      mxFitFunctionML(),
                      mxData(datOpenMx, "raw"))
fit2openMx <- mxRun(mod2OpenMx)
summary(fit2openMx)
