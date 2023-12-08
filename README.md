---
---README
---

**svcm** is an R-package intended for (relatively) easy specification of models combining structural equation models with sparse relationship matrices. The structural equations model involve unknown parameters that describe dependence across variables. The relationship matrix may for example encode genetic relatedness or different forms of shared environments and describe dependence across individuals. The interface is inspired by **OpenMx**, but uses the **Matrix** package for computation.

### An example

Consider a factor model for individual $j$

$$
\boldsymbol{y}_j = \boldsymbol{\beta}+\boldsymbol{\lambda}\eta_j+\boldsymbol{\delta}_j.
$$

$\boldsymbol{\beta}$ is a vector of intercepts, $\boldsymbol{\lambda}$ is a vector of loadings on the common factor $\eta_j$, and $\boldsymbol{\delta}_j$ a vector of specific factors with diagonal covariance matrix $\bf{\Theta}$. The factor model describes dependence across variables. The common factor is specified with a genetic and environmental component

$$
\eta_j = a_j + e_j.
$$

The genetic components $a_j$ represent a additive polygenic contribution to the common factor, and has covariance matrix $\sigma^2_a\bf{A}$. $\bf{A}$ is a sparse symmetric matrix describing genetic relatedness among individuals. When $\bf{A}$ cannot be arranged into small homogeneous block, this model cannot be expressed in common software for structural equation models. The environmental component is independent across individuals with covariance matrix $\sigma^2_e\bf{I}$. The model could also be written for the whole sample as

$$
\bf{Y} = \boldsymbol{x\beta}'+
\boldsymbol{\eta}\boldsymbol{\lambda}'+ \bf{\Delta}.
$$Across variables and individuals, the model implied covariance matrix can be described as

$$
(\boldsymbol{\lambda\lambda'}\sigma^2_a)\otimes\bf{A}+
(\boldsymbol{\lambda\lambda'}\sigma^2_e+\bf{\Theta})\otimes\bf{I}.
$$

### Simulate data

``` r
# Genetic relationships
N = 1000
R = rsparsematrix(N, N,
                  density = 0.01,
                  symmetric = T,
                  rand.x = function(n) sample(1 / c(1, 4, 8, 16), n, T))
diag(R) = 1
# Genetic values
a <- MASS::mvrnorm(1, rep(0, N), 2 * R)
# Environmental values
e <- rnorm(N, 0, sqrt(2))
# Common factors
eta <- a + e
# Loadings
l <- c(1, 0.5, 0.5, 0.8)
# Intercepts
b = c(2, 2, 4, 4)
# specific factors
del = matrix(rnorm(N * 4), N, 4)
# Data
x = rep(1, N)
Y = x %*% t(b) + eta %*% t(l) + del
```

### Fit the model
