---
editor: visual
title: README
toc-title: Table of contents
---

## svcm

**svcm** is an R-package intended for (relatively) easy specification of
models combining structural equation models with sparse relationship
matrices. The structural equations model involve unknown parameters that
describe dependence across variables. The relationship matrix may for
example encode genetic relatedness or different forms of shared
environments and describe dependence across individuals. The interface
is inspired by **OpenMx**, but uses the **Matrix** package for
computation.

### An example

Consider a factor model for individual $j$

$$
\boldsymbol{y}_j = \boldsymbol{\beta}+\boldsymbol{\lambda}\eta_j+\boldsymbol{\delta}_j.
$$

$\boldsymbol{\beta}$ is a vector of intercepts, $\boldsymbol{\lambda}$
is a vector of loadings on the common factor $\eta_j$, and
$\boldsymbol{\delta}_j$ a vector of specific factors with diagonal
covariance matrix $\bf{\Theta}$. The factor model describes dependence
across variables. The common factor is specified with a genetic and
environmental component $$
\eta_j = a_j + e_j.
$$ The genetic components $a_j$ represent a additive polygenic
contribution to the common factor, and has covariance matrix
$\sigma^2_a\bf{A}$. $\bf{A}$ is a sparse symmetric matrix describing
genetic relatedness among individuals. When $\bf{A}$ cannot be arranged
into small homogeneous block, this model cannot be expressed in common
software for structural equation models. The environmental component is
independent across individuals with covariance matrix
$\sigma^2_e\bf{I}$. The model could also be written for the whole sample
as

$$
\bf{Y} = \boldsymbol{x+\beta}'+
\boldsymbol{\eta}\boldsymbol{\lambda}'+ \bf{\Delta}.
$$ Across variables and individuals, the model implied covariance matrix
can be described as

$$
(\boldsymbol{\lambda\lambda'}\sigma^2_a)\otimes\bf{A}+
(\boldsymbol{\lambda\lambda'}\sigma^2_e+\bf{\Theta})\otimes\bf{I}.
$$

## Simulate data

::: cell
``` {.r .cell-code}
# Genetic relationships
str1 <- matrix(c(1, 1/4, 1/4, 1/4, 1, 1/2, 1/4, 1/2, 1), 3, 3)
str2 <- matrix(c(1, 1/2, 1/8, 1/8, 1/2, 1, 1/8, 1/8, 1/8, 1/8, 1, 1/2, 1/8, 1/8, 1/2, 1), 4, 4)
# Genetic relationship matrix
A <- Matrix::bdiag(Matrix::bdiag(replicate(100, str1, simplify = FALSE)),
                   Matrix::bdiag(replicate(100, str2, simplify = FALSE)))
L <- chol(A)
N <- nrow(A)
# Genetic values
a <- sqrt(2) * t(L) %*% rnorm(N)
# Environmental values
e <- rnorm(N, 0, sqrt(2))
# Common factors
eta <- a + e
# Loadings
l <- c(1, 0.5, 0.5, 0.8)
# Intercepts
b = c(2, 2, 4, 4)
# specific factors
del <- matrix(rnorm(N * 4), N, 4)
# Data
x <- rep(1, N)
Y <- x %*% t(b) + eta %*% t(l) + del
```
:::

## Fit the model

Below is an example of one way to fit this model using **svcm**. `pm()`
is used to define parameters as matrices. `svc()` is used to define
structured variance components. Notice that the expressions use
parameter matrices to define the model implied covariance structure
across variables. The argument `R =` is used to supply a relationship
matrix describing the dependence across individuals. Notice also that
`mc()` is used to define the mean structure. There is no limit to how
many `svc()` and `mc()` objects that can be supplied, they will be added
together to form the total mean and covariance structure.

::: cell
``` {.r .cell-code}
library(svcm)
# Relationship matrix for environmental deviations
I <- Matrix::Diagonal(N)
m1 <- svcm(
  #Data
  Y,
  # Parameters
  pm(nrow = 4, ncol = 1, labels = paste0("l", 1:4), free = c(F, T, T, T), values = c(1, 0.5, 0.5, 0.5), name = "l"),
  pm(nrow = 1, ncol = 1, labels = "Sa1", free = T, values = 1, name = "Sa"),
  pm(nrow = 1, ncol = 1, labels = "Se1", free = T, values = 1, name = "Se"),
  pm(nrow = 4, ncol = 4, labels = sapply(1:4, function(x) paste0("th", 1:4, x)), free = diag(T, 4), values = diag(4), name = "TH"),
  pm(nrow = 4, ncol = 1, labels = paste0("b", 1:4), free = T, values = 0, name = "b"),
  # Variance components
  svc(l %*% Sa %*% t(l), R = A),
  svc(l %*% Se %*% t(l) + TH, R = I),
  # Mean components
  mc(b, X = x)
)
```
:::

After the model is defined, it can be fitted using `fit_svcm()`.

::: cell
``` {.r .cell-code}
m1 <- fit_svcm(m1, se = TRUE, control = list(trace = 10))
```

::: {.cell-output .cell-output-stdout}

    iter: objective:    l2  l3  l4  Sa1 Se1 th11    th22    th33    th44    b1  b2  b3  b4  
      0:     21954.323: 0.500000 0.500000 0.500000  1.00000  1.00000  1.00000  1.00000  1.00000  1.00000  0.00000  0.00000  0.00000  0.00000
     10:     10587.914: 0.930493  1.11513  1.44804  2.14166  1.43436  1.88293  1.07970  1.14968 0.845133 0.00712701 0.204433  1.58118  1.01976
     20:     9490.3724: 0.500211 0.511787 0.825085  2.82894  1.25056  1.28862 0.938962 0.963374 0.806252  2.05039  1.96042  4.00035  4.05661
     30:     9485.3146: 0.509562 0.516190 0.807101  2.79282  1.39661  1.22944 0.940035 0.983219 0.838223  2.12898  1.98763  4.06667  4.06862
     40:     9484.1170: 0.504963 0.510909 0.804589  2.59758  1.55659  1.23057 0.940176 0.991343 0.834149  2.12677  1.99732  4.07754  4.06784
     50:     9482.5998: 0.506133 0.515189 0.806423  2.28824  1.81498  1.22994 0.944333 0.982058 0.835053  2.12808  1.99310  4.06955  4.06884
     60:     9481.8772: 0.506071 0.514296 0.805584  1.97229  2.10655  1.21986 0.942956 0.986212 0.837098  2.13376  1.99829  4.07241  4.07446
     70:     9481.8117: 0.504671 0.516250 0.806331  1.86204  2.19059  1.22719 0.941902 0.981451 0.835377  2.13053  1.99424  4.06933  4.07016
     80:     9481.8000: 0.506014 0.515066 0.806494  1.84607  2.20805  1.22751 0.942163 0.981533 0.836768  2.12938  1.99418  4.06973  4.06906
:::

::: {.cell-output .cell-output-stderr}
    Computing standard errors.
:::
:::

After the models is fitted, it can be inspected with `summary()`.

::: cell
``` {.r .cell-code}
summary(m1)
```

::: {.cell-output .cell-output-stdout}
    Log likelihood       Deviance            AIC            BIC 
         -4740.900       9481.800       9507.800       9584.986 

    Fitted parameters:
         Estimate Std. Error z value Pr(>|z|)
    l2     0.5060     0.0232   21.78   <0.001
    l3     0.5151     0.0236   21.79   <0.001
    l4     0.8065     0.0288   27.99   <0.001
    Sa1    1.8462     0.5067    3.64   <0.001
    Se1    2.2078     0.4699    4.70   <0.001
    th11   1.2275     0.1133   10.84   <0.001
    th22   0.9422     0.0581   16.23   <0.001
    th33   0.9815     0.0604   16.26   <0.001
    th44   0.8367     0.0750   11.16   <0.001
    b1     2.1293     0.0971   21.94   <0.001
    b2     1.9942     0.0575   34.67   <0.001
    b3     4.0697     0.0586   69.43   <0.001
    b4     4.0690     0.0786   51.75   <0.001
:::
:::
