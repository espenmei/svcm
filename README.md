svcm
================

**svcm** is an R-package intended for (relatively) easy specification of
models combining structural equation models with sparse relationship
matrices. The structural equations model involve unknown parameters that
describe dependence across variables. The relationship matrix may for
example encode genetic relatedness or different forms of shared
environments and describe dependence across individuals. The interface
is inspired by **OpenMx**, but uses the **Matrix** package for
computation.

## An example

Consider a factor model for individual $j$

$$
\boldsymbol{y}_j = \boldsymbol{\beta}+\boldsymbol{\lambda}\eta_j+\boldsymbol{\delta}_j.
$$

$\boldsymbol{\beta}$ is a vector of intercepts, $\boldsymbol{\lambda}$
is a vector of loadings on the common factor $\eta_j$, and
$\boldsymbol{\delta}_j$ a vector of specific factors with diagonal
covariance matrix $\mathbf{\Theta}$. The factor model describes
dependence across variables. The common factor is specified with a
genetic and environmental component

$$
\eta_j = a_j + e_j.
$$

The genetic components $a_j$ represent an additive polygenic
contribution to the common factor and has covariance matrix
$\sigma^2_a\bf{A}$. $\bf{A}$ is a sparse symmetric matrix describing
genetic relatedness among individuals. When $\bf{A}$ cannot be arranged
into small homogeneous blocks, this model cannot be expressed in common
software for structural equation models. The environmental component is
independent across individuals with covariance matrix
$\sigma^2_e\bf{I}$. The model could also be written for the whole sample
as

$$
\bf{Y} = \boldsymbol{x\beta}'+
\boldsymbol{\eta}\boldsymbol{\lambda}'+ \bf{\Delta}.
$$

Across variables and individuals, the covariance matrix implied by this
model can be described as $$
(\boldsymbol{\lambda\lambda'}\sigma^2_a)\otimes\bf{A}+
(\boldsymbol{\lambda\lambda'}\sigma^2_e+\bf{\Theta})\otimes\bf{I}.
$$

### Simulate data

Below is a simulation of data from the model.

``` r
# Genetic relationships
str1 <- matrix(c(1, 1/4, 1/4, 1/4, 1, 1/2, 1/4, 1/2, 1), 3, 3)
str2 <- matrix(c(1, 1/2, 1/8, 1/8, 1/2, 1, 1/8, 1/8, 1/8, 1/8, 1, 1/2, 1/8, 1/8, 1/2, 1), 4, 4)
# Genetic relationship matrix
A <- Matrix::bdiag(Matrix::bdiag(replicate(100, str1, simplify = FALSE)),
                   Matrix::bdiag(replicate(100, str2, simplify = FALSE)))
L <- chol(A)
N <- nrow(A)
a <- sqrt(2) * t(L) %*% rnorm(N)
e <- rnorm(N, 0, sqrt(2))
eta <- a + e
l <- c(1, 0.5, 0.5, 0.8)
b = c(2, 2, 4, 4)
del <- matrix(rnorm(N * 4), N, 4)
x <- rep(1, N)
Y <- x %*% t(b) + eta %*% t(l) + del
```

### Fit the model

There are several ways to fit this model using **svcm**, but it will
likely be somewhat similar to the example below. `pm()` is used to
define matrices of parameters, $\boldsymbol{\lambda}$,
$\boldsymbol{\beta}$, $\sigma^2_a$, $\sigma^2_e$, and $\bf \Theta$.
`svc()` is used to define structured variance components, the terms
defining the covariance matrix. Notice that the first argument is an
expression defining the model implied covariance structure across
variables as a function of parameters. The second argument `R =` is used
to supply a relationship matrix describing the dependence across
individuals. `mc()` is used to define the mean structure. The first
argument is an expression involving parameters and the second argument
`X =` gives a model matrix of covariates. There is no limit to how many
`pm`, `svc` and `mc` objects that can be supplied, they will be added
together to form the total mean and covariance structure.

``` r
library(svcm)
# Environmental relationship matrix
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

After the model is defined, it can be fitted using `fit_svcm()`.

``` r
m1 <- fit_svcm(m1, se = TRUE, control = list(trace = 10))
```


    iter: objective:    l2  l3  l4  Sa1 Se1 th11    th22    th33    th44    b1  b2  b3  b4  
      0:     22379.702: 0.500000 0.500000 0.500000  1.00000  1.00000  1.00000  1.00000  1.00000  1.00000  0.00000  0.00000  0.00000  0.00000
     10:     10210.760: 0.795091 0.881317  1.34942  2.29318  1.25989  1.73976  1.01528  1.09899  1.04357 0.247248 0.703761  2.36945  1.71903
     20:     9481.3160: 0.466511 0.565542 0.814708  2.71064  1.42126 0.990330 0.930488 0.974211 0.964322  1.86920  2.04426  3.97195  3.87796
     30:     9464.2776: 0.500275 0.527344 0.823104  2.64707  1.37795  1.05123 0.889825 0.985402  1.01873  2.01884  2.06592  4.06584  4.04830
     40:     9463.4903: 0.501289 0.526475 0.825544  2.68233  1.34173  1.04999 0.881956 0.986979  1.00640  2.08199  2.09473  4.09840  4.09535
     50:     9463.3519: 0.500755 0.529274 0.824981  2.74268  1.29960  1.05438 0.886651 0.985399  1.01059  2.09701  2.10172  4.10675  4.10861
     60:     9463.3308: 0.500697 0.529418 0.824869  2.77324  1.27485  1.05401 0.886323 0.985679  1.01208  2.10311  2.10504  4.10964  4.11386
     70:     9463.3252: 0.500276 0.529175 0.824742  2.80185  1.25388  1.05376 0.885822 0.985003  1.01201  2.10654  2.10699  4.11187  4.11695

    Computing standard errors.

After the models is fitted, it can be inspected with `summary()`.

``` r
summary(m1)
```

    Log likelihood       Deviance            AIC            BIC 
         -4731.662       9463.325       9489.325       9566.511 

    Fitted parameters:
         Estimate Std. Error z value Pr(>|z|)
    l2     0.5004     0.0223   22.42   <0.001
    l3     0.5292     0.0236   22.45   <0.001
    l4     0.8247     0.0287   28.70   <0.001
    Sa1    2.8035     0.5233    5.36   <0.001
    Se1    1.2527     0.4306    2.91   0.0036
    th11   1.0534     0.1032   10.20   <0.001
    th22   0.8863     0.0548   16.17   <0.001
    th33   0.9851     0.0610   16.14   <0.001
    th44   1.0120     0.0807   12.54   <0.001
    b1     2.1074     0.1007   20.93   <0.001
    b2     2.1074     0.0586   35.99   <0.001
    b3     4.1122     0.0619   66.48   <0.001
    b4     4.1176     0.0855   48.13   <0.001

## In general

**svcm** can fit models where the covariance matrix can be described as

$$
\bf{V}=\sum_{i=1}^{q}\mathbf \Sigma_i \otimes \mathbf R_i.
$$

The terms in the sum are specified by adding objects returned from
`svc()` to the model. $\mathbf \Sigma_i$ has dimension according to the
number of columns of $\bf{Y}$ and describes dependence among variables.
$\mathbf \Sigma_i$ is a function of parameters defined with `pm()` and
can be specified flexibly using general functions available in **R**.
$\mathbf R_i$ describe dependence across individuals and must be sparse
symmetric matrices with known values.

The means are described with a similar structure

$$
\bf M =\sum_{i=1}^{p} \mathbf X_i \mathbf B_i^{'}.
$$

For the means, the terms in the sum are specified by adding objects
returned from `mc()` to the model. $\mathbf B_i^{'}$ has rows equal to
the number of columns of $\bf Y$, and columns equal to the columns of
$\mathbf X_i$.

The data is described with the normal distribution $$
vec(\bf Y)\sim\mathcal{N}(vec(M),V).
$$
