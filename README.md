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
model can be described as

$$
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
      0:     21278.520: 0.500000 0.500000 0.500000  1.00000  1.00000  1.00000  1.00000  1.00000  1.00000  0.00000  0.00000  0.00000  0.00000
     10:     10162.318: 0.792366 0.812250  1.20572  2.33724  1.36836  2.00308 0.911362 0.980087 0.601323 -0.140794 0.348083  2.43695  1.49357
     20:     9464.0026: 0.536531 0.528660 0.815410  2.76205  1.57875  1.23264 0.967970 0.906137 0.848656  1.46795  1.64950  3.74448  3.51519
     30:     9438.4227: 0.509054 0.499686 0.802022  2.51326  1.58941  1.16397 0.953687 0.945138 0.866413  1.82941  1.85522  3.94124  3.81643
     40:     9436.5278: 0.505799 0.502774 0.796290  2.24901  1.80736  1.12947 0.950205 0.943716 0.886813  1.88930  1.88995  3.97082  3.86238
     50:     9436.0623: 0.504551 0.504001 0.797636  2.09106  1.93504  1.12729 0.952978 0.940317 0.886875  1.89121  1.89125  3.97327  3.86504
     60:     9435.7612: 0.506152 0.503231 0.797753  1.91021  2.07479  1.12302 0.953801 0.943535 0.889764  1.89201  1.89158  3.97282  3.86487
     70:     9435.6783: 0.505112 0.502485 0.798160  1.80524  2.18601  1.12564 0.954498 0.941676 0.889545  1.89308  1.89183  3.97412  3.86565
     80:     9435.6745: 0.505145 0.502728 0.798128  1.77579  2.20829  1.12502 0.954498 0.941549 0.890174  1.89324  1.89149  3.97358  3.86560

    Computing standard errors.

After the models is fitted, it can be inspected with `summary()`.

``` r
summary(m1)
```

    Log likelihood       Deviance            AIC            BIC 
         -4717.837       9435.674       9461.674       9538.860 

    Fitted parameters:
         Estimate Std. Error z value Pr(>|z|)
    l2     0.5051     0.0234   21.63   <0.001
    l3     0.5027     0.0233   21.54   <0.001
    l4     0.7981     0.0290   27.57   <0.001
    Sa1    1.7758     0.5175    3.43   <0.001
    Se1    2.2083     0.4832    4.57   <0.001
    th11   1.1250     0.1099   10.23   <0.001
    th22   0.9545     0.0586   16.30   <0.001
    th33   0.9415     0.0579   16.26   <0.001
    th44   0.8902     0.0761   11.70   <0.001
    b1     1.8932     0.0954   19.84   <0.001
    b2     1.8915     0.0572   33.04   <0.001
    b3     3.9736     0.0569   69.81   <0.001
    b4     3.8656     0.0778   49.71   <0.001

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

The data is described with the normal distribution

$$
vec(\bf Y)\sim\mathcal{N}(vec(M),V).
$$
