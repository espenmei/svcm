README
================

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
environmental component

$$
\eta_j = a_j + e_j.
$$

The genetic components $a_j$ represent a additive polygenic contribution
to the common factor, and has covariance matrix $\sigma^2_a\bf{A}$.
$\bf{A}$ is a sparse symmetric matrix describing genetic relatedness
among individuals. When $\bf{A}$ cannot be arranged into small
homogeneous block, this model cannot be expressed in common software for
structural equation models. The environmental component is independent
across individuals with covariance matrix $\sigma^2_e\bf{I}$. The model
could also be written for the whole sample as

$$
\bf{Y} = \boldsymbol{x\beta}'+
\boldsymbol{\eta}\boldsymbol{\lambda}'+ \bf{\Delta}.
$$

Across variables and individuals, the covariance matrix implied by this
model can be described as $$
\bf{V}=
(\boldsymbol{\lambda\lambda'}\sigma^2_a)\otimes\bf{A}+
(\boldsymbol{\lambda\lambda'}\sigma^2_e+\bf{\Theta})\otimes\bf{I}.
$$ The mean vector implied by the model can be described as $$
\boldsymbol{m} = vec(\boldsymbol{x\beta'}).
$$

Finally, the data is described with the normal distribution $$
vec(Y) \sim \mathcal{N}(\boldsymbol{m},\bf{V}).
$$

## Simulate data

``` r
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

## Fit the model

Below is an example of one way to fit this model using **svcm**. `pm()`
is used to define parameters as matrices. `svc()` is used to define
structured variance components. Notice that the expressions use
parameter matrices to define the model implied covariance structure
across variables. The argument `R =` is used to supply a relationship
matrix describing the dependence across individuals. Notice also that
`mc()` is used to define the mean structure. There is no limit to how
many `pm`, `svc` and `mc` objects that can be supplied, they will be
added together to form the total mean and covariance structure.

``` r
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

After the model is defined, it can be fitted using `fit_svcm()`.

``` r
m1 <- fit_svcm(m1, se = TRUE, control = list(trace = 10))
```


    iter: objective:    l2  l3  l4  Sa1 Se1 th11    th22    th33    th44    b1  b2  b3  b4  
      0:     22106.472: 0.500000 0.500000 0.500000  1.00000  1.00000  1.00000  1.00000  1.00000  1.00000  0.00000  0.00000  0.00000  0.00000
     10:     10215.346: 0.836416 0.783884  1.32430  2.30893  1.30135  1.73818  1.11008  1.13220 0.883441 0.236877 0.586525  2.40862  1.73297
     20:     9470.8411: 0.515797 0.493898 0.812449  2.38003  1.56757  1.00878 0.978941 0.994438 0.917176  1.68961  1.78277  3.83229  3.71360
     30:     9445.9735: 0.516696 0.471899 0.811573  2.22885  1.71383 0.977803 0.997554  1.00774 0.933840  2.00202  1.97173  4.02065  3.99984
     40:     9444.0691: 0.515197 0.475552 0.814056  2.19358  1.79692 0.963652 0.979031 0.997702 0.936235  2.07450  2.01172  4.06147  4.07139
     50:     9443.9058: 0.513130 0.471452 0.812789  2.17720  1.82667 0.954312 0.979261 0.996704 0.941872  2.09843  2.02484  4.07148  4.09293
     60:     9443.8766: 0.513794 0.472560 0.812841  2.17247  1.84118 0.958335 0.976913 0.996158 0.940334  2.10768  2.02915  4.07649  4.09982
     70:     9443.8753: 0.513817 0.472672 0.812627  2.17101  1.84343 0.958403 0.976832 0.996103 0.940211  2.10994  2.03026  4.07737  4.10141
     80:     9443.8751: 0.513790 0.472664 0.812676  2.17016  1.84520 0.958447 0.976831 0.996075 0.940197  2.11099  2.03080  4.07788  4.10225

    Computing standard errors.

After the models is fitted, it can be inspected with `summary()`.

``` r
summary(m1)
```

    Log likelihood       Deviance            AIC            BIC 
         -4721.938       9443.875       9469.875       9547.061 

    Fitted parameters:
         Estimate Std. Error z value Pr(>|z|)
    l2     0.5138     0.0231   22.28   <0.001
    l3     0.4727     0.0229   20.63   <0.001
    l4     0.8127     0.0291   27.95   <0.001
    Sa1    2.1702     0.4952    4.38   <0.001
    Se1    1.8452     0.4427    4.17   <0.001
    th11   0.9584     0.1056    9.08   <0.001
    th22   0.9768     0.0596   16.39   <0.001
    th33   0.9961     0.0594   16.76   <0.001
    th44   0.9402     0.0804   11.70   <0.001
    b1     2.1110     0.0965   21.88   <0.001
    b2     2.0308     0.0591   34.37   <0.001
    b3     4.0779     0.0565   72.12   <0.001
    b4     4.1022     0.0812   50.55   <0.001

## In general

In general, **svcm** can fit models where the covariance matrix can be
described as $$
\bf{V}=\sum_{i=1}^{q}\mathbf{\Sigma}_i \otimes \mathbf{R}_i.
$$ The terms in the sum are specified by adding objects returned from
`svc()` to the model. $\mathbf \Sigma_i$ has dimension according to the
number of columns of $\bf{Y}$ and describes dependence among variables.
$\mathbf \Sigma_i$ is a function of parameters defined with `pm()` and
can be specified flexibly using general functions available in **R**.
$\mathbf{R}_i$ describe dependence across individuals and must be sparse
symmetric matrices with known values.

The means are described with a similar structure $$
\bf{M}=\sum_{i=1}^{p} \mathbf{X}_i \mathbf{B'}.
$$ For the means, the terms in the sum are specified by adding objects
returned from `mc()` to the model.
