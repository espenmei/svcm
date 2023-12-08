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
\bf{Y} = \boldsymbol{x+\beta}'+
\boldsymbol{\eta}\boldsymbol{\lambda}'+ \bf{\Delta}.
$$

Across variables and individuals, the model implied covariance matrix
can be described as

$$
(\boldsymbol{\lambda\lambda'}\sigma^2_a)\otimes\bf{A}+
(\boldsymbol{\lambda\lambda'}\sigma^2_e+\bf{\Theta})\otimes\bf{I}.
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
many `svc()` and `mc()` objects that can be supplied, they will be added
together to form the total mean and covariance structure.

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
      0:     21904.122: 0.500000 0.500000 0.500000  1.00000  1.00000  1.00000  1.00000  1.00000  1.00000  0.00000  0.00000  0.00000  0.00000
     10:     10155.045: 0.827019 0.777851  1.34596  2.25117  1.20461  1.68512 0.993091  1.19810 0.983719 0.256173 0.613845  2.48003  1.70539
     20:     9431.7154: 0.537764 0.456334 0.785880  2.63279  1.26836 0.905099 0.968953  1.05421  1.06679  1.85597  1.88869  3.97220  3.88008
     30:     9419.0540: 0.519254 0.448742 0.844468  2.48543  1.16431  1.02036 0.909013  1.03865  1.00729  2.02755  1.97771  4.03884  4.01481
     40:     9418.5412: 0.523799 0.451872 0.839308  2.43718  1.16234  1.03105 0.902872  1.03938  1.01187  2.05927  1.99590  4.05825  4.04711
     50:     9418.5070: 0.525198 0.452832 0.840405  2.41241  1.17384  1.03536 0.903379  1.03813  1.00820  2.06877  2.00046  4.06201  4.05406
     60:     9418.5016: 0.525543 0.452490 0.840704  2.39881  1.18186  1.03573 0.903324  1.03771  1.00818  2.07186  2.00215  4.06344  4.05670
     70:     9418.5006: 0.525653 0.452579 0.840814  2.39005  1.18737  1.03639 0.903219  1.03769  1.00828  2.07358  2.00306  4.06422  4.05805
     80:     9418.5006: 0.525690 0.452591 0.840847  2.38910  1.18781  1.03658 0.903193  1.03773  1.00829  2.07374  2.00316  4.06429  4.05819

    Computing standard errors.

After the models is fitted, it can be inspected with `summary()`.

``` r
summary(m1)
```

    Log likelihood       Deviance            AIC            BIC 
         -4709.250       9418.501       9444.501       9521.686 

    Fitted parameters:
         Estimate Std. Error z value Pr(>|z|)
    l2     0.5257     0.0248   21.18   <0.001
    l3     0.4526     0.0246   18.39   <0.001
    l4     0.8409     0.0324   25.98   <0.001
    Sa1    2.3890     0.4703    5.08   <0.001
    Se1    1.1879     0.3914    3.03   0.0024
    th11   1.0366     0.1069    9.70   <0.001
    th22   0.9032     0.0569   15.88   <0.001
    th33   1.0377     0.0611   16.98   <0.001
    th44   1.0083     0.0847   11.91   <0.001
    b1     2.0738     0.0949   21.84   <0.001
    b2     2.0032     0.0581   34.50   <0.001
    b3     4.0643     0.0550   73.89   <0.001
    b4     4.0582     0.0823   49.34   <0.001
