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
      0:     21670.779: 0.500000 0.500000 0.500000  1.00000  1.00000  1.00000  1.00000  1.00000  1.00000  0.00000  0.00000  0.00000  0.00000
     10:     10604.141: 0.867177  1.12006  1.44169  2.11895  1.40718  1.83659  1.16146  1.17543 0.879168 -0.0608543 0.285527  1.56786 0.980783
     20:     9592.7116: 0.447792 0.520538 0.849948  2.48137  1.37718  1.10695  1.05572  1.03906 0.881716  1.63752  1.79455  3.76144  3.70677
     30:     9573.9094: 0.475462 0.547149 0.842655  2.47775  1.46748  1.16306  1.05576  1.04676 0.861814  1.90598  1.90675  3.93347  3.88661
     40:     9572.7623: 0.479138 0.539592 0.841644  2.48820  1.51011  1.16558  1.04557  1.04026 0.858474  1.96786  1.93459  3.95975  3.93634
     50:     9572.6268: 0.476961 0.537382 0.839826  2.51350  1.51146  1.16084  1.04814  1.04115 0.856988  1.98873  1.94416  3.97218  3.95696
     60:     9572.6058: 0.477483 0.537880 0.840572  2.53468  1.50046  1.16511  1.04749  1.03951 0.856532  1.99697  1.94716  3.97684  3.96328
     70:     9572.6035: 0.477073 0.537653 0.840460  2.55120  1.49044  1.16473  1.04718  1.03974 0.856807  1.99926  1.94828  3.97804  3.96522
     80:     9572.6035: 0.477063 0.537632 0.840431  2.55201  1.48979  1.16466  1.04717  1.03977 0.856862  1.99956  1.94843  3.97820  3.96541

    Computing standard errors.

After the models is fitted, it can be inspected with `summary()`.

``` r
summary(m1)
```

    Log likelihood       Deviance            AIC            BIC 
         -4786.302       9572.604       9598.604       9675.789 

    Fitted parameters:
         Estimate Std. Error z value Pr(>|z|)
    l2     0.4771     0.0235   20.28   <0.001
    l3     0.5376     0.0243   22.11   <0.001
    l4     0.8404     0.0292   28.73   <0.001
    Sa1    2.5520     0.5182    4.92   <0.001
    Se1    1.4898     0.4406    3.38   <0.001
    th11   1.1647     0.1092   10.66   <0.001
    th22   1.0472     0.0624   16.77   <0.001
    th33   1.0398     0.0641   16.23   <0.001
    th44   0.8569     0.0784   10.93   <0.001
    b1     1.9996     0.1001   19.97   <0.001
    b2     1.9484     0.0583   33.42   <0.001
    b3     3.9782     0.0625   63.68   <0.001
    b4     3.9654     0.0844   46.96   <0.001
