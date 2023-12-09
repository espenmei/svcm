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

Consider a factor model for individual $j$ $$
\boldsymbol{y}_j = \boldsymbol{\beta}+\boldsymbol{\lambda}\eta_j+\boldsymbol{\delta}_j.
$$ $\boldsymbol{\beta}$ is a vector of intercepts,
$\boldsymbol{\lambda}$ is a vector of loadings on the common factor
$\eta_j$, and $\boldsymbol{\delta}_j$ a vector of specific factors with
diagonal covariance matrix $\mathbf{\Theta}$. The factor model describes
dependence across variables. The common factor is specified with a
genetic and environmental component $$
\eta_j = a_j + e_j.
$$ The genetic components $a_j$ represent an additive polygenic
contribution to the common factor and has covariance matrix
$\sigma^2_a\bf{A}$. $\bf{A}$ is a sparse symmetric matrix describing
genetic relatedness among individuals. When $\bf{A}$ cannot be arranged
into small homogeneous blocks, this model cannot be expressed in common
software for structural equation models. The environmental component is
independent across individuals with covariance matrix
$\sigma^2_e\bf{I}$. The model could also be written for the whole sample
as $$
\bf{Y} = \boldsymbol{x\beta}'+
\boldsymbol{\eta}\boldsymbol{\lambda}'+ \bf{\Delta}.
$$ Across variables and individuals, the covariance matrix implied by
this model can be described as $$
(\boldsymbol{\lambda\lambda'}\sigma^2_a)\otimes\bf{A}+
(\boldsymbol{\lambda\lambda'}\sigma^2_e+\bf{\Theta})\otimes\bf{I}.
$$ \### Simulate data Below is a simulation of data from the model.

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
      0:     21751.783: 0.500000 0.500000 0.500000  1.00000  1.00000  1.00000  1.00000  1.00000  1.00000  0.00000  0.00000  0.00000  0.00000
     10:     10560.806: 0.960753  1.04119  1.38934  2.15289  1.39648  1.79359  1.15944  1.11584 0.956805 -0.0782487 0.212850  1.69550  1.09147
     20:     9559.0309: 0.547862 0.524641 0.824809  2.52754  1.41086  1.03242  1.11766 0.971392  1.11562  1.76800  1.87849  3.86330  3.83325
     30:     9550.8107: 0.540981 0.528732 0.801443  2.53855  1.44212  1.01467  1.03186 0.955463  1.05967  1.90518  1.96784  3.94577  3.94770
     40:     9550.3330: 0.540997 0.530662 0.797523  2.65134  1.39314  1.01258  1.03468 0.951310  1.05766  1.93630  1.98750  3.96454  3.96422
     50:     9550.1566: 0.542701 0.530032 0.795922  2.75924  1.31475  1.01076  1.02960 0.954726  1.06008  1.94274  1.98852  3.96630  3.97486
     60:     9550.1127: 0.542638 0.529317 0.795763  2.83432  1.26367  1.00936  1.03089 0.954207  1.05947  1.94831  1.99187  3.96927  3.97841
     70:     9550.1078: 0.542634 0.529409 0.795278  2.86711  1.23993  1.00959  1.03072 0.954356  1.06017  1.94882  1.99219  3.96963  3.97882

    Computing standard errors.

After the models is fitted, it can be inspected with `summary()`.

``` r
summary(m1)
```

    Log likelihood       Deviance            AIC            BIC 
         -4775.054       9550.108       9576.108       9653.293 

    Fitted parameters:
         Estimate Std. Error z value Pr(>|z|)
    l2     0.5425     0.0239   22.68   <0.001
    l3     0.5293     0.0232   22.86   <0.001
    l4     0.7954     0.0290   27.42   <0.001
    Sa1    2.8723     0.5491    5.23   <0.001
    Se1    1.2353     0.4465    2.77   0.0057
    th11   1.0091     0.1061    9.51   <0.001
    th22   1.0309     0.0638   16.16   <0.001
    th33   0.9545     0.0594   16.07   <0.001
    th44   1.0598     0.0827   12.82   <0.001
    b1     1.9486     0.1011   19.27   <0.001
    b2     1.9921     0.0637   31.27   <0.001
    b3     3.9696     0.0618   64.19   <0.001
    b4     3.9787     0.0841   47.32   <0.001

## In general

**svcm** can fit models where the covariance matrix can be described as
$$
\bf{V}=\sum_{i=1}^{q}\mathbf \Sigma_i \otimes \mathbf R_i.
$$ The terms in the sum are specified by adding objects returned from
`svc()` to the model. $\mathbf \Sigma_i$ has dimension according to the
number of columns of $\bf{Y}$ and describes dependence among variables.
$\mathbf \Sigma_i$ is a function of parameters defined with `pm()` and
can be specified flexibly using general functions available in **R**.
$\mathbf R_i$ describe dependence across individuals and must be sparse
symmetric matrices with known values.

The means are described with a similar structure $$
\bf M =\sum_{i=1}^{p} \mathbf X_i \mathbf B_i^{'}.
$$

For the means, the terms in the sum are specified by adding objects
returned from `mc()` to the model. $\mathbf B_i^{'}$ has rows equal to
the number of columns of $\bf Y$, and columns equal to the columns of
$\mathbf X_i$.

The data is described with the normal distribution $$
vec(\bf Y)\sim\mathcal{N}(vec(M),V).
$$
