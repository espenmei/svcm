# pm() errors without a name

    Code
      pm(2, 2, c("a", "b", "c", "d"), TRUE, 1)
    Condition
      Error in `pm()`:
      ! A name is required.

# ic() errors without a name

    Code
      ic(A %*% t(A))
    Condition
      Error in `ic()`:
      ! A name is required.

# const() errors without a name

    Code
      const(Matrix::Diagonal(5))
    Condition
      Error in `const()`:
      ! A name is required.

# svc() errors when R is not a Matrix

    Code
      svc(A, R = matrix(c(1, 0, 0, 1), 2, 2))
    Condition
      Error in `svc()`:
      ! inherits(R, "Matrix") is not TRUE

# svc() errors when R is not symmetric

    Code
      svc(A, R = R_asym)
    Condition
      Error in `svc()`:
      ! Matrix::isSymmetric(R) is not TRUE

# mc() warns when X is rank-deficient

    Code
      mc(B, X = X_rankdef)
    Condition
      Warning in `mc()`:
      The matrix X is not of full rank.
    Output
      Mean model:
      B
      
      Dimension of X:
      5 x 2 

# svcm() errors when pm is missing

    Code
      svcm(d$Y, svc(P, R = d$R), mc(U, X = d$X))
    Condition
      Error in `svcm()`:
      ! At least one pm, svc and mc object must be supplied.

# svcm() errors when svc is missing

    Code
      svcm(d$Y, d$P, d$U, mc(U, X = d$X))
    Condition
      Error in `svcm()`:
      ! At least one pm, svc and mc object must be supplied.

# svcm() errors when mc is missing

    Code
      svcm(d$Y, d$P, d$U, svc(P, R = d$R))
    Condition
      Error in `svcm()`:
      ! At least one pm, svc and mc object must be supplied.

