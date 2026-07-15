# svcm 0.6.0

## New features

* Added precomputed sparsity templates for the marginal covariance matrix
  `V = sum_i V_i` when all variance components are fixed `svc(..., R=)` terms.
  This avoids repeated sparse pattern unions in `expected_cov()` by refilling
  only numeric values each iteration.

* Reuse of the symbolic Cholesky factorization in `objective()`. Because the
  covariance sparsity pattern is fixed, the fill-reducing permutation and
  elimination structure are computed once in `svcm()` and each objective
  evaluation only redoes the numeric factorization via `Matrix::update()`.

## Testing

* Added `test-vtemplate.R` with correctness and regression checks for the
  `V` template and symbolic-Cholesky-reuse paths, including equality to the
  naive summation / full-factorization paths during covariance evaluation and
  model fitting, and guards that neither optimization slows objective().

## Documentation

* Documented the two complementary `svc()` forms: the separable
  `svc(Sigma, R = R)` form declaring a Kronecker term
  (fast, fixed sparsity pattern) and the general `svc(expr)` form for flexible
  covariance structures such as `Z (G %x% A) Z' + E %x% I`. Clarified in the
  `svc()` help, the README, and code comments why the marginal-covariance and
  symbolic-Cholesky precomputations apply only when every term is separable,
  and that mixed models fall back to correct-but-general summation.

# svcm 0.5.0

## New features

* Added `const()` constructor for supplying external constants (e.g. sparse
  matrices) to model expressions. Constants are explicitly declared alongside
  the other model components (`pm`, `svc`, `mc`, `ic`) and are available by
  name in all expressions evaluated in the model's computing environment.
  They are accessible post-fit via `mod$consts`.

## Bug fixes and improvements

* Model expressions involving `t()`, `%*%`, and other Matrix functions now
  resolve correctly without requiring `Matrix` to be explicitly attached.
  The model's computing environment now inherits from the svcm package
  namespace, making all imported Matrix functions available automatically.

* Fixed missing `importFrom` declarations for `AIC`, `BIC`, `logLik`,
  `nlminb`, `vcov` (stats) and `as` (methods).

* Added `R (>= 4.1.0)` dependency to reflect use of the native pipe `|>`
  and anonymous function shorthand `\()`.

## Documentation

* Added a `\donttest` example to `svcm()` demonstrating the genetic factor
  model from the README in two equivalent forms: the implicit `R=` Kronecker
  syntax and the explicit `const()` + `kronecker()` syntax.

## Testing

* Added a formal test suite (151 tests across five files) covering
  constructors, numerical parameter recovery, equivalence with lme4 and
  lavaan (complete and missing data), and delta-method SE equivalence
  between direct and Cholesky parameterisations.
