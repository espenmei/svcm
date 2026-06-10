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
