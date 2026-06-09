# fit_svcm() errors on non-svcm input

    Code
      fit_svcm(list())
    Condition
      Error in `fit_svcm()`:
      ! Only objects of type svcm are accepted.

# fit_svcm() warns when model is already fitted

    Code
      invisible(fit_svcm(fit))
    Condition
      Warning in `fit_svcm()`:
      This model has already been fitted.

# logLik() errors when model has not been fitted

    Code
      logLik(mod)
    Condition
      Error in `.assert_fitted_svcm()`:
      ! The model has not been fitted.

