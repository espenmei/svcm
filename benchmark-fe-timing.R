# Timing: profiled fe() vs explicit pm() + mc() on a model with many fixed
# effects. Run with: Rscript benchmark-fe-timing.R
suppressMessages(devtools::load_all("."))

set.seed(20240917)

# Dimensions: I responses, each with k covariates => I * k fixed effects.
J <- 400   # individuals (rows of Y)
I <- 5     # responses (columns of Y)
k <- 25    # covariates per response (including intercept)
n_fe <- I * k

# Design matrix shared across responses: intercept + (k - 1) covariates.
Z <- matrix(rnorm(J * (k - 1)), J, k - 1)
X <- cbind(1, Z)

# True coefficients: I x k, one row per response.
B_true <- matrix(rnorm(I * k), I, k)

# Variance: shared random intercept (rank-1 across responses) + residual noise.
eta <- rnorm(J, 0, sqrt(1.0))
eps <- matrix(rnorm(J * I, 0, sqrt(0.5)), J, I)
Y   <- X %*% t(B_true) + matrix(rep(eta, I), J, I) + eps

R <- Matrix::Diagonal(J)

var_part <- function() list(
  pm(I, 1, rep("lp",  I), TRUE,          1,          "LP"),
  pm(I, I, rep("lth", I), diag(TRUE, I), diag(1, I), "LTH"),
  ic(tcrossprod(LP),  name = "P"),
  ic(tcrossprod(LTH), name = "TH"),
  svc(P + TH, R = R)
)

# Explicit model: free B (I x k) optimized jointly => n_fe extra parameters.
mod_explicit <- do.call(svcm, c(
  list(Y), var_part(),
  list(pm(I, k, sapply(1:k, \(j) paste0("b", 1:I, "_", j)), TRUE, 0, "B"),
       mc(B, X = X))
))

# Profiled model: the same fixed effects handled by fe(X).
mod_profiled <- do.call(svcm, c(
  list(Y), var_part(),
  list(fe(X))
))

cat(sprintf("Model: J=%d individuals, I=%d responses, k=%d covariates\n", J, I, k))
cat(sprintf("Fixed effects: %d   (optimized parameters: explicit=%d, profiled=%d)\n\n",
            n_fe, length(theta(mod_explicit)), length(theta(mod_profiled))))

t_explicit <- system.time(
  fit_explicit <- fit_svcm(mod_explicit,
                           control = list(iter.max = 500, eval.max = 700))
)
t_profiled <- system.time(fit_profiled <- fit_svcm(mod_profiled))

ll_e <- as.numeric(logLik(fit_explicit))
ll_p <- as.numeric(logLik(fit_profiled))

cat(sprintf("Explicit  pm()+mc():  %6.2f s   iters=%4d   conv=%d   logLik=%.4f\n",
            t_explicit[["elapsed"]], fit_explicit$opt$iterations,
            fit_explicit$opt$convergence, ll_e))
cat(sprintf("Profiled  fe():       %6.2f s   iters=%4d   conv=%d   logLik=%.4f\n",
            t_profiled[["elapsed"]], fit_profiled$opt$iterations,
            fit_profiled$opt$convergence, ll_p))
cat(sprintf("\nSpeed-up: %.1fx   (logLik difference: %.2e)\n",
            t_explicit[["elapsed"]] / t_profiled[["elapsed"]],
            abs(ll_e - ll_p)))
