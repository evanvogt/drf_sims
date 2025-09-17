# inside your session after generating one dataset (or before simulation)
set.seed(1)
n <- 2000
dat_res <- generate_csh_data(1, n, return_truth = TRUE, seed = 1)
dat <- dat_res$dataset

# Inspect mean scale used in data generation (if you return them)
# If you don't return per-subject scale in dataset, compute same as in function:
W <- dat$W
X1 <- dat$X1; X2 <- dat$X2; X3 <- dat$X3
params <- survival_scenario_params[survival_scenario_params$scenario==1, ]

log_shape1 <- log(params$shape1_base) + params$b1_shape1 * X1 + params$b2_shape1 * X2
# shape1 <- exp(log_shape1)   # not needed for this check

# scale formula from your code:
log_scale1 <- log(params$scale1_base) + W * (params$ATE1 * params$bW_scale1 + params$HTE1 * params$b3_scale1 * X3)
scale1 <- exp(log_scale1)

mean(scale1[W==0])
mean(scale1[W==1])
mean(scale1[W==1]) / mean(scale1[W==0])   # should be ≈ exp(bW) if no HTE


median_weibull <- function(scale, shape) scale * (log(2)^(1/shape))

# baseline analytic median
median_base <- median_weibull(params$scale1_base, params$shape1_base)
# analytic treated median if bW applied: (scale * exp(bW))
median_treated_analytic <- median_weibull(params$scale1_base * exp(params$bW_scale1[1]), params$shape1_base) 
median_base; median_treated_analytic; median_base - median_treated_analytic

# empirical medians from data (cause==1)
median(dat$Y[dat$D == 1 & dat$W==0], na.rm=TRUE)
median(dat$Y[dat$D == 1 & dat$W==1], na.rm=TRUE)


# If generate_csh_data returns truth:
truth <- dat_res$truth
mean(truth$rmst_treat_cs1) - mean(truth$rmst_control_cs1)

# Empirical RMST (simple numerical integration by arm on cause-specific CIF)
# You can also estimate RMST by integrating Kaplan-Meier / CIF estimates; but quick check:
aggregate(Y ~ W, data = subset(dat, D==1), FUN = median)  # quick check
