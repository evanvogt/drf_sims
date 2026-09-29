###############
# title: continuous DGMs with missingness - CI example
###############
# Same scenarios and machinery as missing/continuous/, reached through
# R/dgm_scenarios.R and R/missingness.R.
#
# This study used to call the mechanisms "AUX" / "AUX-Y"; since 2026-09-28 the
# only names are MAR / MNAR-Y0 / MNAR-tau (R/dgm_scenarios.R). Its parameter
# grid only ever runs MAR, so the unobserved-U branches never fire here - which
# is also why its divergent p1 definition (it did not subtract the U term) made
# no difference to any result it produced. It does run on continuous_missing's
# correlated covariates (X01-X03 auxiliary), since 2026-09-28.
#
# The old header claimed the imputation count was reduced "from 50 to 20" to
# keep CI generation tractable, but the code set n_imp <- 50 and the completion
# message printed "<n_imp> 50 imputed datasets". 50 is what actually ran.

source(here::here("R", "missingness.R"))

generate_continuous_scenario_data <- function(scenario, n, mech, return_truth = TRUE) {
  generate_scenario_data(scenario, n, set = "continuous_missing",
                         return_truth = return_truth, mech = mech)
}

introduce_missingness_continuous <- function(data, type, prop, mech, U = NULL) {
  introduce_missingness(data, type, prop, mech, U = U)
}

handle_missingness_continuous <- function(data, method) {
  handle_missingness(data, method, n_imp = 50)
}

generate_and_process_continuous_data <- function(scenario, n, return_truth = TRUE,
                                                 type, prop, mech, method) {
  generate_and_process_data(scenario, n, set = "continuous_missing",
                            return_truth = return_truth, type = type,
                            prop = prop, mech = mech, method = method,
                            n_imp = 50)
}

get_continuous_oracle_info <- function(scenario, bW) {
  get_oracle_info(scenario, bW, set = "continuous_missing")
}
