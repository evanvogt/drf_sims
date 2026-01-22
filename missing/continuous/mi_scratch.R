##########
# Title: multiple imputation procedure
##########

path <- here()
# generate cts scenario data and add the missingness
source(here("utils.R"))
source(here("missing/continuous/cts_miss_dgms.R"))
source(here("missing/continuous/cts_miss_models.R"))

data_result <- generate_continuous_scenario_data(
  scenario = 4,
  n = 500,
  mech = "AUX"
)

miss_dataset <- introduce_missingness_continuous(
  data = data_result$dataset,
  type = "both",
  prop = 0.3,
  mech = "AUX",
  U = U
)


predMat <- matrix(1, nrow = ncol(data_result$dataset), ncol = ncol(data_result$dataset))
diag(predMat) <- 0
predMat[,which(colnames(data_result$dataset) %in% c("Y", "W"))] <- 0
colnames(data_result$dataset)

imputation <- mice(miss_dataset, m = 10, predictorMatrix = predMat)

data_mi <- lapply(seq_len(10), function(i) {
  complete(imputation, i)
})


# ok now we added the MI imputation to the DGM file, need to figure out how the modelling changes

# can I just wrap the usual function in a lapply?
fmla_info <- get_continuous_oracle_info(scenario, data_result$bW)

short_mi <- list(data_mi[[1]], data_mi[[2]], data_mi[[3]])
result_list <- future_map(short_mi, function(data) {
  run_all_cate_methods(
    data = data,
    n_folds = n_folds,
    workers = 3,
    sl_lib = NULL,
    fmla_info = NULL
  )
}, .options = furrr_options(seed = TRUE))


# use the result list to do the pooling
combine_mi <- function(res_list, model) {
  res <- list()
  # combine point estimates
  tau_list <- lapply(res_list, function(x) x[[model]][["tau"]])
  tau_mat <- do.call(cbind, tau_list)
  tau <- rowMeans(tau_mat)
  res$tau <- tau
  
  # combine variance estimates if they exist
  var_list <- lapply(res_list, function(x) x[[model]][["variance"]])
  var_mat <- do.call(cbind, var_list)
  tot_var <- NULL
  if (!is.null(var_mat)) {
    w_var <- rowMeans(var_mat)
    b_var <- apply(tau_mat, 1, var)
    tot_var <- w_var + (1 + 1/length(res_list))*b_var
    res$variance <- tot_var
  }
  
  return(res)
}

results <- list()
results$causal_forest <- combine_mi(result_list, "causal_forest")
results$dr_random_forest <- combine_mi(result_list, "dr_random_forest")
results$dr_semi_oracle <- combine_mi(result_list, "dr_semi_oracle")


# now need to think about how the bootstrapping is gonna work 

