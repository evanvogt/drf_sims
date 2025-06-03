extract_results <- function(result_type, scens, samplesizes, models) {
  result_list <- setNames(vector("list", length(scens)), scens)
  
  for (scenario in scens) {
    result_list[[scenario]] <- setNames(vector("list", length(samplesizes)), as.character(samplesizes))
    
    for (n in samplesizes) {
      result_list[[scenario]][[as.character(n)]] <- setNames(vector("list", length(models)), models)
      
      for (model in models) {
        dir_path <- file.path("live/results", scenario, n, model)
        
        if (dir.exists(dir_path)) {
          files <- list.files(dir_path, pattern = "res_sim", full.names = TRUE)
          
          if (length(files) != 1000 & length(files) != 0) {
            message(paste0(scenario, " ", n, " ", model, " is incomplete, check which files need rerunning"))
          }
          
          if (length(files) == 1000) {
            results <- lapply(seq_len(1000), function(i) {
              res <- readRDS(file.path(dir_path, paste0("res_sim_", i, ".RDS")))
              res[[result_type]]
            })
            result_list[[scenario]][[as.character(n)]][[model]] <- results
          }
        }
      }
    }
  }
  
  return(result_list)
}