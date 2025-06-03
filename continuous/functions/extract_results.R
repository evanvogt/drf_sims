extract_results <- function(result_type, scens, samplesizes, models) {
  result_list <- setNames(vector("list", length(scens)), scens)
  
  for (scenario in scens) {
    result_list[[scenario]] <- setNames(vector("list", length(samplesizes)), as.character(samplesizes))
    
    for (n in samplesizes) {
      result_list[[scenario]][[as.character(n)]] <- setNames(vector("list", length(models)), models)
      
      for (model in models) {
        dir_path <- file.path("live/results/continuous", scenario, n, model)
        
        if (dir.exists(dir_path)) {
          files <- list.files(dir_path, pattern = result_type, full.names = TRUE)
        }
        
        if (is.null(files)) {
          print(paste0(scenario, "_", n, " ", model, " ", result_type, " is not there - check !"))
        }
        
        if (length(files) > 1) {
          print(paste0(scenario, "_", n, " ", model, " ", result_type, " multiple files - check directory and remove non-relevant files"))
        }
        
        if (length(files) == 1) {
          res <- readRDS(files)
          result_list[[scenario]][[as.character(n)]][[model]] <- res
        }
        
      }
    }
  }
  
  return(result_list)
}