######################
# title: model for logistic regression analysis
# date started: 11/02/25
# date finished:
# author: Ellie Van Vogt
#####################


#libraries
library(stats)
library(dplyr)
library(doParallel)

#paths
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# arguments
args <- commandArgs(trailingOnly = TRUE)
scenario <- as.character(args[1])
n <- as.numeric(args[2])

# load in the data
datasets <- readRDS(paste0(c("live/data/", scenario, "_", n, ".rds"), collapse = ""))
datasets <- lapply(datasets, `[[`, 1) # just want the data not the truth


# cores and folds
# not crossfitting this because I don't think it makes sense?
# available cores seems to not be working atm :()
n_cores <- 10 #floor(parallelly::availableCores() *0.9)



logistic_analysis <- function(data) {
  X <- as.matrix(data[, -c(1:2)])
  Y <- data$Y
  W <- data$W
  
  covs <- colnames(X)
  
  get_interaction <- function(cov) {
    fmla <- as.formula(paste0(c("Y ~ W*", cov), collapse = ""))
    int <- glm(fmla, family = binomial(link = "logit"), data)
    ests <- coef(summary(int))[4,]
    res <- c(cov, ests)
    names(res) <- c("name", "est", "std.err", "z.value", "p.value")
    return(res)
  }
  
  interactions <- bind_rows(lapply(covs, get_interaction))
  
  # adjustments for multiple testing
  interactions <- interactions %>% 
    mutate(
      p.BH = p.adjust(p.value, method = "BH"),
      p.Bonf = p.adjust(p.value, method = "bonferroni"))
  
  return(interactions)
}

# do the logistic regressions in parallel:
t0 <- Sys.time()
results <- mclapply(datasets, logistic_analysis, mc.cores = n_cores)
t1 <- Sys.time()
print(t1-t0)

#save the results
saveRDS(results, paste0(c("live/results/", scenario, "/", n, "/Logistic/", "interactions.RDS"), collapse = ""))
