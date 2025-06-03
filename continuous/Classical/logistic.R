######################
# title: model for logistic regression analysis
# date started: 11/02/25
# date finished:
# author: Ellie Van Vogt
#####################
set.seed(1998)

#libraries
library(stats)
library(dplyr)
library(furrr)

#paths
path <- "/rds/general/user/evanvogt/projects/nihr_drf_simulations"
setwd(path)

# arguments
args <- commandArgs(trailingOnly = TRUE)
scenario <- as.character(args[1])
n <- as.numeric(args[2])
n_cores <- 5

oldplan <- plan(multicore, workers = n_cores)

# load in the data
datasets <- readRDS(paste0(c("live/data/continuous/", scenario, "_", n, ".RDS"), collapse = ""))
datasets <- lapply(datasets, `[[`, 1) # just want the data not the truth


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
results <- future_map(datasets, logistic_analysis, .progress = T, .options = furrr_options(seed = T))
t1 <- Sys.time()
print(t1-t0)

plan(oldplan)

#save the results
saveRDS(results, paste0(c("live/results/continuous/", scenario, "/", n, "/Logistic/", "interactions.RDS"), collapse = ""))
