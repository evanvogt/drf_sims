require(parallel)
setup_rng_stream <- function(trial_num, seed = 1998) {
  RNGkind("L'Ecuyer-CMRG")
  
  set.seed(seed)

  for(i in 1:trial_num) {
    .Random.seed <<- parallel::nextRNGStream(.Random.seed)
  }
}
