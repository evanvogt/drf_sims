##########
# title: extract nuisance / pseudo-value diagnostics - competing risk
##########
# Reduces surv_all.RDS to two small objects for plotting the saved nuisances
# (W.hat, pseudo.hat*, po) of every pseudo-value DR-learner arm, which nothing
# else in the pipeline reads:
#   nuisance_run_summary.RDS   one row per (rho, scenario, censoring, run, arm,
#                              estimand, quantity), all runs
#   nuisance_indiv_sample.RDS  one row per person, for N_RUNS_SAMPLE runs per
#                              combo
#
# HPC only: surv_all.RDS is ~729MB serialised. Parallelised over the 28
# (rho, scenario, censoring) combos (see the note above future_map()).

library(here)
library(dplyr)
library(tidyr)
library(purrr)
source(here("competing_risk/surv_config.R"))

# One future::multisession worker per task below; keep the jobscript's PBS
# ncpus/ompthreads in step with this (jobscripts/surv_nuisance_extract.sh).
# There are 28 (rho, scenario, censoring) combos, so this can go as high as 28
# given enough cores/memory; workers <= 1 runs sequentially instead.
workers <- 7

N_RUNS_SAMPLE <- 20 # runs sampled per (rho, scenario, censoring) combo for the individual-level object

# The nuisance arms nested under sim_res$nuisances$<base_learner>$<arm>, and
# the estimands each arm carries - see all_cate_surv_models() in surv_models.R.
# rf has arms whole_oob/whole_scf/cvps_scf, sl has whole/cvps, rsf (the
# randomForestSRC DR-learner) has oob/scf; base learner is kept as its own
# column so arm names don't collide across them, and arm_label =
# paste(base_learner, arm) is the single axis that distinguishes all 7.
BASE_LEARNERS <- c("rf", "sl", "rsf")
ESTIMANDS <- c("RMTL1", "RMTL2", "RMSTc")
QUANTITIES <- c("po", "pseudo.hat", "pseudo0.hat", "pseudo.hat.cf", "W.hat")

#' One run's nuisance leaves, flattened to a long tibble
#'
#' One row per (individual, base_learner, arm, estimand), wide on the 5
#' nuisance quantities (po, pseudo.hat, pseudo0.hat, pseudo.hat.cf, W.hat) so
#' they stay aligned per individual - needed for e.g. a po-vs-W.hat scatter.
extract_nuisance_long <- function(sim_res) {
  W <- sim_res$data$W
  n <- length(W)

  map_dfr(BASE_LEARNERS, function(bl) {
    arms <- sim_res$nuisances[[bl]]
    map_dfr(names(arms), function(arm) {
      estimand_list <- arms[[arm]]
      map_dfr(intersect(names(estimand_list), ESTIMANDS), function(est) {
        cell <- estimand_list[[est]]
        tibble(
          base_learner = bl,
          arm = arm,
          arm_label = paste(bl, arm, sep = "_"),
          estimand = est,
          id = seq_len(n),
          W = W,
          po = cell$po,
          pseudo.hat = cell$pseudo.hat,
          pseudo0.hat = cell$pseudo0.hat,
          pseudo.hat.cf = cell$pseudo.hat.cf,
          W.hat = cell$W.hat
        )
      })
    })
  })
}

#' Collapse one run's long nuisance tibble to summary stats per quantity
#'
#' prop_extreme flags |value - median| > 3 x IQR within the same (arm,
#' estimand, quantity) cell - only meaningful for po, which the AIPW
#' correction term can blow up when W.hat sits near the trim bounds.
#' prop_at_lo/prop_at_hi flag exact equality to the 0.05/0.95 trim_ps()
#' bounds (R/utils.R) - only meaningful for W.hat.
summarise_nuisance_run <- function(long) {
  long %>%
    pivot_longer(all_of(QUANTITIES), names_to = "quantity", values_to = "value") %>%
    group_by(base_learner, arm, arm_label, estimand, quantity) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      iqr = IQR(value, na.rm = TRUE),
      min = min(value, na.rm = TRUE),
      max = max(value, na.rm = TRUE),
      prop_extreme = mean(abs(value - median) > 3 * iqr, na.rm = TRUE),
      prop_at_lo = mean(abs(value - 0.05) < 1e-8, na.rm = TRUE),
      prop_at_hi = mean(abs(value - 0.95) < 1e-8, na.rm = TRUE),
      .groups = "drop"
    )
}

#' One (rho, scenario, censoring) combo's runs -> run-summary + individual sample
#'
#' Takes just this combo's own results list (list(run=, result=) entries),
#' not the full collected object - see the note above future_map() below for
#' why that's the piece worth parallelising over.
#'
#' Each combo samples its own N_RUNS_SAMPLE run ids independently - under
#' future_map(..., furrr_options(seed = TRUE)) every task gets its own
#' statistically sound RNG stream, not a shared, script-reproducible one, so
#' which runs land in nuisance_indiv_sample.RDS will differ run to run. That's
#' fine for an exploratory diagnostic; set workers <- 1 first if a fixed
#' set.seed() reproducible sample is ever needed.
process_combo <- function(results_list, rho, scenario, censoring) {
  run_ids <- map_int(results_list, "run")
  sampled_ids <- sample(run_ids, size = min(N_RUNS_SAMPLE, length(run_ids)))

  per_run <- map(results_list, function(entry) {
    long <- extract_nuisance_long(entry$result)

    summary <- summarise_nuisance_run(long) %>%
      mutate(rho = rho, scenario = scenario, censoring = censoring,
             run = entry$run, .before = 1)

    indiv <- if (entry$run %in% sampled_ids) {
      mutate(long, rho = rho, scenario = scenario, censoring = censoring,
             run = entry$run, .before = 1)
    } else {
      NULL
    }

    list(summary = summary, indiv = indiv)
  })

  message(
    "  finished rho=", rho, " scenario=", scenario, " censoring=", censoring,
    " (", length(run_ids), " runs)"
  )

  list(
    summary = bind_rows(map(per_run, "summary")),
    indiv = bind_rows(map(per_run, "indiv"))
  )
}

run_combo <- function(args) {
  process_combo(args$results_list, args$rho, args$scenario, args$censoring)
}

message("Reading surv_all.RDS...")
all_results_df <- readRDS(file.path(study$res_path, "surv_all.RDS"))

# Bundle each combo's own slice into its own list element BEFORE handing it to
# future_map(): future_map(.x, .f) ships .x[[i]] to worker i, so mapping over
# this pre-sliced list means each worker only receives its own ~1/28th of the
# collected data. Writing future_map(seq_len(nrow(all_results_df)), function(i)
# all_results_df$results[[i]] ...) instead would make every worker export the
# ENTIRE all_results_df as a captured global (its results list-column holds
# every run's data) - the opposite of what parallelising here is for.
combo_args <- pmap(
  list(all_results_df$results, all_results_df$rho, all_results_df$scenario,
       all_results_df$censoring),
  function(results_list, rho, scenario, censoring) {
    list(results_list = results_list, rho = rho, scenario = scenario,
         censoring = censoring)
  }
)

message(
  "Extracting nuisances from ", nrow(all_results_df), " (rho, scenario, censoring) combos, ",
  sum(lengths(all_results_df$results)), " runs total, ", workers, " worker(s)..."
)

if (workers > 1) {
  require(future)
  require(furrr)
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  processed <- furrr::future_map(
    combo_args,
    run_combo,
    .options = furrr::furrr_options(seed = TRUE)
  )
} else {
  processed <- map(combo_args, run_combo)
}

nuisance_run_summary <- bind_rows(map(processed, "summary"))
nuisance_indiv_sample <- bind_rows(map(processed, "indiv"))

saveRDS(nuisance_run_summary, file.path(study$res_path, "nuisance_run_summary.RDS"))
saveRDS(nuisance_indiv_sample, file.path(study$res_path, "nuisance_indiv_sample.RDS"))

message(
  "Done: ", nrow(nuisance_run_summary), " run-summary rows, ",
  nrow(nuisance_indiv_sample), " individual-level rows."
)
