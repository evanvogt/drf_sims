##########
# title: appendix figures - the true CATE structure of the DGMs
##########
# Appendix figures of the true CATE each DGM puts into a trial population at
# rho = 0.5 (dgm_population.R; no simulation output is read). Continuous /
# binary tau is drawn as tau(x) - ATE. The binary modifiers are the continuous
# ones with the sign reversed and scaled by RD_SCALE, so the binary row mirrors
# the continuous one.
#
# Writes to ../results/thesis_figures/dgm/:
#   dgm_cate_x4.png           tau - ATE against X4, coloured by X3
#   dgm_cate_dist.png         the distribution of tau - ATE
#   dgm_cate_sc3_surface.png  scenario 3's tau - ATE over (X4, X5) by X3
#   dgm_m0_cate.png           the control mean (risk) against tau - ATE
#   dgm_surv_cate.png         competing risks: the subdistribution and net RMST
#                             CATE of each event by X3
#   dgm_surv_cate_points.png  the same, as jittered points
#
#   cd results_processing/thesis_figures; Rscript dgm.R   # a few minutes

library(here)
library(patchwork)
source(here("R", "figures.R"))
source(here("results_processing", "dgm_population.R"))

# paths
path <- here()
fig_path <- file.path(dirname(path), "results", "thesis_figures", "dgm")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

OUTCOME_LABELS <- c(continuous = "Continuous", binary = "Binary")
# draws per panel in the scatter figures; the distributions use them all
N_PLOT <- 4000
X3_LABELS <- c(`0` = "X3 = 0", `1` = "X3 = 1")
Y_LAB <- "τ(x) − ATE\n(continuous: mean difference; binary: risk difference)"

# ---- truths -------------------------------------------------------------------

cv <- dgm_covariates()
pop <- cate_population(cv)
summ <- cate_summary(pop, cv)

label_cts_bin <- function(df) {
  df %>%
    mutate(
      outcome = label_factor(outcome, OUTCOME_LABELS),
      scenario = label_factor(scenario, SS_SCENARIO_LABELS)
    )
}

set.seed(DGM_SEED)
plot_ids <- sample(cv$id, N_PLOT)
pts <- pop %>%
  filter(id %in% plot_ids) %>%
  left_join(cv, by = "id") %>%
  mutate(X3 = label_factor(X3, X3_LABELS)) %>%
  label_cts_bin()
summ_l <- label_cts_bin(summ)

# ---- 1. tau - ATE against X4 --------------------------------------------------

ggplot(pts, aes(X4, tau_c, colour = X3)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 0.4, alpha = 0.4) +
  facet_grid(outcome ~ scenario, scales = "free_y") +
  drf_scale() +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "X4", y = Y_LAB, colour = NULL) +
  drf_theme() +
  theme(legend.position = "bottom")
save_fig("dgm_cate_x4.png", fig_path)

# ---- 2. distribution of tau - ATE ---------------------------------------------
# a violin of a constant (scenario 1) would be drawn with a bandwidth of 1, so
# the constant scenarios get the boxplot's line alone

dist <- label_cts_bin(pop)
varying <- dist %>%
  group_by(outcome, scenario) %>%
  filter(sd(tau_c) > 0) %>%
  ungroup()
sd_labels <- summ_l %>%
  left_join(
    dist %>% group_by(outcome) %>% summarise(top = max(tau_c), .groups = "drop"),
    by = "outcome"
  ) %>%
  mutate(label = sprintf("SD = %.*f", ifelse(outcome == "Binary", 3, 2), sd_tau))

ggplot(dist, aes(scenario, tau_c)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_violin(data = varying, scale = "width", fill = "grey85", colour = "grey40") +
  geom_boxplot(width = 0.12, outlier.shape = NA, coef = 0) +
  geom_text(data = sd_labels, aes(y = top, label = label), vjust = -0.4, size = 3) +
  facet_grid(outcome ~ ., scales = "free_y") +
  # the violin layer, trained first, has no Null, which would otherwise go last
  scale_x_discrete(limits = unname(SS_SCENARIO_LABELS)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = NULL, y = Y_LAB) +
  drf_theme()
save_fig("dgm_cate_dist.png", fig_path)

# ---- 3. scenario 3 over (X4, X5) ------------------------------------------------
# one fill scale per outcome, so one plot each, stacked

surface_grid <- expand.grid(X4 = seq(-3, 3, length.out = 121),
                            X5 = seq(-3, 3, length.out = 121), X3 = 0:1)
contour_pts <- cv %>%
  filter(id %in% sample(id, 20000)) %>%
  mutate(X3 = label_factor(X3, X3_LABELS))

surface_plot <- function(outcome) {
  cate_grid(outcome, 3, surface_grid) %>%
    mutate(X3 = label_factor(X3, X3_LABELS),
           outcome = label_factor(outcome, OUTCOME_LABELS)) %>%
    ggplot(aes(X4, X5)) +
    geom_raster(aes(fill = tau_c)) +
    geom_density_2d(data = contour_pts, colour = "grey20", linewidth = 0.3,
                    bins = 6) +
    facet_grid(outcome ~ X3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0) +
    coord_equal(expand = FALSE, xlim = c(-3, 3), ylim = c(-3, 3)) +
    labs(fill = "τ(x) − ATE") +
    drf_theme() +
    theme(panel.spacing.x = unit(1, "lines"))
}

(surface_plot("continuous") / surface_plot("binary"))
save_fig("dgm_cate_sc3_surface.png", fig_path, width = 18, height = 17)

# ---- 4. control mean against tau - ATE ------------------------------------------
# one plot per outcome: m0 is on each outcome's own scale

m0_plot <- function(outcome_label, x_lab) {
  d <- filter(pts, outcome == outcome_label)
  lab <- filter(summ_l, outcome == outcome_label) %>%
    mutate(label = ifelse(is.na(cor_m0_tau), "τ constant",
                          sprintf("r = %.2f", cor_m0_tau)))
  ggplot(d, aes(m0, tau_c)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 0.4, alpha = 0.3, colour = "grey30") +
    geom_text(data = lab, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.1, vjust = 1.4, size = 3) +
    facet_grid(outcome ~ scenario, scales = "free_y") +
    labs(x = x_lab, y = "τ(x) − ATE") +
    drf_theme()
}

(m0_plot("Continuous", "Control mean m0(x)") /
    m0_plot("Binary", "Control risk m0(x)"))
save_fig("dgm_m0_cate.png", fig_path)

# ---- 5. competing risks: subdistribution and net RMST CATE by X3 ------------------
# Both on the RMST scale, so they read the same way round - positive: treatment
# delays the event. The subdistribution RMST CATE is the RMTL CATE negated

SURV_ESTIMAND_PLOT_LABELS <- c(sub_RMST = "Subdistribution RMST (−RMTL)",
                               net_RMST = "Net RMST (cause-specific)")
EVENT_LABELS <- c(`Event 1` = "Event of interest", `Event 2` = "Competing event")

spop <- surv_population(cv, estimands = names(SURV_ESTIMAND_PLOT_LABELS)) %>%
  mutate(
    estimand = label_factor(estimand, SURV_ESTIMAND_PLOT_LABELS),
    event = label_factor(event, EVENT_LABELS),
    scenario = label_factor(scenario, SURV_SCENARIO_PLOT_LABELS),
    X3 = label_factor(X3, X3_LABELS)
  )
surv_varying <- spop %>%
  group_by(scenario, estimand, event, X3) %>%
  filter(sd(tau) > 1e-8) %>%
  ungroup()

surv_dodge <- position_dodge(width = 0.8)
# what the violin and points versions share: each estimand's ATE as a dotted
# line under the data, the X3 group means over every draw on top, and the
# facets, scales and theme
surv_ate <- geom_hline(
  data = summarise(spop, ate = mean(tau), .by = c(scenario, estimand, event)),
  aes(yintercept = ate, colour = estimand),
  linetype = "dotted", linewidth = 0.5, show.legend = FALSE
)
surv_mean <- stat_summary(data = spop, aes(fill = estimand), fun = mean,
                          geom = "point", position = surv_dodge, size = 1.8,
                          shape = 21, colour = "black", stroke = 0.3)
surv_common <- list(
  facet_grid(event ~ scenario, scales = "free_y",
             labeller = labeller(scenario = function(x) sub(", ", ",\n", x))),
  drf_scale(),
  scale_fill_paletteer_d("rcartocolor::Safe"),
  labs(x = NULL, y = "RMST CATE, 28-day horizon (days)", colour = NULL,
       fill = NULL),
  drf_theme(),
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())
)

ggplot(spop, aes(X3, tau, fill = estimand)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  surv_ate +
  geom_violin(data = surv_varying, position = surv_dodge, scale = "width",
              width = 0.75, alpha = 0.7, colour = "grey30", linewidth = 0.25) +
  surv_mean +
  surv_common
save_fig("dgm_surv_cate.png", fig_path, width = 33, height = 15)

# ---- 5b. the same, as jittered points ---------------------------------------------
# N_PLOT draws, not all DGM_N_POP_SURV: the tails show as single points, and a
# constant CATE (no effect on that event's own hazard) as a band at its value,
# where the violin version has only the mean

set.seed(DGM_SEED)
surv_pts <- filter(spop, id %in% sample(unique(spop$id), N_PLOT))

ggplot(surv_pts, aes(X3, tau, colour = estimand)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  surv_ate +
  geom_point(position = position_jitterdodge(jitter.width = 0.6,
                                             dodge.width = 0.8,
                                             seed = DGM_SEED),
             size = 0.3, alpha = 0.08) +
  surv_mean +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)),
         fill = "none") +
  surv_common
save_fig("dgm_surv_cate_points.png", fig_path, width = 33, height = 15)
