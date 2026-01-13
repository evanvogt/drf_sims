# Combine into one tidy data frame for easier plotting/analysis
ci_df <- bind_rows(
  lapply(names(ci_length), function(scn) {
    lapply(names(ci_length[[scn]]), function(sz) {
      lapply(names(ci_length[[scn]][[sz]]), function(mod_name) {
        data.frame(
          scenario = scn,
          size = sz,
          model = mod_name,
          length = as.vector(ci_length[[scn]][[sz]][[mod_name]])
        )
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
)


res_ci <- ci_df %>%
  filter(scenario %in% paste0("scenario_", c(1, 3, 8, 9))) %>%
  filter(model != "T_RF") %>%
  mutate(
    scenario = factor(
      recode(scenario,
             `scenario_1` = "Null",
             `scenario_3` = "Simple",
             `scenario_8` = "Complex",
             `scenario_9` = "Non-linear"),
      levels = c("Null", "Simple", "Complex", "Non-linear")
    ),
    n = factor(
      recode(size,
             `250` = "250",
             `500` = "500",
             `1000` = "1000"),
      levels = c("250", "500", "1000")
    ),
    model = factor(
      recode(model,
             `CF` = "Causal forest",
             `DR_RF` = "DR-RandomForest",
             `DR_oracle` = "DR-oracle"),
      levels = c("Causal forest", "DR-RandomForest", "DR-oracle")
    )
  ) %>%
  select(-size)

res_ci %>%
  ggplot(aes(x = n, y = length, color = model, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~scenario, nrow = 2) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "CI length",
       x = "Sample size")


#####################
coverage_df <- bind_rows(
  lapply(names(cate_coverage), function(scn) {
    lapply(names(cate_coverage[[scn]]), function(sz) {
      lapply(names(cate_coverage[[scn]][[sz]]), function(mod_name) {
        data.frame(
          scenario = scn,
          size = sz,
          model = mod_name,
          length = as.vector(cate_coverage[[scn]][[sz]][[mod_name]])
        )
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
)

res_coverage <- coverage_df %>%
  filter(scenario %in% paste0("scenario_", c(1, 3, 8, 9))) %>%
  filter(model != "T_RF") %>%
  mutate(
    scenario = factor(
      recode(scenario,
             `scenario_1` = "Null",
             `scenario_3` = "Simple",
             `scenario_8` = "Complex",
             `scenario_9` = "Non-linear"),
      levels = c("Null", "Simple", "Complex", "Non-linear")
    ),
    n = factor(
      recode(size,
             `250` = "250",
             `500` = "500",
             `1000` = "1000"),
      levels = c("250", "500", "1000")
    ),
    model = factor(
      recode(model,
             `CF` = "Causal forest",
             `DR_RF` = "DR-RandomForest",
             `DR_oracle` = "DR-oracle"),
      levels = c("Causal forest", "DR-RandomForest", "DR-oracle")
    )
  ) %>%
  select(-size)

res_coverage %>%
  ggplot(aes(x = n, y = length, color = model, fill = model)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~scenario, nrow = 2) +
  scale_fill_paletteer_d("rcartocolor::Safe") +
  scale_color_paletteer_d("rcartocolor::Safe") +
  theme_minimal() +
  labs(y = "Coverage",
       x = "Sample size")
