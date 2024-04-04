library(lavaan.bingof)
library(lavaan)
library(tidyverse)
library(furrr)
plan(multisession, workers = parallel::detectCores() - 2)

B <- 250
power_sim <- function(i = 1, samp_size = 1000, model_no = 1) {
  
  pop <- make_population(model_no, H1 = TRUE, seed = 31324, Sigma2_attr = TRUE)
  set.seed(NULL)
  
  dat1 <- gen_data_bin_srs(pop, n = samp_size)
  dat2 <- gen_data_bin_clust(pop, n = samp_size)
  dat3 <- gen_data_bin_strcl(pop, n = samp_size)
  
  fit1 <- sem(
    txt_mod(model_no),
    data = dat1,
    estimator = "PML",
    std.lv = TRUE
  )
  fit2 <- sem(
    txt_mod(model_no),
    data = dat2,
    estimator = "PML",
    std.lv = TRUE,
    sampling.weights = "wt"
  )
  fit3 <- sem(
    txt_mod(model_no),
    data = dat3,
    estimator = "PML",
    std.lv = TRUE,
    sampling.weights = "wt"
  )

  bind_rows(
    all_tests(fit1, sim = i) |> mutate(sampling = "SRS"),
    all_tests(fit2, sim = i) |> mutate(sampling = "Cluster"),
    all_tests(fit3, sim = i) |> mutate(sampling = "Strat-clust"),
  )
}

res <-
  expand_grid(
    samp_size = c(500, 1000, 2500, 3750, 5000),
    model_no = 1:5
  ) |>
  mutate(res = list(NA))

for (k in seq_len(nrow(res))) {
  the_model_no <- res$model_no[[k]]
  the_samp_size <- res$samp_size[[k]]
  
  cli::cli_alert_info("Running model = {the_model_no} / samp_size = {the_samp_size}")
  res$res[[k]] <- future_map(
    .x = 1:B,
    .f = \(x) power_sim(i = x, samp_size = the_samp_size, model_no = the_model_no),
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
  )
  cat("\n")
}
save(res, file = "res_power.RData")

plot_df <-
  res |>
  unnest(res) |>
  unnest(res) |>
  summarise(
    rej_rate = mean(pval < 0.05),
    crit = sd(pval < 0.05) / sqrt(B),
    .by = c(samp_size, model_no, name, sampling)
  ) |>
  mutate(
    sampling = factor(sampling, levels = c("SRS", "Cluster", "Strat-clust")),
    name = factor(name, levels = levels(lavaan.bingof::res_srs_type1$name)),
    model_no = factor(model_no, labels = levels(lavaan.bingof::res_srs_type1$sim))
  )

p_power_samp <-
  ggplot(plot_df, aes(samp_size, rej_rate)) +
  geom_line(aes(col = sampling), linewidth = 0.6) +
  facet_grid(model_no ~ name) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top", plot.margin = unit(c(-8,0,0,0), "pt"),
        legend.box.spacing = unit(0, "pt")) +
  labs(x = "Sample size (n)", y = "Power", col = "Sampling method", fill = "Sampling method")
  
p_power_tests <-
  ggplot(plot_df, aes(samp_size, rej_rate)) +
  geom_line(aes(col = name), position = position_jitter(w = 100, h = 0),
            linewidth = 0.5) +
  facet_grid(model_no ~ sampling) +
  # scale_colour_manual(
  #   values = c(rep("#30123BFF", 3), "#7A0403FF",  "#FABA39FF", "#1AE4B6FF")
  # ) +
  # scale_linetype_manual(
  #   values = c("dashed", "dotted", "dotdash", rep("solid", 3))
  # ) +
  # ggsci::scale_color_ucscgb() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top", plot.margin = unit(c(-8,0,0,0), "pt"),
        legend.box.spacing = unit(0, "pt")) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(x = "Sample size (n)", y = "Power", col = "Sampling method", linetype = "Sampling method")

save(p_power_samp, p_power_tests, file = "R/p_power.RData")
