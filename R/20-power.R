library(lavaan.bingof)
library(lavaan)
library(tidyverse)
library(furrr)
plan(multisession, workers = parallel::detectCores() - 2)

B <- 100
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
    samp_size = c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000),
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


res |>
  unnest(res) |>
  unnest(res) |>
  summarise(
    rej_rate = mean(pval < 0.05),
    crit = sd(pval < 0.05) / sqrt(B),
    .by = c(samp_size, model_no, name, sampling)
  ) |>
  mutate(sampling = factor(sampling, levels = c("SRS", "Cluster", "Strat-clust"))) |>
  ggplot(aes(samp_size, rej_rate)) +
  # geom_line(aes(col = sampling)) +
  geom_smooth(aes(col = sampling), se = FALSE) +
  # geom_ribbon(aes(ymin = rej_rate - crit, ymax = rej_rate + crit, fill = sampling), alpha = 0.2) +
  facet_grid(model_no ~ name) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  labs(x = "Sample size (n)", y = "Power", col = "Sampling method", fill = "Sampling method")
  
