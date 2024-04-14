library(lavaan.bingof)
library(tidyverse)
# pak::pkg_install("davidsjoberg/ggsankey")
library(ggsankey)
theme_set(theme_sankey())

dat <- gen_data_bin(model_no = 1, n = 1000, seed = 123) %>%
  mutate(across(starts_with("y"), function(x) as.numeric(x == 1)))

# Create bivariate moments
id <- combn(5, 2)
biv_dat <- NULL
for (k in seq_len(ncol(id))) {
  i <- id[1, k]
  j <- id[2, k]
  biv_dat[[k]] <- dat[[i]] * dat[[j]]
}
biv_dat <- as.data.frame(biv_dat)

# Create vector of name arrangements (for fill purposes as well)
uni_bi_y <- c(paste0("y", 5:1), rev(paste0("y", id[1, ], id[2, ])))
zy <-
  tibble(y = sort(uni_bi_y),
         z = rev(paste0("z", sprintf("%02d", 1:15)))) %>%
  arrange(z)
names(biv_dat) <-
  zy %>%
  slice(match(paste0("y", id[1, ], id[2, ]), y)) %>%
  pull(z)

get_uni2 <- function(x) {
  if (is.na(x)) return(NA)
  y <- zy %>%
    slice(match(x, z)) %>%
    pull(y)
  x <- paste0("y", str_split(gsub("y", "", y), "")[[1]])
  zy %>%
    slice(match(x, y)) %>%
    pull(z)
}

plot_df <-
  dat %>%
  mutate(pattern = apply(across(everything()), 1, paste0, collapse = ""),
         data = "All\n(n x p)") %>%
  bind_cols(biv_dat) %>%
  
  # Flow from response patterns to bivariate marginals (positive responses)
  pivot_longer(starts_with("z"), names_to = "bivariate",
               values_to = "bi_vals") %>%
  mutate(bivariate = case_when(bi_vals == 0 ~ NA, TRUE ~ bivariate)) %>%
  
  # Flow from bivariate to univariate positive marginals
  rowwise() %>%
  mutate(univariate = list(get_uni2(bivariate))) %>%
  unnest_longer(univariate) %>%
  
  # Create the sankey data frame
  make_long(data, pattern, bivariate, univariate) %>%
  mutate(#col = my_col_fn(node),
    x = factor(x, labels = c("Multivariate\nBernoulli Data",
                             "Response\nPatterns", "Bivariate\nMoments",
                             "Univariate\nMoments")))

p_tmp <-
  ggplot(plot_df, aes(x = x, next_x = next_x, node = node,
                      next_node = next_node, fill = node)) +
  geom_sankey(flow.alpha = 0.5, space = 300, width = 0.25) +
  geom_sankey_text(aes(label = node), space = 300) +
  scale_fill_viridis_d(option = "turbo") +
  theme(legend.position = "none") +
  labs(x = NULL)

my_labeller <- function(x) {
  if (grepl("z", x)) {
    zy %>%
      slice(match(x, z)) %>%
      pull(y)
  } else {
    x
  }
}

# Extract data from plot builder
pb <- ggplot_build(p_tmp)
lab_dat <- pb$data[[3]] %>%
  as_tibble %>%
  rowwise() %>%
  mutate(labz = case_when(freq < 800 ~ NA, TRUE ~ my_labeller(node)))

p <-
  ggplot(plot_df, aes(x = x, next_x = next_x, node = node,
                      next_node = next_node,
                      fill = node)) +
  geom_sankey(flow.alpha = 0.5, space = 300, width = 0.3) +
  geom_label(data = lab_dat, aes(x, y, label = labz), fill = "gray80",
             alpha = 0.5, col = "gray20", inherit.aes = FALSE, size = 2.5) +
  scale_fill_viridis_d(option = "turbo") +
  scale_colour_viridis_d(option = "turbo") +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL)

p_mult_bern_fig <- p
save(p_mult_bern_fig, file = "R/p_mult_bern_fig.RData")
