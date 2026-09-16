library(AlphaSimR)
library(dplyr)
library(purrr)
library(tidyr)
library(here)

source(here("r_scripts", "03_metrics.R"))


# 1. Sweep parameters
# Background and Dense genomic at 50, 100, 200 and 400 plants per generation
# 10 plants are selected at every size, so the selected proportion falls from 20% to 2.5%

pop_sizes <- c(50, 100, 200, 400)
seeds <- 2001:2050
schemes <- c("Background", "Dense genomic")

sweep_grid <- expand_grid(seed = seeds, scheme = schemes, pop_size = pop_sizes) %>%
  mutate(replicate = match(seed, seeds))


# 2. Run full population-size sweep

sweep_runs <- pmap_dfr(
  sweep_grid,
  \(seed, scheme, pop_size, replicate) {
    run_scheme(scheme, fnd_single, pop_size, 10, 6, seed) %>%
      mutate(replicate = replicate, .after = seed)
  }
)

saveRDS(sweep_runs, here("data", "processed", "popsize_sweep.rds"))

sweep_runs %>%
  count(scheme, pop_size, generation)


# 3. Summarise population-size effects

sweep_summary <- sweep_runs %>%
  group_by(scheme, pop_size, generation) %>%
  summarise(
    recovery_ibd = mean(recovery_ibd),
    recovery_marker = mean(recovery_marker),
    linkage_drag_cM = mean(linkage_drag_cM),
    genotyping_cost = mean(genotyping_cost),
    .groups = "drop"
  )

sweep_thresholds <- sweep_runs %>%
  group_by(replicate, seed, scheme, pop_size) %>%
  summarise(gen_99 = first_gen(generation, recovery_ibd, 0.99), .groups = "drop")

sweep_thresholds %>%
  group_by(scheme, pop_size) %>%
  summarise(reached_99 = sum(!is.na(gen_99)), mean_gen_99 = mean(gen_99, na.rm = TRUE), .groups = "drop")


# 4. Cost to reach 99% recovery
# Genotyping cost is summed up to and including the first generation at 99% IBD recovery

cost_to_99 <- sweep_runs %>%
  group_by(replicate, seed, scheme, pop_size) %>%
  arrange(generation) %>%
  mutate(cumulative_cost = cumsum(genotyping_cost)) %>%
  filter(recovery_ibd >= 0.99) %>%
  slice_head(n = 1) %>%
  ungroup()

cost_to_99 %>%
  group_by(scheme, pop_size) %>%
  summarise(mean_gen_99 = mean(generation), mean_cost_99 = mean(cumulative_cost), .groups = "drop")


# 5. Save sweep summaries

saveRDS(sweep_summary, here("data", "processed", "popsize_summary.rds"))

saveRDS(sweep_thresholds, here("data", "processed", "popsize_thresholds.rds"))

saveRDS(cost_to_99, here("data", "processed", "cost_to_99.rds"))
