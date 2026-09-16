library(AlphaSimR)
library(dplyr)
library(purrr)
library(here)

source(here("r_scripts", "03_metrics.R"))

schemes <- c("No marker", "Foreground", "Background", "Dense genomic")


# 1. Five-replicate profiling
# A short run at 100 plants and 3 generations to check run time before the full replication

profile_start <- Sys.time()

profile_runs <- 101:105 %>%
  map_dfr(\(s) map_dfr(schemes, \(x) run_scheme(x, fnd_single, 100, 10, 3, s)))

Sys.time() - profile_start

profile_runs %>%
  group_by(scheme, generation) %>%
  summarise(recovery_ibd = mean(recovery_ibd), linkage_drag_cM = mean(linkage_drag_cM), .groups = "drop")


# 2. Run replicated simulations
# 50 seeded replicates per scheme, 100 plants per backcross generation, 10 selected, 6 generations

seeds <- 1001:1050

replicate_runs <- seeds %>%
  map_dfr(\(s) map_dfr(schemes, \(x) run_scheme(x, fnd_single, 100, 10, 6, s))) %>%
  mutate(replicate = match(seed, seeds)) %>%
  relocate(replicate, .after = seed)

saveRDS(replicate_runs, here("data", "processed", "replicate_runs.rds"))

replicate_runs %>%
  count(scheme, generation)


# 3. Generations to recovery thresholds
# First generation at which mean IBD recovery of the selected plants reaches 95% and 99%

recovery_thresholds <- replicate_runs %>%
  group_by(replicate, seed, scheme) %>%
  summarise(
    gen_95 = first_gen(generation, recovery_ibd, 0.95),
    gen_99 = first_gen(generation, recovery_ibd, 0.99),
    .groups = "drop"
  )

recovery_thresholds %>%
  group_by(scheme) %>%
  summarise(
    reached_95 = sum(!is.na(gen_95)),
    mean_gen_95 = mean(gen_95, na.rm = TRUE),
    reached_99 = sum(!is.na(gen_99)),
    mean_gen_99 = mean(gen_99, na.rm = TRUE),
    .groups = "drop"
  )


# 4. Convergence of generations-to-99%
# Running mean over replicates for the two schemes that reach 99% in every replicate

convergence_99 <- recovery_thresholds %>%
  filter(scheme %in% c("Background", "Dense genomic")) %>%
  arrange(scheme, replicate) %>%
  group_by(scheme) %>%
  mutate(running_mean_99 = cummean(gen_99)) %>%
  ungroup()

convergence_99 %>%
  group_by(scheme) %>%
  slice_tail(n = 10) %>%
  select(scheme, replicate, running_mean_99) %>%
  print(n = Inf)


# 5. Save recovery and convergence results

saveRDS(recovery_thresholds, here("data", "processed", "recovery_thresholds.rds"))

saveRDS(convergence_99, here("data", "processed", "convergence_99.rds"))
