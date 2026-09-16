library(AlphaSimR)
library(dplyr)
library(purrr)
library(here)

source(here("r_scripts", "03_metrics.R"))


# 1. Add second target locus
# TargetQTL2 sits at the centre of chromosome 5, unlinked to TargetQTL on chromosome 1
# The recurrent parent carries 0 and the donor carries 1 on both haplotypes
# The chip seed matches 03_metrics.R so the pyramiding chip is fixed between runs

target2_name <- "TargetQTL2"
target2_chr <- 5
target2_pos <- max(SP$genMap[[target2_chr]]) / 2

pyramid_map <- addSegSite(
  tracked_founder_map,
  siteName = target2_name,
  chr = target2_chr,
  mapPos = target2_pos,
  haplo = matrix(c(0L, 0L, 1L, 1L), ncol = 1)
)

SP_pyr <- SimParam$new(pyramid_map)
SP_pyr$setTrackRec(TRUE)

SP_pyr$restrSegSites(excludeSnp = c(target_name, target2_name), overlap = FALSE)

set.seed(founders$seed)
SP_pyr$addSnpChip(nSnpPerChr = 1000, name = "BackgroundChip")

pyramid_founders <- newPop(pyramid_map, simParam = SP_pyr)

fnd_pyramid <- list(
  recurrent = pyramid_founders[1],
  donor = pyramid_founders[2],
  SP = SP_pyr,
  targets = c(target_name, target2_name),
  target_chr = c(target_chr, target2_chr)
)

# Stop if the founders are not 0 and 2 at both targets

founder_target_geno <- pullMarkerGeno(pyramid_founders, markers = fnd_pyramid$targets, simParam = SP_pyr)

stopifnot(all(founder_target_geno[1, ] == 0), all(founder_target_geno[2, ] == 2))


# 2. Run replicated pyramiding simulations
# Foreground selection requires a donor allele at both targets
# 50 replicates, 100 plants per generation, 10 selected, 6 generations

pyramid_schemes <- c("Background", "Dense genomic")
pyramid_seeds <- 3001:3050

pyramid_runs <- pyramid_seeds %>%
  map_dfr(\(s) map_dfr(pyramid_schemes, \(x) run_scheme(x, fnd_pyramid, 100, 10, 6, s))) %>%
  mutate(replicate = match(seed, pyramid_seeds)) %>%
  relocate(replicate, .after = seed)

pyramid_runs %>%
  count(scheme, generation)


# 3. Summarise pyramiding performance
# Drag is reported separately at TargetQTL (linkage_drag_cM) and TargetQTL2 (linkage_drag_cM_2)

pyramid_summary <- pyramid_runs %>%
  group_by(scheme, generation) %>%
  summarise(
    recovery_marker = mean(recovery_marker),
    recovery_ibd = mean(recovery_ibd),
    n_carriers = mean(n_carriers),
    linkage_drag_cM = mean(linkage_drag_cM),
    linkage_drag_cM_2 = mean(linkage_drag_cM_2),
    .groups = "drop"
  )

pyramid_summary

pyramid_thresholds <- pyramid_runs %>%
  group_by(replicate, seed, scheme) %>%
  summarise(gen_99 = first_gen(generation, recovery_ibd, 0.99), .groups = "drop")

pyramid_thresholds %>%
  group_by(scheme) %>%
  summarise(reached_99 = sum(!is.na(gen_99)), mean_gen_99 = mean(gen_99, na.rm = TRUE), .groups = "drop")


# 4. Population size needed for one and two targets
# Observed carrier rates come from the 100-plant single-target sweep and the pyramiding runs
# Required size is the smallest population that gives at least 10 carriers
# with 99% probability under a binomial draw at the expected carrier rate

required_pop <- function(p, n_selected = 10, prob = 0.99) {
  n <- n_selected

  while (pbinom(n_selected - 1, n, p, lower.tail = FALSE) < prob) {
    n <- n + 1
  }

  n
}

single_rate <- readRDS(here("data", "processed", "popsize_sweep.rds")) %>%
  filter(pop_size == 100) %>%
  summarise(rate = mean(n_carriers / pop_size)) %>%
  pull(rate)

pyramid_rate <- pyramid_runs %>%
  summarise(rate = mean(n_carriers / pop_size)) %>%
  pull(rate)

pyramid_popsize <- tibble(
  n_targets = c(1, 2),
  expected_rate = c(0.5, 0.25),
  observed_rate = c(single_rate, pyramid_rate)
) %>%
  mutate(required_pop_size = map_dbl(expected_rate, required_pop))

pyramid_popsize


# 5. Save pyramiding results

saveRDS(pyramid_runs, here("data", "processed", "pyramid_runs.rds"))

saveRDS(pyramid_summary, here("data", "processed", "pyramid_summary.rds"))

saveRDS(pyramid_thresholds, here("data", "processed", "pyramid_thresholds.rds"))

saveRDS(pyramid_popsize, here("data", "processed", "pyramid_popsize.rds"))
