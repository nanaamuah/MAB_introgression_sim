library(AlphaSimR)
library(dplyr)
library(tibble)
library(here)

source(here("r_scripts", "02_schemes.R"))

founders <- readRDS(here("data", "raw", "founder_data.rds"))

target_name <- founders$target_name


# 1. Rebuild founders with recombination tracking
# The saved SimParam has no tracking, so the two founders are imported into a new one
# addSnpChip picks SNPs at random; the founder seed fixes the chip between runs

founder_haplo <- pullSegSiteHaplo(c(founders$recurrent, founders$donor), simParam = founders$SP)

gen_map <- getGenMap(founders$founder_seqs)

tracked_founder_map <- importHaplo(haplo = founder_haplo, genMap = gen_map, ploidy = 2)

SP <- SimParam$new(tracked_founder_map)
SP$setTrackRec(TRUE)

SP$restrSegSites(excludeSnp = target_name, overlap = FALSE)

set.seed(founders$seed)
SP$addSnpChip(nSnpPerChr = 1000, name = "BackgroundChip")

tracked_founders <- newPop(tracked_founder_map, simParam = SP)

recurrent <- tracked_founders[1]
donor <- tracked_founders[2]


# 2. True IBD recurrent-parent recovery
# Share of loci on both haplotypes that descend from recurrent-parent founder haplotypes

ibd_recovery <- function(pop, recurrent, SP) {
  rec_ids <- unique(as.vector(pullIbdHaplo(recurrent, simParam = SP)))
  pop_ibd <- pullIbdHaplo(pop, simParam = SP)

  vapply(
    seq_len(nInd(pop)),
    function(i) mean(pop_ibd[c(2 * i - 1, 2 * i), , drop = FALSE] %in% rec_ids),
    numeric(1)
  )
}


# 3. Locate target locus

target_chr <- which(vapply(SP$genMap, function(x) target_name %in% names(x), logical(1)))

if (length(target_chr) != 1) {
  stop("TargetQTL could not be uniquely located.")
}


# 4. Linkage drag
# On each side of the target, the nearest locus homozygous for recurrent-parent IBD marks the segment end
# Drag is the distance between the two ends in cM; a chromosome end is used when no such locus exists

linkage_drag <- function(pop, recurrent, SP, target_chr, target_name) {
  ibd <- pullIbdHaplo(pop, simParam = SP)
  rec_ids <- unique(as.vector(pullIbdHaplo(recurrent, simParam = SP)))

  chr_map <- SP$genMap[[target_chr]]
  chr_names <- names(chr_map)
  target_idx <- which(chr_names == target_name)

  if (length(target_idx) != 1) {
    stop("Target locus could not be uniquely located.")
  }

  chr_cols <- match(chr_names, colnames(ibd))

  vapply(
    seq_len(nInd(pop)),
    function(i) {
      hap <- ibd[c(2 * i - 1, 2 * i), chr_cols, drop = FALSE]

      is_recurrent <- matrix(hap %in% rec_ids, nrow = nrow(hap), ncol = ncol(hap))

      donor_present <- colSums(is_recurrent) < 2

      left_idx <- which(seq_along(chr_map) < target_idx & !donor_present)
      right_idx <- which(seq_along(chr_map) > target_idx & !donor_present)

      left_pos <- if (length(left_idx)) chr_map[max(left_idx)] else min(chr_map)
      right_pos <- if (length(right_idx)) chr_map[min(right_idx)] else max(chr_map)

      as.numeric((right_pos - left_pos) * 100)
    },
    numeric(1)
  )
}


# 5. Run one replicate of a scheme
# The F1 is backcrossed to the recurrent parent for n_generations; selection follows 02_schemes.R
# Recovery, drag and cost are recorded on the selected plants of every generation
# With two targets, drag at the second target is stored as linkage_drag_cM_2

run_scheme <- function(scheme, fnd, pop_size, n_selected, n_generations, seed) {
  set.seed(seed)

  current_pop <- makeCross2(fnd$recurrent, fnd$donor, matrix(c(1, 1), ncol = 2), simParam = fnd$SP)

  results <- vector("list", n_generations)

  for (gen in seq_len(n_generations)) {
    bc_pop <- make_bc(fnd$recurrent, current_pop, pop_size, fnd$SP)

    carriers <- get_carriers(bc_pop, fnd$targets, fnd$SP)

    selected <- select_plants(scheme, bc_pop, carriers, n_selected, fnd$recurrent, fnd$donor, fnd$SP)

    current_pop <- bc_pop[selected]

    marker <- rp_share(current_pop, fnd$recurrent, fnd$donor, fnd$SP)
    ibd <- ibd_recovery(current_pop, fnd$recurrent, fnd$SP)

    res <- tibble(
      scheme = scheme,
      seed = seed,
      generation = gen,
      pop_size = pop_size,
      n_selected = n_selected,
      n_carriers = length(carriers)
    ) %>%
      bind_cols(
        genotyping_cost(
          scheme, pop_size, length(carriers), length(fnd$targets),
          fnd$recurrent, fnd$donor, fnd$SP
        )
      ) %>%
      mutate(
        recovery_marker = mean(marker),
        recovery_ibd = mean(ibd),
        recovery_gap = mean(marker - ibd),
        target_retention = length(get_carriers(current_pop, fnd$targets, fnd$SP)) / nInd(current_pop)
      )

    for (k in seq_along(fnd$targets)) {
      drag_col <- if (k == 1) "linkage_drag_cM" else paste0("linkage_drag_cM_", k)
      drag <- linkage_drag(current_pop, fnd$recurrent, fnd$SP, fnd$target_chr[k], fnd$targets[k])
      res[[drag_col]] <- mean(drag)
    }

    results[[gen]] <- res
  }

  bind_rows(results)
}


# 6. First generation reaching a recovery threshold
# Returns NA when the threshold is not reached within the simulated generations

first_gen <- function(generation, recovery, threshold) {
  hit <- generation[recovery >= threshold]

  if (length(hit)) min(hit) else NA_integer_
}


# 7. Founder set for the single-target runs

fnd_single <- list(
  recurrent = recurrent,
  donor = donor,
  SP = SP,
  targets = target_name,
  target_chr = target_chr
)
