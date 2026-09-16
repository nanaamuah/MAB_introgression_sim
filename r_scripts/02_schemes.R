library(AlphaSimR)
library(dplyr)
library(tibble)


# 1. Backcross and carrier helpers
# Each backcross plant has the recurrent parent as one parent and a random selected plant as the other
# A carrier holds at least one donor allele at every target locus

make_bc <- function(recurrent, current_pop, pop_size, SP) {
  plan <- cbind(rep(1, pop_size), sample(seq_len(nInd(current_pop)), pop_size, replace = TRUE))

  makeCross2(recurrent, current_pop, plan, simParam = SP)
}

get_carriers <- function(pop, target_names, SP) {
  geno <- pullMarkerGeno(pop, markers = target_names, simParam = SP)

  which(rowSums(geno > 0) == length(target_names))
}


# 2. Recurrent-parent allele share
# Only loci that differ between the two founders are used
# dense = FALSE reads the SNP chip; dense = TRUE reads every segregating site

poly_loci <- function(recurrent, donor, SP, dense = FALSE) {
  pull_geno <- if (dense) pullSegSiteGeno else pullSnpGeno

  which(pull_geno(recurrent, simParam = SP)[1, ] != pull_geno(donor, simParam = SP)[1, ])
}

rp_share <- function(pop, recurrent, donor, SP, dense = FALSE) {
  pull_geno <- if (dense) pullSegSiteGeno else pullSnpGeno
  poly <- poly_loci(recurrent, donor, SP, dense)

  if (!length(poly)) {
    stop("No polymorphic loci between founders.")
  }

  rp <- pull_geno(recurrent, simParam = SP)[1, poly]
  geno <- pull_geno(pop, simParam = SP)[, poly, drop = FALSE]

  # Count recurrent-parent alleles at each locus (0, 1 or 2)
  rp_count <- sweep(geno, 2, rp, FUN = function(x, r) ifelse(r == 2, x, 2 - x))

  rowMeans(rp_count / 2)
}


# 3. Selection step for the four schemes
# No marker: random plants, target allele not checked
# Foreground: random plants among carriers
# Background: carriers ranked on recurrent-parent share at polymorphic chip SNPs
# Dense genomic: carriers ranked on recurrent-parent share at all polymorphic segregating sites

select_plants <- function(scheme, bc_pop, carriers, n_selected, recurrent, donor, SP) {
  if (scheme == "No marker") {
    return(sample(seq_len(nInd(bc_pop)), n_selected))
  }

  if (length(carriers) < n_selected) {
    stop("Insufficient target carriers.")
  }

  if (scheme == "Foreground") {
    return(sample(carriers, n_selected))
  }

  if (!scheme %in% c("Background", "Dense genomic")) {
    stop("Unknown scheme.")
  }

  recovery <- rp_share(bc_pop[carriers], recurrent, donor, SP, dense = scheme == "Dense genomic")

  carriers[order(recovery, decreasing = TRUE)[seq_len(n_selected)]]
}


# 4. Genotyping cost per generation
# Foreground scores the target loci on every plant
# Background scores the polymorphic loci on carriers only
# Cost is counted as data points (plants x markers)

genotyping_cost <- function(scheme, pop_size, n_carriers, n_targets, recurrent, donor, SP) {
  n_poly <- if (scheme %in% c("Background", "Dense genomic")) {
    length(poly_loci(recurrent, donor, SP, dense = scheme == "Dense genomic"))
  } else {
    0
  }

  n_fore <- if (scheme == "No marker") 0 else n_targets

  tibble(
    markers_scored = n_fore + n_poly,
    genotyping_cost = pop_size * n_fore + n_carriers * n_poly
  )
}
