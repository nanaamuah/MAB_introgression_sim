library(AlphaSimR)
library(here)


# 1. Generate founder genome
# runMacs gives 10 maize-like chromosomes of about 2 Morgans each
# The seed is saved with the founders

seed <- 50
set.seed(seed)

founder_seqs <- runMacs(nInd = 100, nChr = 10, segSites = 2000, species = "MAIZE")


# 2. Add target locus
# The target sits at the centre of chromosome 1
# The first 50 individuals carry allele 0 and the last 50 carry allele 1

gen_map <- getGenMap(founder_seqs)

chr1_map <- gen_map$pos[gen_map$chr == 1]
target_pos <- mean(range(chr1_map))

target_haplo <- matrix(c(rep(0L, 100), rep(1L, 100)), ncol = 1)

founder_seqs <- addSegSite(
  founder_seqs,
  siteName = "TargetQTL",
  chr = 1,
  mapPos = target_pos,
  haplo = target_haplo
)


# 3. Configure simulation
# TargetQTL is kept off the SNP chip so background selection cannot see it
# The trait is not used later; newPop and makeDH draw phenotypes for it,
# so removing it would change the random draws that produced the saved founders

SP <- SimParam$new(founder_seqs)

SP$restrSegSites(excludeSnp = "TargetQTL", overlap = FALSE)

SP$importTrait(markerNames = "TargetQTL", addEff = 1, name = "TargetTrait", varE = 1)

SP$addSnpChip(nSnpPerChr = 1000, name = "BackgroundChip")


# 4. Create contrasting inbred founders
# Doubled haploids make every line fully inbred
# The recurrent parent has genotype 0 and the donor genotype 2 at TargetQTL
# Among those pairs, the two lines with the largest SNP distance are kept

base_pop <- newPop(founder_seqs, simParam = SP)
inbred_pop <- makeDH(base_pop, simParam = SP)

target_geno <- pullMarkerGeno(inbred_pop, markers = "TargetQTL", simParam = SP)[, 1]

recurrent_candidates <- which(target_geno == 0)
donor_candidates <- which(target_geno == 2)

if (!length(recurrent_candidates) || !length(donor_candidates)) {
  stop("Contrasting founders could not be identified at TargetQTL.")
}

snp_geno <- pullSnpGeno(inbred_pop, simParam = SP)

genomic_distance <- as.matrix(dist(snp_geno, method = "manhattan"))

candidate_distance <- genomic_distance[recurrent_candidates, donor_candidates, drop = FALSE]

best_pair <- which(candidate_distance == max(candidate_distance), arr.ind = TRUE)[1, ]

recurrent_idx <- recurrent_candidates[best_pair[1]]
donor_idx <- donor_candidates[best_pair[2]]

recurrent_parent <- inbred_pop[recurrent_idx]
donor <- inbred_pop[donor_idx]


# 5. Validate and save
# Stop if the founders are not 0 and 2 at TargetQTL or if TargetQTL is on the chip

founders <- c(recurrent_parent, donor)

founder_target_geno <- pullMarkerGeno(founders, markers = "TargetQTL", simParam = SP)[, 1]

stopifnot(
  founder_target_geno[1] == 0,
  founder_target_geno[2] == 2,
  !"TargetQTL" %in% colnames(pullSnpGeno(recurrent_parent, simParam = SP))
)

saveRDS(
  list(
    SP = SP,
    recurrent = recurrent_parent,
    donor = donor,
    founder_seqs = founder_seqs,
    seed = seed,
    target_name = "TargetQTL",
    target_chr = 1,
    target_pos = target_pos
  ),
  here("data", "raw", "founder_data.rds")
)
