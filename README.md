# Marker-assisted backcross simulation

This repository simulates marker-assisted backcrossing in AlphaSimR on a maize-like genome. It compares four schemes for moving a donor allele into an elite recurrent parent. The comparison covers the generations needed to reach 95% and 99% recurrent parent genome recovery, the donor segment left around the target (linkage drag), and the genotyping cost. It also tests how population size and a second target locus change these results.

The full write-up is in `report.qmd`. All numbers in the report are read from `data/processed/`.

## Schemes

Every scheme backcrosses to the recurrent parent with 100 plants per generation and keeps 10 plants.

- No marker: 10 plants are chosen at random; the target allele is not checked.
- Foreground: 10 plants are chosen at random among plants carrying the donor allele at the target.
- Background: carriers are ranked on the share of recurrent parent alleles at chip SNPs (1,000 per chromosome) that differ between the two founders.
- Dense genomic: carriers are ranked on the same share at every segregating site that differs between the founders.

Dense genomic here means selection on genomic similarity to the recurrent parent at all polymorphic sites. Selection on a genomic relationship matrix or on genomic estimated breeding values is a different definition and was not run.

## Simulation setup

- Genome: `runMacs` maize history, 10 chromosomes of about 2 Morgans, 2,000 segregating sites per chromosome.
- Founders: two doubled haploid lines, genotype 0 (recurrent parent) and 2 (donor) at the target, chosen for the largest SNP distance.
- Target: centre of chromosome 1, excluded from the SNP chip.
- Recovery: marker-estimated (chip) and true identity-by-descent (IBD) values, with recombination tracking on.
- Drag: distance in cM between the nearest loci homozygous for recurrent parent IBD on each side of the target.
- Cost: genotype data points; the target is scored on every plant and polymorphic loci on carriers only.
- Replicates: 50 per scheme; seeds 1001 to 1050 (schemes), 2001 to 2050 (population size), 3001 to 3050 (pyramiding).
- Population size sweep: 50, 100, 200 and 400 plants for Background and Dense genomic, 10 plants kept at each size.
- Pyramiding: a second target at the centre of chromosome 5; carriers must hold both donor alleles.

The population needed to keep 10 carriers with 99% probability is 33 plants for one target and 70 plants for two targets. These two values come from the binomial distribution and do not depend on the simulation seeds.

## Repository structure

```
MAB_introgression_sim/
├─ README.md
├─ report.qmd              # write-up; reads data/processed and outputs/figures
├─ references.bib
├─ _quarto.yml
├─ renv.lock
├─ pro_details.txt         # project design
├─ r_log.txt               # console output of the logged run
├─ baseline_figures/       # output of the unchanged template
├─ r_scripts/
│  ├─ 00_reproduce_temp.R  # Bančič et al. (2025) traitIntrogression.R, unchanged
│  ├─ 01_founders.R        # founder genome, target locus, contrasting founders
│  ├─ 02_schemes.R         # backcross, carrier, selection and cost functions
│  ├─ 03_metrics.R         # tracked founders, IBD recovery, drag, run_scheme()
│  ├─ 04_run_reps.R        # 50 replicates of the four schemes
│  ├─ 05_popsize_sweep.R   # population size and cost
│  ├─ 06_pyramiding.R      # two target loci
│  └─ 07_visualizations.R  # figures
├─ data/
│  ├─ raw/                 # founder_data.rds (founders, SimParam, seed)
│  └─ processed/           # per-generation results and summaries (.rds)
└─ outputs/figures/
```

## How to run

```r
renv::restore()

source("r_scripts/01_founders.R")        # only if founder_data.rds is to be rebuilt
source("r_scripts/04_run_reps.R")
source("r_scripts/05_popsize_sweep.R")
source("r_scripts/06_pyramiding.R")      # needs popsize_sweep.rds from step 05
source("r_scripts/07_visualizations.R")
```

`02_schemes.R` and `03_metrics.R` hold functions and are sourced by scripts 04 to 06. Render the report with `quarto render report.qmd`.

## Reproducibility

Each replicate sets its own seed, and the SNP chip is drawn with the founder seed (50). The run recorded in `r_log.txt` was made before the chip seed was added, so chip-based values (Background and marker-estimated recovery) will differ slightly on re-run. The report reads every value from the output files, so the text follows the new run.

## Not included

- A telomeric target position. A central target is expected to leave more drag, so the drag values partly reflect this choice.
- Two linked targets in repulsion. The population sizes needed grow quickly and the result depends on recombination between the two loci, which needs its own analysis.
- Genotyping error, field costs and polygenic background traits.

## References

Bančič, J., Greenspoon, P., Gaynor, R. C., & Gorjanc, G. (2025). Plant breeding simulations with AlphaSimR. *Crop Science*, 65(1), e21312. https://doi.org/10.1002/csc2.21312

Gaynor, R. C., Gorjanc, G., & Hickey, J. M. (2021). AlphaSimR: an R package for breeding program simulations. *G3 Genes|Genomes|Genetics*, 11(2), jkaa017. https://doi.org/10.1093/g3journal/jkaa017
