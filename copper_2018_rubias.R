#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Copper River Sockeye Inseason with `rubias` 2018 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This inseason analysis serves two purposes:
#   1) Analyze inseason test fishery samples against Tylers 46 population baseline.
#   2) Analyze inseason test fishery samples against the coastwide 473 population baseline.
#      To verify that there is low/no non-local harvest of sockeye.
#   3) Analyze inseason test fishery samples against a master baseline (46 copper + 473 coastwide).
#      To verify that there is low/no non-local harvest of sockeye.
# 
# The baseline update was described in `PWScopper46baseline_rubias.R`


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Load baseline and objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
setwd("V:/Analysis/2_Central/Sockeye/PWSCopper/Copper Inseason 2018/")

library(tidyverse)
library(lubridate)
library(rubias)
library(ggthemes)
library(gridExtra)
library(ggpubr)
# library(conflicted)

source("C:/Users/krshedd//R/Functions.GCL.R")
load_objects(path = "Objects")
load_objects(path = "rubias/baseline")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### SW26 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Find replace 26

sillys <- "SCDVTF18"
sw = 26
sillys_strata <- paste0(sillys, "_SW", sw)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Read in mixture genotypes ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(LocusControl)
CreateLocusControl.GCL(markersuite = "Sockeye2011_96SNPs", username = "krshedd", password = password)
LOKI2R.GCL(sillyvec = sillys, username = "krshedd", password = password)
save_sillys(sillyvec = sillys, path = "Raw genotypes")
# load_sillys(path = "Raw genotypes", sillyvec = sillys)

sapply(sillys, function(silly) {get(paste0(silly, ".gcl"))$n} )  # 190
sapply(sillys, function(silly) {table(get(paste0(silly, ".gcl"))$attributes$CAPTURE_DATE, useNA = "always")} )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Stratify mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tissue_table <- read_csv(file = "Tissue inventory/GEN_SAMPLED_FISH_TISSUE_SW26.csv")
table(tissue_table$CAPTURE_DATE)
n_samp <- 244  # How many tissues were built for this sampling event???

# IDs
samp_dates <- unique(get(paste0(sillys, ".gcl"))$attributes$CAPTURE_DATE)[1]

SW26_IDs <- AttributesToIDs.GCL(silly = sillys, attribute = "CAPTURE_DATE", matching = samp_dates)
SW26_IDs <- list(na.omit(SW26_IDs))
names(SW26_IDs) <- sillys

# Pool
PoolCollections.GCL(collections = sillys, loci = loci96, IDs = SW26_IDs, newname = sillys_strata)
SCDVTF18_SW26.gcl$n ## 190
table(SCDVTF18_SW26.gcl$attributes$CAPTURE_DATE)

strata_date_range <- unique(SCDVTF18_SW26.gcl$attributes$CAPTURE_DATE)

# Save
# dir.create("Raw genotypes/Strata")
save_sillys(sillyvec = sillys_strata, path = "Raw genotypes/Strata/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sillys_strata_2018_SW26_n <- matrix(data = NA, nrow = length(sillys_strata), ncol = 4, 
                                    dimnames = list(sillys_strata, c("Genotyped", "Missing", "Duplicate", "Final")))

#### Check loci
## Get sample size by locus
original_sillys_strata_2018_SW26_n_locus <- SampSizeByLocus.GCL(sillyvec = sillys_strata, loci = loci96)
min(original_sillys_strata_2018_SW26_n_locus)  ## 167
round(apply(original_sillys_strata_2018_SW26_n_locus, 1, function(locus) {min(locus) / max(locus)}), 2)  # 0.88

original_sillys_strata_percent_locus <- apply(original_sillys_strata_2018_SW26_n_locus, 1, function(row) {row / max(row)} )  # 0.89
which(apply(original_sillys_strata_percent_locus, 2, min) < 0.8)  # no re-runs!

# Genotpying percentage by locus and strata
require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(original_sillys_strata_percent_locus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares

#### Check individuals
### Initial
## Get number of individuals per silly before removing missing loci individuals
sillys_strata_2018_SW26_n[, "Genotyped"] <- sapply(paste0(sillys_strata, ".gcl"), function(x) get(x)$n)

### Missing
## Remove individuals with >20% missing data
sillys_strata_2018_SW26_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = sillys_strata, proportion = 0.8)
save_objects(objects = "sillys_strata_2018_SW26_MissLoci", path = "Objects")

## Get number of individuals per silly after removing missing loci individuals
sillys_strata_2018_SW26_n[, "Missing"] <- sillys_strata_2018_SW26_n[, "Genotyped"] - 
  sapply(paste0(sillys_strata, ".gcl"), function(x) get(x)$n)

### Duplicate
## Check within collections for duplicate individuals.
sillys_strata_2018_SW26_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = sillys_strata, loci = loci96, quantile = NULL, minproportion = 0.95)
detach("package:reshape", unload = TRUE); library(tidyverse)
sillys_strata_2018_SW26_DuplicateCheckReportSummary <- sapply(sillys_strata, function(x) sillys_strata_2018_SW26_DuplicateCheck95MinProportion[[x]]$report)
sillys_strata_2018_SW26_DuplicateCheckReportSummary
save_objects(objects = "sillys_strata_2018_SW26_DuplicateCheckReportSummary", path = "Objects")

## Remove duplicate individuals
sillys_strata_2018_SW26_RemovedDups <- RemoveDups.GCL(sillys_strata_2018_SW26_DuplicateCheck95MinProportion)

### Final
sillys_strata_2018_SW26_n[, "Final"] <- sapply(paste0(sillys_strata, ".gcl"), function(x) get(x)$n)
## Get number of individuals per silly after removing duplicate individuals
sillys_strata_2018_SW26_n[, "Duplicate"] <- sillys_strata_2018_SW26_n[, "Genotyped"] - sillys_strata_2018_SW26_n[, "Missing"] - sillys_strata_2018_SW26_n[, "Final"]
sillys_strata_2018_SW26_n

save_objects(objects = "sillys_strata_2018_SW26_n", path = "Objects")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Combine loci

# for loci91 for PWSCopper46Pops
CombineLoci.GCL(sillyvec = sillys_strata, markerset = c("One_MHC2_251", "One_MHC2_190"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = sillys_strata, markerset = c("One_Cytb_26", "One_CO1", "One_Cytb_17"), delim=".", update = TRUE)

# for loci91_v2 for Coastwide482Pops
CombineLoci.GCL(sillyvec = sillys_strata, markerset = c("One_MHC2_190", "One_MHC2_251"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = sillys_strata, markerset = c("One_CO1", "One_Cytb_17", "One_Cytb_26"), delim=".", update = TRUE)

# for loci89 for KMA473Pops
CombineLoci.GCL(sillyvec = sillys_strata, markerset = c("One_Tf_ex10-750", "One_Tf_ex3-182"), delim = ".", update = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### KMA473 MSA with rubias ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create mixtures
dir.create("rubias/mixture/kma473")
copper_2018_SW26_loci89.mix <- create_rubias_mixture(sillyvec = sillys_strata, loci = loci89, path = "rubias/mixture/kma473")
save_objects(objects = "copper_2018_SW26_loci89.mix", path = "rubias/mixture/kma473")
# load_objects(path = "rubias/mixture")
str(copper_2018_SW26_loci89.mix, give.attr = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Run MSA
dir.create("rubias/output/kma473")
copper_2018_SW26_loci89.out <- run_rubias_mixture(reference = kma473_loci89.base,
                                                  mixture = copper_2018_SW26_loci89.mix, 
                                                  group_names = Groups15,
                                                  gen_start_col = 5, 
                                                  method = "MCMC", 
                                                  reps = 25000, 
                                                  burn_in = 5000, 
                                                  pb_iter = 100,
                                                  sample_int_Pi = 10,
                                                  path = "rubias/output/kma473")
str(copper_2018_SW26_loci89.out, give.attr = FALSE, max.level = 2)
copper_2018_SW26_loci89.out$mixing_proportions %>% 
  group_by(repunit) %>% 
  summarize(rho = sum(pi))  # Copper = 99%


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Coastwide482 MSA with rubias ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create mixtures
dir.create("rubias/mixture/coastwide482")
copper_2018_SW26_loci91_v2.mix <- create_rubias_mixture(sillyvec = sillys_strata, loci = loci91_v2, path = "rubias/mixture/coastwide482")
save_objects(objects = "copper_2018_SW26_loci91_v2.mix", path = "rubias/mixture/coastwide482")
# load_objects(path = "rubias/mixture")
str(copper_2018_SW26_loci91_v2.mix, give.attr = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Run MSA
dir.create("rubias/output/coastwide482")
copper_2018_SW26_loci91_v2.out <- run_rubias_mixture(reference = copper482_loci91_v2.base,
                                                     mixture = copper_2018_SW26_loci91_v2.mix, 
                                                     group_names = PWSCopper9Groups_pub,
                                                     gen_start_col = 5, 
                                                     method = "MCMC", 
                                                     reps = 25000, 
                                                     burn_in = 5000, 
                                                     pb_iter = 100,
                                                     sample_int_Pi = 10,
                                                     path = "rubias/output/coastwide482")
str(copper_2018_SW26_loci91_v2.out, give.attr = FALSE, max.level = 2)
copper_2018_SW26_loci91_v2.out$mixing_proportions %>% 
  group_by(repunit) %>% 
  summarize(rho = sum(pi))  # Other < 1%

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### MSA with rubias ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create mixtures
# dir.create("rubias/mixture")
copper_2018_SW26.mix <- create_rubias_mixture(sillyvec = sillys_strata, loci = loci91, path = "rubias/mixture")
save_objects(objects = "copper_2018_SW26.mix", path = "rubias/mixture")
# load_objects(path = "rubias/mixture")
str(copper_2018_SW26.mix, give.attr = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Run MSA
copper_2018_SW26.out <- run_rubias_mixture(reference = copper46_loci91.base,
                                           mixture = copper_2018_SW26.mix, 
                                           group_names = PWSCopper8Groups_pub,
                                           gen_start_col = 5, 
                                           method = "MCMC", 
                                           reps = 25000, 
                                           burn_in = 5000, 
                                           pb_iter = 100,
                                           sample_int_Pi = 10,
                                           path = "rubias/output")
str(copper_2018_SW26.out, give.attr = FALSE, max.level = 2)
# save.image("copper_2018_SW26.RData")
# load("rubias/output/copper_2018_SW26.RData")
copper_2018_SW26.out$mixing_proportions %>% 
  group_by(repunit) %>% 
  summarize(rho = sum(pi))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## rubais output
copper_2018_SW26.sum <- custom_combine_rubias_output(rubias_output = copper_2018_SW26.out, group_names = PWSCopper8Groups_pub, bias_corr = FALSE)
str(copper_2018_SW26.sum)

# Save as R object via `dput`
# dir.create("Estimates objects")

# Write out as .csv file
# dir.create("Estimates tables")
copper_2018_SW26_dates.sum <- copper_2018_SW26.sum %>%
  separate(col = mixture_collection, into = c("year", "stat_week"), sep = "_SW") %>%  # split mixture into silly and stat_week
  separate(col = year, into = c("trash", "year"), sep = "SCDVTF") %>%
  select(-trash) %>%
  mutate(year = as.integer(year) + 2000) %>%  # year as integer
  mutate(sample_date = paste(month(strata_date_range), day(strata_date_range), sep = "/")) %>%  # sample date
  rename(reporting_group = repunit) %>%  # rename to reporting group
  select(year, stat_week, sample_date, reporting_group, mean, sd, median, `5%`, `95%`, `P=0`)  # order


save_objects(objects = "copper_2018_SW26_dates.sum", path = "Estimates objects")
write_csv(x = copper_2018_SW26_dates.sum, path = "Estimates tables/CopperRiver_Sockeye_2018_SW26_Estimates.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create report ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dget data from all SWs
copper_2018_dates.sum <- bind_rows(lapply(list.files(path = "Estimates objects", pattern = "2018_SW", full.names = TRUE), dget)) %>% 
  mutate(stat_week = factor(x = stat_week, levels = unique(stat_week)))
save_objects(objects = "copper_2018_dates.sum", path = "Estimates objects")

# Save all output for Rmd
save.image("copper_2018_SW26.RData")
