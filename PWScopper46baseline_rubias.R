# Quick look with `rubias` to determine if PWS/Copper groups are highly identifiable and see if can use Chignik/Port Moller 24 SNPs
load("V:/Analysis/2_Central/Sockeye/PWSCopper/NFWF baseline project/PWSCopper2013SockeyeBaseline.RData")
require(tidyverse)
require(rubias)
source("C:/Users/krshedd/R/Functions.GCL.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Baseline update ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tyler created a 62 collection, 46 population PWS/Copper River baseline
# Here I will do 2 main baseline cleanup tasks
# 1) Create a rubias version of Tyler's baseline with associated objects
# 2) Test the baseline for 96 SNPs, and also 24 SNPs used for Chignik and Port Moller to see if we can reduce the markerset


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Set up directory structure
setwd("V:/Analysis/2_Central/Sockeye/PWSCopper/Copper Inseason 2018/")
dirs <- c("Baseline genotypes", "Estimates objects", "Estimates tables", "Copper baseline test", "Objects", "Raw genotypes", "rubias", "Updates")
sapply(dirs, dir.create)

rubias_dirs <- c("baseline", "mixture", "output")
setwd("V:/Analysis/2_Central/Sockeye/PWSCopper/Copper Inseason 2018/rubias/")
sapply(rubias_dirs, dir.create)
setwd("V:/Analysis/2_Central/Sockeye/PWSCopper/Copper Inseason 2018/")

## Save baseline objects
baseline_objects <- c("LocusControl", "loci91", "loci96", "PWSCopper62Collections", "PWSCopper46Pops", "CommonNames46", "PWSCopper2013SockeyeSampleSizes", "PWSCopper8Groups", "PWSCopper46GroupVec2")
save_objects(objects = baseline_objects, path = "Objects")

colors8 <- readClipboard()  # V:\Analysis\2_Central\Sockeye\PWSCopper\NFWF baseline project\PWS Copper River Sockeye Analysis.xlsx, tab "Names", column Q
save_objects(objects = "colors8", path = "Objects")

PWSCopper8Groups_pub <- c("Southwest PWS", "Main Bay/North PWS", "Copper  Delta", "Upper Copper", "Gulkana", "Tazlina and Outlets", "Klutina", "Lower Copper")
save_objects(objects = "PWSCopper8Groups_pub", path = "Objects")

## Save sillys
save_sillys(sillyvec = PWSCopper46Pops, path = "Baseline genotypes")

## Get Port Moller and Chignik 24 loci sets
loci24_chignik <- dget(file = "V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2018/Baseline 2012/Objects/loci24.txt")
loci24_port_moller <- dget(file = "V:/Analysis/2_Central/Sockeye/Bristol Bay/2013 Baseline/Objects/loci24Theta2.txt")
save_objects(objects = c("loci24_chignik", "loci24_port_moller"), path = "Objects")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Load objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

library(tidyverse)
library(rubias)
source("C:/Users/krshedd/R/Functions.GCL.R")

setwd("V:/Analysis/2_Central/Sockeye/PWSCopper/Copper Inseason 2018/")
load_objects(path = "Objects")
load_sillys(path = "Baseline genotypes/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### rubias testing ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create baseline
# loci91
copper46_loci91.base <- create_rubias_baseline(sillyvec = PWSCopper46Pops, loci = loci91, group_names = PWSCopper8Groups_pub, groupvec = PWSCopper46GroupVec2, baseline_name = "copper46_loci91")
save_objects(objects = "copper46_loci91.base", path = "rubias/baseline/")
load_objects(path = "rubias/baseline/")

# loci24_port_moller
copper46_loci24_port_moller.base <- create_rubias_baseline(sillyvec = PWSCopper46Pops, loci = loci24_port_moller, group_names = PWSCopper8Groups_pub, groupvec = PWSCopper46GroupVec2, baseline_name = "copper46_loci24_port_moller")
save_objects(objects = "copper46_loci24_port_moller.base", path = "rubias/baseline/")
# load_objects(path = "rubias/baseline/")

# loci24_chignik_mod
loci24_chignik_mod <- gsub(pattern = "One_GPDH2.One_GPDH", replacement = "One_GPDH", x = loci24_chignik)
save_objects(objects = "loci24_chignik_mod", path = "Objects")
copper46_loci24_chignik_mod.base <- create_rubias_baseline(sillyvec = PWSCopper46Pops, loci = loci24_chignik_mod, group_names = PWSCopper8Groups_pub, groupvec = PWSCopper46GroupVec2, baseline_name = "copper46_loci24_chignik_mod")
save_objects(objects = "copper46_loci24_chignik_mod.base", path = "rubias/baseline/")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Leave one out testing ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir.create("Baseline test results")

## loci91
copper46_loci91.base_loo <- assess_reference_loo(reference = copper46_loci91.base, gen_start_col = 5, reps = 100, mixsize = 200)

# Summarize to reporting unit level
loo_loci91_out <- copper46_loci91.base_loo %>% 
  mutate(repunit_f = factor(x = repunit, levels = PWSCopper8Groups_pub)) %>% 
  group_by(repunit_scenario, iter, repunit_f) %>% 
  summarise(true_repprop = sum(true_pi), repprop_posterior_mean = sum(post_mean_pi), repu_n = sum(n)) %>% 
  mutate(repu_n_prop = repu_n / sum(repu_n))

# Plot LOO results
ggplot(loo_loci91_out, aes(x = repu_n_prop, y = repprop_posterior_mean, colour = repunit_f)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  scale_color_manual(name = "Reporting Group", values = colors8) +
  facet_wrap(~ repunit_f) +
  xlab("True Reporting Group Proportion") +
  ylab("Posterior Mean Reporting Group Proportion") +
  ggtitle("Lynn Canal Leave-one-out Test Results 91 loci")
ggsave(filename = "Baseline test results/Leave-one-out_loci91.png", device = "png", width = 6.5, height = 6.5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## loci24_port_moller
copper46_loci24_port_moller.base_loo <- assess_reference_loo(reference = copper46_loci24_port_moller.base, gen_start_col = 5, reps = 100, mixsize = 200)

# Summarize to reporting unit level
loo_loci24_port_moller_out <- copper46_loci24_port_moller.base_loo %>% 
  mutate(repunit_f = factor(x = repunit, levels = PWSCopper8Groups_pub)) %>% 
  group_by(repunit_scenario, iter, repunit_f) %>% 
  summarise(true_repprop = sum(true_pi), repprop_posterior_mean = sum(post_mean_pi), repu_n = sum(n)) %>% 
  mutate(repu_n_prop = repu_n / sum(repu_n))

# Plot LOO results
ggplot(loo_loci24_port_moller_out, aes(x = repu_n_prop, y = repprop_posterior_mean, colour = repunit_f)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  scale_color_manual(name = "Reporting Group", values = colors8) +
  facet_wrap(~ repunit_f) +
  xlab("True Reporting Group Proportion") +
  ylab("Posterior Mean Reporting Group Proportion") +
  ggtitle("Lynn Canal Leave-one-out Test Results 24_port_moller loci")
ggsave(filename = "Baseline test results/Leave-one-out_loci24_port_moller.png", device = "png", width = 6.5, height = 6.5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## loci24_chignik_mod
copper46_loci24_chignik_mod.base_loo <- assess_reference_loo(reference = copper46_loci24_chignik_mod.base, gen_start_col = 5, reps = 100, mixsize = 200)

# Summarize to reporting unit level
loo_loci24_chignik_mod_out <- copper46_loci24_chignik_mod.base_loo %>% 
  mutate(repunit_f = factor(x = repunit, levels = PWSCopper8Groups_pub)) %>% 
  group_by(repunit_scenario, iter, repunit_f) %>% 
  summarise(true_repprop = sum(true_pi), repprop_posterior_mean = sum(post_mean_pi), repu_n = sum(n)) %>% 
  mutate(repu_n_prop = repu_n / sum(repu_n))

# Plot LOO results
ggplot(loo_loci24_chignik_mod_out, aes(x = repu_n_prop, y = repprop_posterior_mean, colour = repunit_f)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  scale_color_manual(name = "Reporting Group", values = colors8) +
  facet_wrap(~ repunit_f) +
  xlab("True Reporting Group Proportion") +
  ylab("Posterior Mean Reporting Group Proportion") +
  ggtitle("Lynn Canal Leave-one-out Test Results 24_chignik_mod loci")
ggsave(filename = "Baseline test results/Leave-one-out_loci24_chignik_mod.png", device = "png", width = 6.5, height = 6.5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Decision point ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Use all 96 SNPs in order to maintain high accuracy and precision
# Kyle Shedd Fri Jun 29 15:46:18 2018


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Combine Tyler's 2013 PWS/Copper with Coastwide 473 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
require(tidyverse)
require(rubias)
source("C:/Users/krshedd/R/Functions.GCL.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get coastwide silly's
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/")

#~~~~~~~~~~~~~~~~~~
# Get Locus Control, objects, and populations
LocusControl <- dget(file = "Objects/LocusControl103.txt")
loci89 <- dget(file = "Objects/loci89.txt")
KMA473Pops <- dget(file = "Objects/KMA473Pops.txt")
KMA473PopsGroupVec15 <- dget(file = "Objects/KMA473PopsGroupVec15.txt")
Groups15 <- dget(file = "Objects/Groups15.txt")
Colors15 <- dget(file = "Objects/Colors15.txt")

load_sillys(path = "Raw genotypes/PostQCPopsloci103", sillyvec = KMA473Pops)
length(objects(pattern = "\\.gcl"))

#~~~~~~~~~~~~~~~~~~
# Create rubias baseline of KMA473
setwd("V:/Analysis/2_Central/Sockeye/PWSCopper/Copper Inseason 2018/")
kma473_loci89.base <- create_rubias_baseline(sillyvec = KMA473Pops, loci = loci89, group_names = Groups15, groupvec = KMA473PopsGroupVec15, baseline_name = "kma473_loci89")
save_objects(objects = "kma473_loci89.base", path = "rubias/baseline/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## rubias loo tests
# kma 473 loci89
kma473_loci89.base_loo <- assess_reference_loo(reference = kma473_loci89.base, gen_start_col = 5, reps = 100, mixsize = 200)

# Summarize to reporting unit level
loo_loci91_v2_out <- kma473_loci89.base_loo %>% 
  mutate(repunit_f = factor(x = repunit, levels = Groups15)) %>% 
  group_by(repunit_scenario, iter, repunit_f) %>% 
  summarise(true_repprop = sum(true_pi), repprop_posterior_mean = sum(post_mean_pi), repu_n = sum(n)) %>% 
  mutate(repu_n_prop = repu_n / sum(repu_n))

# Plot LOO results
ggplot(loo_loci91_v2_out, aes(x = repu_n_prop, y = repprop_posterior_mean, colour = repunit_f)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  scale_color_manual(name = "Reporting Group", values = Colors15) +
  facet_wrap(~ repunit_f) +
  xlab("True Reporting Group Proportion") +
  ylab("Posterior Mean Reporting Group Proportion") +
  ggtitle("Lynn Canal Leave-one-out Test Results coastwide_loci91_v2 loci")
ggsave(filename = "Baseline test results/Leave-one-out_kma_loci89.png", device = "png", width = 6.5, height = 6.5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Drop PWS/Copper sillys
table(KMA473PopsGroupVec15)  # remove 37 PWS Copper sillys (group 14)
sum(KMA473PopsGroupVec15 == which(Groups15 == "PWS Copper"))
Coastwide436Pops <- KMA473Pops[which(KMA473PopsGroupVec15 != which(Groups15 == "PWS Copper"))]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get PWS/Copper silly's
setwd("V:/Analysis/2_Central/Sockeye/PWSCopper/Copper Inseason 2018/")

load_objects(path = "Objects")
load_sillys(path = "Baseline genotypes")

Coastwide482Pops <- c(PWSCopper46Pops, Coastwide436Pops)
PWSCopper482GroupVec <- c(PWSCopper46GroupVec2, rep(9, 436))
PWSCopper9Groups_pub <- c(PWSCopper8Groups_pub, "Other")
colors9 <- c(colors8, "grey20")

CombineLoci.GCL(sillyvec = PWSCopper46Pops, markerset = c("One_MHC2_190", "One_MHC2_251"), delim = ".", update = TRUE)
CombineLoci.GCL(sillyvec = PWSCopper46Pops, markerset = c("One_CO1", "One_Cytb_17", "One_Cytb_26"), delim=".", update = TRUE)

loci91_v2 <- c(loci91[-c(90:91)], "One_MHC2_190.One_MHC2_251", "One_CO1.One_Cytb_17.One_Cytb_26")
save_objects(objects = c("Coastwide482Pops", "PWSCopper482GroupVec", "PWSCopper9Groups_pub", "loci91_v2", "colors9", "loci89", "Groups15"), path = "Objects")

copper482_loci91_v2.base <- create_rubias_baseline(sillyvec = Coastwide482Pops, loci = loci91_v2, group_names = PWSCopper9Groups_pub, groupvec = PWSCopper482GroupVec, baseline_name = "copper482_loci91_v2")
save_objects(objects = "copper482_loci91_v2.base", path = "rubias/baseline/")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## rubias loo tests
# copper 482 loci91_v2
copper482_loci91_v2.base_loo <- assess_reference_loo(reference = copper482_loci91_v2.base, gen_start_col = 5, reps = 100, mixsize = 200)

# Summarize to reporting unit level
loo_loci91_v2_out <- copper482_loci91_v2.base_loo %>% 
  mutate(repunit_f = factor(x = repunit, levels = PWSCopper9Groups_pub)) %>% 
  group_by(repunit_scenario, iter, repunit_f) %>% 
  summarise(true_repprop = sum(true_pi), repprop_posterior_mean = sum(post_mean_pi), repu_n = sum(n)) %>% 
  mutate(repu_n_prop = repu_n / sum(repu_n))

# Plot LOO results
ggplot(loo_loci91_v2_out, aes(x = repu_n_prop, y = repprop_posterior_mean, colour = repunit_f)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  scale_color_manual(name = "Reporting Group", values = colors9) +
  facet_wrap(~ repunit_f) +
  xlab("True Reporting Group Proportion") +
  ylab("Posterior Mean Reporting Group Proportion") +
  ggtitle("Lynn Canal Leave-one-out Test Results coastwide_loci91_v2 loci")
ggsave(filename = "Baseline test results/Leave-one-out_coastwide_loci91_v2.png", device = "png", width = 6.5, height = 6.5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### New Groups, split Mendeltna ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

library(tidyverse)
library(rubias)
source("C:/Users/krshedd/R/Functions.GCL.R")

setwd("V:/Analysis/2_Central/Sockeye/PWSCopper/Copper Inseason 2018/")
load_objects(path = "Objects")
load_sillys(path = "Baseline genotypes/")

sapply(PWSCopper8Groups_pub, function(grp) {i <- which(PWSCopper8Groups_pub == grp); PWSCopper46Pops[PWSCopper46GroupVec2 == i]})

PWSCopper9Groups_46Pops_pub <- c(PWSCopper8Groups_pub[1:5], "Tazlina", "Klutina/Tonsina Outlets", "Klutina Lake", PWSCopper8Groups_pub[8])
PWSCopper46GroupVec9 <- as.numeric(readClipboard())  # V:\Analysis\2_Central\Sockeye\PWSCopper\NFWF baseline project\PWS Copper River Sockeye Analysis.xlsx; tab Names; column AC

sapply(PWSCopper9Groups_46Pops_pub, function(grp) {i <- which(PWSCopper9Groups_46Pops_pub == grp); PWSCopper46Pops[PWSCopper46GroupVec9 == i]})
cbind("Silly" = PWSCopper46Pops, "old" = PWSCopper8Groups_pub[PWSCopper46GroupVec2], "new" = PWSCopper9Groups_46Pops_pub[PWSCopper46GroupVec9])

colors9 <- c(colors8[1:5], "grey20", "grey60", colors8[7:8])

save_objects(objects = c("PWSCopper9Groups_46Pops_pub", "PWSCopper46GroupVec9", "colors9"), path = "Objects")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create baseline
# loci91
copper46_9groups_loci91.base <- create_rubias_baseline(sillyvec = PWSCopper46Pops, loci = loci91, group_names = PWSCopper9Groups_46Pops_pub, groupvec = PWSCopper46GroupVec9, baseline_name = "copper46_9groups_loci91")
save_objects(objects = "copper46_9groups_loci91.base", path = "rubias/baseline/")

## loci91
copper46_9groups_loci91.base_loo <- assess_reference_loo(reference = copper46_9groups_loci91.base, gen_start_col = 5, reps = 100, mixsize = 200)

# Summarize to reporting unit level
loo_loci91_out <- copper46_9groups_loci91.base_loo %>% 
  mutate(repunit_f = factor(x = repunit, levels = PWSCopper9Groups_46Pops_pub)) %>% 
  group_by(repunit_scenario, iter, repunit_f) %>% 
  summarise(true_repprop = sum(true_pi), repprop_posterior_mean = sum(post_mean_pi), repu_n = sum(n)) %>% 
  mutate(repu_n_prop = repu_n / sum(repu_n))

# Plot LOO results
ggplot(loo_loci91_out, aes(x = repu_n_prop, y = repprop_posterior_mean, colour = repunit_f)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  scale_color_manual(name = "Reporting Group", values = colors9) +
  facet_wrap(~ repunit_f) +
  xlab("True Reporting Group Proportion") +
  ylab("Posterior Mean Reporting Group Proportion") +
  ggtitle("Lynn Canal Leave-one-out Test Results 91 loci")
ggsave(filename = "Baseline test results/Leave-one-out_9groups_loci91.png", device = "png", width = 6.5, height = 6.5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Map
ncol_per_pop <- sapply(strsplit(x = PWSCopper46Pops, "\\."), length)
PWSCopper60GroupVec9 <- rep(x = PWSCopper46GroupVec9, times = ncol_per_pop)
writeClipboard(PWSCopper60GroupVec9)
writeClipboard(PWSCopper9Groups_46Pops_pub[PWSCopper60GroupVec9])
writeClipboard(colors9[PWSCopper60GroupVec9])
write.table(x = t(col2rgb(colors9[PWSCopper60GroupVec9])), file = 'clipboard', col.names = FALSE, row.names = FALSE)
t(col2rgb(colors9))

plot.new()
legend("topleft", legend = PWSCopper9Groups_46Pops_pub, fill = colors9, bty = "n", cex = 1.5, title = "Reporting Groups")
