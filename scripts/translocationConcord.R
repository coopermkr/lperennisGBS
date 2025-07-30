#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' L. perennis Seed Exchange Simulation: Concord
#' @date 2025-07-29
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


# Load libraries
library(tidyverse)
library(cgsim)
library(vcfR)
library(adegenet)

#### Test Run
# Load Northeast VCF
nevcf <- read.vcfR("2.stacks/neFiltered.vcf")

# Load full population file and subset to Concord population
conList <- read.delim(file = "2.stacks/nemap.txt", header = FALSE) |>
  rename(id = V1,
         pop = V2) |>
  filter(id %in% colnames(nevcf@gt),
         pop == "CN")

# Subset VCF to only Concord samples and convert to genind
conind <- vcfR2genind(nevcf[,c("FORMAT",conList$id)])

# Average missing data for full set
contab <- tab(conind, freq = TRUE, NA.method = "mean")

# For this first attempt, select 50 random SNPs
lociNames <- sample(locNames(conind),50)

consub <- conind[,loc = lociNames]

# Average missing data
consubtab <- tab(consub, freq = TRUE, NA.method = "mean")

# Load in and subset VCF to AL population
# Load Northeast VCF
nevcf <- read.vcfR("2.stacks/neFiltered.vcf")

# Load full population file and subset to Concord population
alList <- read.delim(file = "2.stacks/nemap.txt", header = FALSE) |>
  rename(id = V1,
         pop = V2) |>
  filter(id %in% colnames(nevcf@gt),
         pop == "AL")

# Subset VCF to only Concord samples and convert to genind
alind <- vcfR2genind(nevcf[,c("FORMAT",alList$id)])

# Average missing data for full set
altab <- tab(alind, freq = TRUE, NA.method = "mean")

# Filter both tab files to only contain the same loci
common <- intersect(colnames(contab), colnames(altab))

contab <- contab[,common]
altab <- altab[,common]

# Run the test simulation
testSim <- effect_sim(pop = consubtab, syears = 5, n = 100)
names(testSim)

save(testSim, file = "4.translocation/testSim.R")
load("4.translocation/testSim.R")
testSim

# Calculate expected and predicted heterozygosity
exp <- gd_ret_exp(out = testSim, ne = "Nc")
exp <- apply(exp, 1, mean)

pred <- gd_ret_pred(out = testSim, gd = "He")
cbind(pred[,1], exp)

#### Predict 100 year diversity in Concord Population

# Code in overlapping generations that increase in fecundity
reprorate <- c(0, 0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1, 1) # Reproduction rates
surrate <- c(0.3, 0.3, rep(0.75, 7), 0) # Pulled from Pavlovic
meth <- "repro_pool"
burn <- 50 # burn-in generations to create stable age distribution
nind <- 70000 # Starting population size
years <- 100 # Number of simulation years
gr <- 0.1 # Population variation rate mean- estimated from census
gm <- "ann"  # growth model
reps <- 5 # Number of replicates to run
nvar <- 0.2 # Variance in initial population size
noff <- 280 # Maximum offspring per individual allowed per year

conDrift <- effect_sim(pop = contab, syears = 1:years, n = nind, n_rep = reps, 
                      growth_model = gm, growth_rate = gr, norm = "norm",
                      nvar = nvar, burnin = burn, noff = noff,
                      mrr = reprorate, frr = reprorate, msr = surrate, fsr = surrate)

save(conDrift, file = "4.translocation/conDrift.Rdata")

## Add random bad years without translocation
bad_years <- sort(sample(1:years, 0.10*years, replace = F)) # Percent of bad years
bad_prop <- 0.1 # Proportion of population removed during bad years

conBad <- effect_sim(pop = contab, syears = 1:years, n = nind, n_rep = reps, 
                    growth_model = gm, growth_rate = gr, norm = "norm",
                    nvar = nvar, burnin = burn, noff = noff,
                    mrr = reprorate, frr = reprorate, msr = surrate, fsr = surrate,
                    # Add bad year parameters
                    cyears = bad_years, cat_prop = bad_prop)

save(conBad, file = "4.translocation/conBad.Rdata")

## Add single translocation from AL in 2006
# Load in and subset VCF to AL population
# Load Northeast VCF
nevcf <- read.vcfR("2.stacks/neFiltered.vcf")

# Set transplant parameters
tsr <- c(0.3, 0.3) # Transplant survival rate- Pulled from Gibbs
trr <- c(0, 0) # Transplant reproduction rates-- assuming seeds
yeffect <- 2 # Number of years transplant survival applies- Pulled from Gibbs/Pavlovic
tr_in <- af(altab) # List of 92-CL allele frequencies
ntran <- 20000 # Number of transplants
ytran <- 4 # Year of augmentation

conTrans <- effect_sim(pop = contab, syears = 1:years, n = nind, n_rep = reps, 
                      growth_model = gm, growth_rate = gr, norm = "norm",
                      nvar = nvar, burnin = burn, noff = noff,
                      mrr = reprorate, frr = reprorate, msr = surrate, fsr = surrate,
                      # Add translocation parameters
                      trans = tr_in, tsim = F, tsr = tsr, tyears = ytran, ntrans = ntran,
                      trans_year = yeffect, trr = trr)

save(conTrans, file = "4.translocation/conTrans.Rdata")

# Add random bad years with translocation
conTransBad <- effect_sim(pop = contab, syears = 1:years, n = nind, n_rep = reps, 
                         growth_model = gm, growth_rate = gr, norm = "norm",
                         nvar = nvar, burnin = burn, noff = noff,
                         mrr = reprorate, frr = reprorate, msr = surrate, fsr = surrate,
                         # Add translocation parameters
                         trans = tr_in, tsim = F, tsr = tsr, tyears = ytran, ntrans = ntran,
                         trans_year = yeffect, trr = trr,
                         # Add bad year parameters
                         cyears = bad_years, cat_prop = bad_prop)

save(conTransBad, file = "4.translocation/conTransBad.Rdata")

