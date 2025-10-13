#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' L. perennis Seed Exchange Simulation: Vermont
#' @date 2025-07-29
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


# Load libraries
library(tidyverse)
library(cgsim)
library(vcfR)
library(adegenet)

#### Population Establishment
# Load Northeast VCF
nevcf <- read.vcfR("2.stacks/neFiltered.vcf")

# Load full population file and subset to Vermont population
vtList <- read.delim(file = "2.stacks/nemap.txt", header = FALSE) |>
  rename(id = V1,
         pop = V2) |>
  filter(id %in% colnames(nevcf@gt),
         pop == "CL")

# Subset VCF to only Vermont 2022 samples and convert to genind
vtind <- vcfR2genind(nevcf[,c("FORMAT",vtList$id)])

# Average missing data for full set
vttab <- tab(vtind, freq = TRUE, NA.method = "mean")

# Load full population file and subset to Albany population
alList <- read.delim(file = "2.stacks/nemap.txt", header = FALSE) |>
  rename(id = V1,
         pop = V2) |>
  filter(id %in% colnames(nevcf@gt),
         pop == "92-CL")

# Subset VCF to only Vermont 1992 samples and convert to genind
alInd <- vcfR2genind(nevcf[,c("FORMAT",alList$id)])

# Average missing data for full set
altab <- tab(alInd, freq = TRUE, NA.method = "mean")

# Filter both tab files to only contain the same loci
common <- intersect(colnames(vttab), colnames(altab))

# Subset the original VT loci to only include those in the AL set
vttab <- vttab[,common]

# Load full population file and subset to Vermont 1992 population
vt92List <- read.delim(file = "2.stacks/nemap.txt", header = FALSE) |>
  rename(id = V1,
         pop = V2) |>
  filter(id %in% colnames(nevcf@gt),
         pop == "92-CL")

# Subset VCF to only Vermont 1992 samples and convert to genind
vt92ind <- vcfR2genind(nevcf[,c("FORMAT",vt92List$id)])

# Average missing data for full set
vt92tab <- tab(vt92ind, freq = TRUE, NA.method = "mean")

# Filter both tab files to only contain the same loci
vtcommon <- intersect(colnames(vttab), colnames(vt92tab))

# Now subset all three so they share only the same loci
vttab <- vttab[,vtcommon]
vt92tab <- vt92tab[,vtcommon]
altab <- altab[,common]
# There are 3738 loci between all three populations

#### Predict 100 year diversity in Vermont Population

# Code in overlapping generations that increase in fecundity
reprorate <- c(0, 0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1, 1) # Reproduction rates
surrate <- c(0.3, 0.3, rep(0.75, 7), 0) # Pulled from Pavlovic
meth <- "repro_pool"
burn <- 100 # burn-in generations to create stable age distribution
nind <- 44 # Starting population size
years <- 100 # Number of simulation years
gr <- c(0.0125,.001) # Population variation rate mean- estimated from census
gm <- "exp"  # growth model
reps <- 5 # Number of replicates to run
nvar <- 0.15 # Variance in initial population size
noff <- 280 # Maximum offspring per individual allowed per year

vtDrift <- effect_sim(pop = vttab, syears = 1:years, n = nind, n_rep = reps, 
                      growth_model = gm, growth_rate = gr, norm = "norm",
                      nvar = nvar, burnin = burn, noff = noff,
                      mrr = reprorate, frr = reprorate, msr = surrate, fsr = surrate)

save(vtDrift, file = "4.translocation/vtDrift.Rdata")

## Add random bad years without translocation
bad_years <- sort(sample(1:years, 0.10*years, replace = F)) # Percent of bad years
bad_prop <- 0.1 # Proportion of population removed during bad years

vtBad <- effect_sim(pop = vttab, syears = 1:years, n = nind, n_rep = reps, 
                    growth_model = gm, growth_rate = gr, norm = "norm",
                    nvar = nvar, burnin = burn, noff = noff,
                    mrr = reprorate, frr = reprorate, msr = surrate, fsr = surrate,
                    # Add bad year parameters
                    cyears = bad_years, cat_prop = bad_prop)

save(vtBad, file = "4.translocation/vtBad.Rdata")

## Add single translocation from 92-VT in 2006
# Set transplant parameters
tsr <- c(1, 0.3) # Transplant survival rate- Pulled from Gibbs
trr <- c(0, 0) # Transplant reproduction rates-- assuming seeds
yeffect <- 2 # Number of years transplant survival applies- Pulled from Gibbs/Pavlovic
tr_in <- af(vt92tab) # List of 92-CL allele frequencies
ntran <- 24 # Number of transplants
ytran <- 4 # Year of augmentation

vt92Trans <- effect_sim(pop = vttab, syears = 1:years, n = nind, n_rep = reps, 
                      growth_model = gm, growth_rate = gr, norm = "norm",
                      nvar = nvar, burnin = burn, noff = noff,
                      mrr = reprorate, frr = reprorate, msr = surrate, fsr = surrate,
                      # Add translocation parameters
                      trans = tr_in, tsim = F, tsr = tsr, tyears = ytran, ntrans = ntran,
                      trans_year = yeffect, trr = trr)

save(vt92Trans, file = "4.translocation/vt92Trans.Rdata")

# Add random bad years with translocation
vt92TransBad <- effect_sim(pop = vttab, syears = 1:years, n = nind, n_rep = reps, 
                         growth_model = gm, growth_rate = gr, norm = "norm",
                         nvar = nvar, burnin = burn, noff = noff,
                         mrr = reprorate, frr = reprorate, msr = surrate, fsr = surrate,
                         # Add translocation parameters
                         trans = tr_in, tsim = F, tsr = tsr, tyears = ytran, ntrans = ntran,
                         trans_year = yeffect, trr = trr,
                         # Add bad year parameters
                         cyears = bad_years, cat_prop = bad_prop)

save(vt92TransBad, file = "4.translocation/vt92TransBad.Rdata")

## Reset transplant source to Albany, keep all other params the same
# Set transplant parameters
tr_in <- af(altab) # List of 92-CL allele frequencies
vtalTrans <- effect_sim(pop = vttab, syears = 1:years, n = nind, n_rep = reps, 
                      growth_model = gm, growth_rate = gr, norm = "norm",
                      nvar = nvar, burnin = burn, noff = noff,
                      mrr = reprorate, frr = reprorate, msr = surrate, fsr = surrate,
                      # Add translocation parameters
                      trans = tr_in, tsim = F, tsr = tsr, tyears = ytran, ntrans = ntran,
                      trans_year = yeffect, trr = trr)

save(vtalTrans, file = "4.translocation/vtalTrans.Rdata")

# Add random bad years with translocation
vtalTransBad <- effect_sim(pop = vttab, syears = 1:years, n = nind, n_rep = reps, 
                         growth_model = gm, growth_rate = gr, norm = "norm",
                         nvar = nvar, burnin = burn, noff = noff,
                         mrr = reprorate, frr = reprorate, msr = surrate, fsr = surrate,
                         # Add translocation parameters
                         trans = tr_in, tsim = F, tsr = tsr, tyears = ytran, ntrans = ntran,
                         trans_year = yeffect, trr = trr,
                         # Add bad year parameters
                         cyears = bad_years, cat_prop = bad_prop)

save(vtalTransBad, file = "4.translocation/vtalTransBad.Rdata")

# Calculate population descriptions
summarize_n(vtTransBad)

# Calculate expected heterozygosity
exp <- gd_ret_exp(out = vtTransBad, ne = "Nc")
exp <- apply(exp, 1, mean)

# Calculate predicted heterozygosity
pred <- gd_ret_pred(out = vtTransBad, gd = "He")
het <- cbind(pred[,1], exp)

# Create summary plot



