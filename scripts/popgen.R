#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' L. perennis Popgen Calculations
#' @date 2024-07-31
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'

library(ggplot2)
library(tidyverse)
library(hierfstat)
library(adegenet)
library(performance)
library(broom)
library(vcfR)

# Load VCF
nevcf <- read.vcfR("2.stacks/neFiltered.vcf")

neList <- read.delim(file = "2.stacks/nemap.txt", header = FALSE) |>
  rename(id = V1,
         pop = V2) |>
  filter(id %in% colnames(nevcf@gt))

# Convert to genind
neind <- vcfR2genind(nevcf, pop = neList$pop, return.alleles = TRUE)

# Calculate locus-level diversity 
lpHier <- genind2hierfstat(neind)

#lp.div <- basic.stats(data = lpHier, diploid = TRUE, digits = 4)

# Calculate mean gene diversities by population
div <- data.frame(Hs = hierfstat::Hs(lpHier),
                  Ho = hierfstat::Ho(lpHier),
                  region = c("Vermont", "New Hampshire", "New Hampshire",
                             "New York", "Vermont",
                             "New Hampshire", "New Hampshire",
                             "Massachusetts", "New York"))

# Calculate mean observed heterozygosities by population
lpHet <- hierfstat::Ho(lpHier)

#write.csv(x = div, file = "3.popgen/neDiversity.csv")

div <- read.csv(file = "3.popgen/neDiversity.csv") |>
  rename(Pop = X)

palReg <- c( "#E2A3C7", "#778da9", "#EC7D10", "#63A46C")


## Model ecosystem area as predictor of diversity
hetMod <- div |>
  filter(!is.na(acres)) |>
  mutate(kacres = acres/1000) |>
  lm(formula = Ho ~ kacres)

check_posterior_predictions(hetMod)
check_normality(hetMod)
check_heteroscedasticity(hetMod)

tidy(hetMod) 
# 0.00253 increase in observed heterozysogisty per 1000 acres p = 0.357


#### Graph diversity
lpDivNe <- ggplot(data = div,
                   mapping = aes(x = Hs, y = Ho, 
                                 color = region)) +
  scale_color_manual(values = palReg) +
  geom_point(size = 4) +
  labs(title = "L. perennis Regional \nHeterozygosity") +
  theme_classic(base_size = 16)+
  geom_segment(aes(x = 0.018, xend = 0.15,
                   y = 0.018, yend = 0.15),
               linewidth = 0.1, linetype = "dashed")

#png(filename= "3.popgen/lpDivNe.png", width = 900, height = 500)
#lpDivNe
#dev.off()

#### Population Differentiation
# Fst between populations
lpDist <- genet.dist(lpHier)

lpDist |>
  as.matrix() |>
  write.csv("3.popgen/Fst.csv")

# IBD Analysis
library(ade4)
library(vegan)

# Load in pairwise genetic and geographic distance matrices
fstDist <- read.csv("3.popgen/Fst.csv") |>
  select(!X) |>
  unname()

geoDist <- read.csv("4.IBD/geoDist.csv") |>
  select(!X) |>
  unname()

# Perform mantel test with vegan on the two matrices
vegan::mantel(as.dist(fstDist), as.dist(geoDist), method = "pearson", na.rm = TRUE)

# r=0.79 with p = 0.001


#### Inbreeding
# Calculate population-level Fis CIs by bootstrapping
inbredHier <- boot.ppfis(lpHier, nboot = 100)

Region = c("Vermont", "New Hampshire", "New Hampshire", "New York",
           "Vermont", "New Hampshire",
           "New Hampshire", "Massachusetts",
           "New York")

inbred <- data.frame(ll = inbredHier$fis.ci[,,1], 
           pop = unique(neList$pop),
           Region = Region) |>
  rename(ll = ll.ll,
         hl = ll.hl) |>
  mutate(Fis = (ll+hl)/2)

#write.csv(inbred, "3.popgen/inbreeding.csv")

inbred <- read_csv("3.popgen/inbreeding.csv") |>
  select(!...1)

palReg <- c( "#E2A3C7", "#778da9", "#EC7D10", "#63A46C")

# Plot Fis with confidence intervals
FisPlot <- ggplot(data = inbred,
                  mapping = aes(x = pop, y = Fis, color = Region)) +
  geom_point(size = 2) +
  geom_errorbar(mapping = aes(ymin = ll, ymax = hl)) +
  theme_classic(base_size = 16) +
  labs(title = "Fis confidence intervals with 1000 bootstraps",
       x = "Population")

#png(filename= "FisPlot.png", width = 700, height = 500)
#FisPlot
#dev.off()

# Just Northeast region Fis with confidence intervals
inbredRegions <- inbredHier |>
  filter(Region %in% c("Massachusetts", "New Hampshire", "New York", "Vermont")) |>
  ggplot(mapping = aes(x = pop, y = Fis, color = Region)) +
  geom_point(size = 2) +
  geom_errorbar(mapping = aes(ymin = ll, ymax = hl)) +
  theme_classic(base_size = 16) +
  labs(title = "Fis confidence intervals with 1000 bootstraps",
       x = "Population")

#png(filename= "FisRegions.png", width = 700, height = 500)
#inbredRegions
#dev.off()

# Model population area against Fis
fisMod <- inbredHier |>
  filter(!is.na(acres)) |>
  mutate(kacres = acres/1000) |>
  lm(formula = Fis ~ kacres)

check_posterior_predictions(fisMod)
check_normality(fisMod)
check_heteroscedasticity(fisMod)

tidy(fisMod) 
# 0.0280 increase in Fis per 1000 acres p = 0.416
