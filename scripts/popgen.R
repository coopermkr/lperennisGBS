#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' L. perennis Popgen Calculations
#' @date 2024-07-31
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'

library(ggplot2)
library(ggpubr)
library(tidyr)
library(hierfstat)
library(adegenet)
library(performance)
library(broom)

#### Calculate locus-level diversity 
lpHier <- genind2hierfstat(lpInd)

#lp.div <- basic.stats(data = lpHier, diploid = TRUE, digits = 4)

# Calculate mean gene diversities by population
div <- data.frame(Hs = hierfstat::Hs(lpHier),
                  Ho = hierfstat::Ho(lpHier),
                  region = c("Vermont", "New Hampshire", "New Hampshire",
                             "New York", "Florida", "Florida", "Mid", "Vermont",
                             "New Hampshire", "Mid", "New Hampshire", "Mid", 
                             "Massachusetts", "Florida", "Mid", "Mid", "New York"))

# Calculate mean observed heterozygosities by population
lpHet <- hierfstat::Ho(lpHier)

write.csv(x = div, file = "lpDiversity.csv")

div <- read.csv(file = "lpDiversity.csv") |>
  mutate(Pop = X,
         Region = c("Vermont", "New Hampshire", "New Hampshire", "New York",
                    "Florida", "Florida", "Mid", "Vermont", "New Hampshire",
                    "Mid", "New Hampshire", "Mid", "Massachusetts",
                    "Florida", "Mid", "Mid", "New York"))

palReg <- c("#C200FB","#E2A3C7", "#DC136C", "#778da9", "#EC7D10", "#63A46C")


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
lpDivAll <- ggplot(data = div,
                   mapping = aes(x = Hs, y = Ho, 
                                 color = Region)) +
  scale_color_manual(values = palReg) +
  geom_point(size = 4) +
  labs(title = "L. perennis Regional \nHeterozygosity") +
  theme_classic(base_size = 16)+
  geom_segment(aes(x = 0.018, xend = 0.15,
                   y = 0.018, yend = 0.15),
               linewidth = 0.1, linetype = "dashed")

png(filename= "lpHetAll.png", width = 900, height = 500)
lpDivAll
dev.off()


palReg <- c("#E2A3C7", "#DC136C", "#778da9", "#EC7D10", "#63A46C")

lpDivReg <- ggplot(data = div |> filter(Hs < 0.1),
                   mapping = aes(x = Hs, y = Ho, 
                                 color = Region)) +
  scale_color_manual(values = palReg) +
  geom_point(size = 4) +
  labs(title = "L. perennis Population Heterozygosity") +
  theme_classic(base_size = 16)+
  geom_segment(aes(x = 0.018, xend = 0.045,
                   y = 0.018, yend = 0.045),
               linewidth = 0.1, linetype = "dashed")

png(filename= "lpDivReg.png", width = 900, height = 500)
lpDivReg
dev.off()

#### Population Differentiation
# Fst between populations
lpDist <- genet.dist(lpHier)

lpDist |>
  write.table("Fst.csv")

# IBD Analysis
library(ade4)
library(vegan)

# Load in pairwise genetic and geographic distance matrices
fstDist <- read.csv("Fst.csv") |>
  select(!X) |>
  unname()

geoDist <- read.csv("geoDist.csv") |>
  select(!X) |>
  unname()

# Perform mantel test with vegan on the two matrices
vegan::mantel(as.dist(fstDist), as.dist(geoDist), method = "pearson", na.rm = TRUE)

# r=0.75 with p = 0.001


#### Inbreeding
# Calculate population-level Fis CIs by bootstrapping
inbredHier <- boot.ppfis(lpHier, nboot = 100) |>
  cbind(inbred$fis.ci, unique(popList$pop)) |>
  rename(pop = `unique(popList$pop)`) |>
  mutate(Fis = (ll+hl)/2)

write.csv(inbredHier, "inbreeding.csv")

inbredHier <- read.csv("inbreeding.csv") |>
  mutate(Region = c("Vermont", "New Hampshire", "New Hampshire", "New York",
                    "Florida", "Florida", "Mid", "Vermont", "New Hampshire",
                    "Mid", "New Hampshire", "Mid", "Massachusetts",
                    "Florida", "Mid", "Mid", "New York"))

palReg <- c("#C200FB","#E2A3C7", "#DC136C", "#778da9", "#EC7D10", "#63A46C")

# Plot Fis with confidence intervals
FisPlot <- ggplot(data = inbredHier,
                  mapping = aes(x = pop, y = Fis, color = Region)) +
  geom_point(size = 2) +
  geom_errorbar(mapping = aes(ymin = ll, ymax = hl)) +
  theme_classic(base_size = 16) +
  labs(title = "Fis confidence intervals with 1000 bootstraps",
       x = "Population")

png(filename= "FisPlot.png", width = 700, height = 500)
FisPlot
dev.off()

# Just Northeast region Fis with confidence intervals
inbredRegions <- inbredHier |>
  filter(Region %in% c("Massachusetts", "New Hampshire", "New York", "Vermont")) |>
  ggplot(mapping = aes(x = pop, y = Fis, color = Region)) +
  geom_point(size = 2) +
  geom_errorbar(mapping = aes(ymin = ll, ymax = hl)) +
  theme_classic(base_size = 16) +
  labs(title = "Fis confidence intervals with 1000 bootstraps",
       x = "Population")

png(filename= "FisRegions.png", width = 700, height = 500)
inbredRegions
dev.off()

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
