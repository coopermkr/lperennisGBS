#'''''''''''''''''''''''''''''''''''''
#' Vortex Parameters
#' @author Cooper Kimball-Rhines
#' @date 2026-01-13
#'''''''''''''''''''''''''''''''''''''


# Load libraries
library(tidyverse)
library(performance)
library(broom)
library(emmeans)
library(glmmTMB)
library(DHARMa)
library(grid)
library(gridExtra)

meta <- tibble(Population = c("AL", "MO", "SA", "AB", "ABAL", "ABMO", "ABSA", "CN", "CNAL", "CNMO", "CNSA", "HK", "HKAL", "HKMO", "HKSA"),
               nPops = str_length(Population)/2,
               size = c("large", "small", "medium", "small", "small-large", "small-small", "small-medium",
                        "large", "large-large", "small-large", "medium-large",
                        "medium", "medium-large", "small-medium", "medium-medium"),
               risk = as.factor(c("N", "Y", "N", "Y", "N-Y", "Y-Y", "N-Y", "N", "N-N", "N-Y", "N-N", "N", "N-N", "N-Y", "N-N")),
               Fst = c(NA, NA, NA, NA, 0.021, 0.021, 0.021, NA, 0.017, 0.018, 0.018, NA, 0.02, 0.021, 0.021))

#### Load in translocation vortex runs ####
tNe <- read_csv("6.vortexAnalysis/transNeCollated.csv",
                         col_names = TRUE) |>
  na.omit()

ggplot(tNe, mapping = aes(x = Year, y = Ne, color = Population)) +
  geom_point() + theme_light()

# Filter out points greater than census size (Ne miscalculation, most likely infinite values)
filtNe <- tNe |> filter(Ne < 72)

ggplot(filtNe, mapping = aes(x = Year, y = Ne, color = Population)) +
  geom_point() + theme_light()
# There doesn't appear to be much of a pattern, but maybe when we combine them as groups

ggplot(filtNe, mapping = aes(x = Ne)) +
  geom_histogram() + theme_light()

#### Is dual better than mono? ####


filtNeMeta <-  merge(meta, filtNe, by = "Population")

# Is one population better than two?
monoDual <- lmer(data = filtNeMeta, formula = Ne ~ nPops + Year + (1|Population))

check_model(monoDual) # The fit isn't very good
tidy(monoDual) # Number of populations is not significant

# Retry as a gamma model
library(glmmTMB)
library(DHARMa)

monoDualGLM <- glm(data = filtNeMeta, formula = Ne ~ Population + Year,
                   family = Gamma(link = "log"))

check_model(monoDualGLM) # Fit is way better

monoDualRes <- simulateResiduals(monoDualGLM, n = 500)

plotQQunif(monoDualRes,
           testUniformity = FALSE,
           testOutliers = FALSE,
           testDispersion = FALSE)
#Looks pretty good

#Inference
summary(monoDualGLM)
tidy(monoDualGLM)
r2(monoDualGLM) # Terrible R2

neEM <- emmeans(monoDualGLM, specs = ~ nPops) |>
  contrast(method = "pairwise") |>
  confint() # No significance

# What about average Ne over the entire simulation?
neSum <- filtNeMeta |>
  group_by(Replicate, Population, nPops) |>
  summarize(meanNe = mean(Ne),
            sdNe = sd(Ne))

# Check histogram-- is a gamma still necessary?
ggplot(neSum, mapping = aes(x = meanNe)) +
  geom_histogram() + theme_light() # Looks pretty normally distributed

meanModel <- lm(data = neSum, formula = meanNe ~ nPops)

check_model(meanModel) # Some extremely influential points but should be ok
tidy(meanModel) # No significance

meanEM <- emmeans(meanModel, specs = ~ nPops) |>
  contrast(method = "pairwise") |>
  confint() # No significance

ggplot(data = filtNeMeta, mapping = aes(x = nPops, y = Ne, group = nPops)) +
  geom_boxplot() + theme_light()

#### Does donor pop size matter? ####
dualFounders <- filtNeMeta |> filter(nPops > 1)

dualFounders |> group_by(size) |> count() # LL, MM, and SS have fewer replicates than the rest

# Check distribution-- probably going to need a gamma again
ggplot(dualFounders, mapping = aes(x = Ne)) +
  geom_histogram() # Yes need gamma again

dualGLM <- glm(data = dualFounders, formula = Ne ~ size + Year,
                   family = Gamma(link = "log"))

# check model fit
check_model(dualGLM) # Looks good

dualRes <- simulateResiduals(dualGLM, n = 500)

plotQQunif(dualRes,
           testUniformity = FALSE,
           testOutliers = FALSE,
           testDispersion = FALSE) # Looks good

# Analyze model coefficients
tidy(dualGLM)

dualEM <- emmeans(dualGLM, specs = ~ size) |>
  contrast(method = "pairwise") |>
  confint() # No significance

dualEM |> plot() + geom_vline(xintercept = 0) + theme_light() # No significance!

#### Does Fst affect Ne? ####
fstGLM <- glm(data = dualFounders, formula = Ne ~ Fst + Year,
               family = Gamma(link = "log"))

check_model(fstGLM)
tidy(fstGLM) # Nope!- this is continuous so no EMmeans


## How quickly does the boost of having a second population decay?
# Bin years by 10s
yearBin <- function(cutoff) {
  
  filtNeMeta |> 
    filter(Year > cutoff,
           Year < cutoff+10) |>
    group_by(nPops) |>
    summarize(aveNe = mean(Ne),
              sdNe = sd(Ne),
              bin = cutoff)
}

bins <- as.vector(seq(from = 0, to = 150, by = 10))
binnedNe <- map(.x = bins, .f = yearBin) |> list_rbind()

ggplot(binnedNe, mapping = aes(x = bin, y = aveNe, group = nPops)) +
  geom_point()


#### With program output ####
suppSumm <- read_csv("6.vortexAnalysis/suppSummary.csv")

# join with metadata
suppMeta <- merge(meta, suppSumm, by = "Population") # Should also filter out metapopulations

# Is there a difference between one and two-pop founders?
npopGLM <- glm(data = suppMeta, formula = meanNe ~ nPops,
                   family = Gamma(link = "log"))

check_model(npopGLM) # 

npopRes <- simulateResiduals(npopGLM, n = 500)

plotQQunif(npopRes,
           testUniformity = FALSE,
           testOutliers = FALSE,
           testDispersion = FALSE)
#Looks pretty good

#Inference
summary(npopGLM)
tidy(npopGLM)

npopEM <- emmeans(npopGLM, specs = ~ nPops) |>
  contrast(method = "pairwise") |>
  confint() # No significance

# Plot population vs. meanNe
popNePlot <- suppMeta |>
  filter(nPops == 1) |>
  mutate(Population = as.factor(Population),
         Population = fct_recode(Population, "NH3" = "AB","NY1" = "AL","NH1" = "CN","NH2-22" = "HK", 
                           "MA1" = "MO","NY2" = "SA"),
         Population = factor(Population, levels = c("MA1", "NH1", "NH2-22", "NH3", "NY1", "NY2"))) |>
  ggplot(suppMeta, mapping = aes(x = Population, y = meanNe)) +
  geom_boxplot() + theme_light(base_size = 16) +
  labs(x = "Seed Donor Population", y = "Mean Effective Population Size (Ne)",
       subtitle = "Single-Source Founded Population Simulations", tag = "(A)")

# Plot population size vs. meanNe
meanNePlot <- suppMeta |>
  mutate(size = as.factor(size),
         size = factor(size, levels = c("large", "medium", "small", "large-large", "medium-large", "small-large",
                       "medium-medium", "small-medium", "small-small")),
         size = fct_recode(size, "Large" = "large","Medium" = "medium","Small" = "small","Large-Large" = "large-large", 
                           "Large-Medium" = "medium-large","Large-Small" = "small-large",
                          "Medium-Medium" = "medium-medium","Medium-Small" = "small-medium","Small-Small" = "small-small")) |>
ggplot(suppMeta, mapping = aes(x = size, y = meanNe)) +
  geom_boxplot() + theme_light(base_size = 16) +
  labs(x = "Seed Donor Population Size", y = "Mean Effective Population Size (Ne)",
       subtitle = "Source Size Effect on Founded Population Simulations", tag = "(B)") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))
  
jpeg(filename= "6.vortexAnalysis/MeanNe.jpg", width = 8, height = 10, units = "in", res = 100)
grid.arrange(popNePlot, meanNePlot, ncol = 1)
dev.off()


# Analyze differences between mono-founder populations first
suppMono <- suppMeta |> filter(nPops == 1)

ggplot(suppMono, mapping = aes(x = size, y = meanNe)) +
  geom_boxplot() + theme_light()

ggplot(suppMono, mapping = aes(x = meanNe)) +
  geom_histogram() + theme_light()

# Fit GLM for meanNe
suppMonoGLM <- glm(data = suppMono, formula = meanNe ~ size,
                   family = Gamma(link = "log"))

check_model(suppMonoGLM) # Fit is way better

monoDualRes <- simulateResiduals(suppMonoGLM, n = 500)

plotQQunif(suppMonoGLM,
           testUniformity = FALSE,
           testOutliers = FALSE,
           testDispersion = FALSE)
#Looks pretty good

#Inference
summary(suppMonoGLM)
tidy(suppMonoGLM)

monoEM <- emmeans(suppMonoGLM, specs = ~ size) |>
  contrast(method = "pairwise") |>
  confint() # No significance

monoEM |> plot() + geom_vline(xintercept = 0)

kruskal.test(data = suppMono, meanNe ~ size)


### Differences between dual-founder populations
suppDuo <- suppMeta |> filter(nPops > 1)

ggplot(suppDuo, mapping = aes(x = size, y = meanNe)) +
  geom_boxplot() + theme_light() + geom_point()

ggplot(suppDuo, mapping = aes(x = meanNe)) +
  geom_histogram() + theme_light()

# Fit GLM for meanNe
suppDuoGLM <- glm(data = suppDuo, formula = meanNe ~ size,
                   family = Gamma(link = "log"))

check_model(suppDuoGLM) # Good fit, though there are some weirdly influential observations

suppDuoRes <- simulateResiduals(suppDuoGLM, n = 500)

plotQQunif(suppDuoRes,
           testUniformity = FALSE,
           testOutliers = FALSE,
           testDispersion = FALSE) # Seems ok, not that many data points
#Looks pretty good

#Inference
summary(suppDuoGLM)
tidy(suppDuoGLM)

suppDuoEM <- emmeans(suppDuoGLM, specs = ~ size) |>
  contrast(method = "pairwise", adjust = "tukey") # No significance

suppDuoEM |> plot() + geom_vline(xintercept = 0)

anova(suppDuoGLM)
kruskal.test(data = suppDuo, meanNe ~ size)

# Model at-risk status effect on Ne
meanNePlot <- suppMeta |>
  ggplot(suppMeta, mapping = aes(x = risk, y = meanNe)) +
  geom_boxplot() + theme_light(base_size = 16) +
  labs(x = "Seed Donor At-Risk Status", y = "Mean Effective Population Size (Ne)") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

# Construct and check model
riskGLM <- glm(data = suppMeta, formula = meanNe ~ risk,
                   family = Gamma(link = "log"))

check_model(riskGLM) # Correct shape
r2(riskGLM) # Fit is bad

riskRes <- simulateResiduals(riskGLM, n = 500)

plotQQunif(riskRes,
           testUniformity = FALSE,
           testOutliers = FALSE,
           testDispersion = FALSE)
#Looks pretty good

#Inference
summary(riskGLM)

riskEM <- emmeans(riskGLM, specs = ~ risk) |>
  contrast(method = "pairwise") |>
  confint() # No significance

riskEM |> plot() + geom_vline(xintercept = 0)


