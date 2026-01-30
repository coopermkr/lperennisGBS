#'''''''''''''''''''''''''''''''''''''
#' Vortex Parameters
#' @author Cooper Kimball-Rhines
#' @date 2026-01-13
#'''''''''''''''''''''''''''''''''''''

# Load libraries
library(tidyverse)
library(glmmTMB)
library(performance)
library(broom)
library(modelr)

#### Minimum Viable Population Size ####
# Sensitivity tests were run for 300 iterations for 140 years (~40 generations)
# With factorial parameters for catastrophes at 0.5% increments and pop sizes/K up from 1-250
# Load in replicate sensitivity tests
mvp1 <- read_csv("5.translocation/VOutput/MVP1.csv")
mvp2 <- read_csv("5.translocation/VOutput/MVP2.csv")
mvp3 <- read_csv("5.translocation/VOutput/MVP3.csv")
mvp4 <- read_csv("5.translocation/VOutput/MVP4.csv")
mvp5 <- read_csv("5.translocation/VOutput/MVP5.csv")

mvp <- rbind(mvp1, mvp2, mvp3, mvp4, mvp5) |>
  rename(pcat = SV1,
         N0 = SV2) |>
  mutate(nExtinct = PE*nRuns)

# Map the parameter space- factorial
ggplot(mvp, mapping = aes(x = N0, y = pcat)) +
  geom_point()

# Plot starting population size vs. extinction probability
ggplot(mvp, mapping = aes(x = N0, y = PE)) +
  geom_point() + facet_wrap(vars(pcat)) + theme_light()

ggplot(mvp, mapping = aes(x = N0, y = PE, weight = nRuns, color = pcat)) +
  geom_point() + geom_hline(yintercept = 0.01, linetype = "dashed") +
  theme_light()

# Plot catastrophe probability vs. extinction probability
ggplot(mvp, mapping = aes(x = pcat, y = PE)) +
  geom_point()

mvp |> select(nRuns, nExtinct, PE, N0, pcat) |> filter(PE > 0.01) |> arrange(desc(N0))

# Extinction probability at a catastrophe level predicted by Frankham:
mvp |> filter(pcat == 5) |>
  ggplot(mapping = aes(x = N0, y = PE)) +
  geom_point() + geom_hline(yintercept = 0.01, linetype = "dashed") +
  theme_light()

ggplot(mvp, mapping = aes(x = PE)) +
  geom_histogram() + theme_light()

#### Calculating MVP ####
# Fit a model to the replicates predicting extinction probability
# Quasibinomial should be correct to deal with the PE = 1 and 0

mvpFilt <- filter(mvp, N0 < 150)

mvpGlm <- glm(formula = PE ~ N0 + pcat,
              family = binomial(link = "logit"),
              data = mvp,
              weights = nRuns)

check_model(mvpGlm) # This model is bad and shouldn't be trusted
summary(mvpGlm)

#Inference: coefficients, fit, and dispersion parameter
tidy(mvpGlm)
glance(mvpGlm)
summary(mvpGlm)

# Plot line over data points
mvp |>
  ggplot(mapping = aes(x = N0, y = PE, weight = nRuns)) +
  geom_point() + geom_hline(yintercept = 0.01, linetype = "dashed") +
  theme_light() +
  stat_smooth(method = "glm", method.args= list(family = binomial))

# Pull out MVP directly from simulations
minCat <- mvp |>
  filter(PE >= 0.01) |>
  group_by(pcat, Population) %>%
  summarize(mvpSize = max(N0))

mvp |> filter(PE >= 0.01, pcat == 5) |> arrange(desc(N0)) |> select(N0, PE, `N-extant`)
# For Frankham's catastrophe = 5%: MVP=107
# For severe drought catastrophe = 2%: MVP=65

catPlot <- ggplot(minCat, mapping = aes(x = pcat, y = mvpSize)) +
  geom_point() + theme_light()

# Fit pcat as predictor of MVP
pcatMod <- lm(mvpSize ~ pcat,
              data = minCat)

check_model(pcatMod) # Good

# Assess model
r2(pcatMod)
tidy(pcatMod) # Significant coefficient: each increment of pcat increases MVP by 7

# Plot the line, calculate means and standard deviations
mvpPlot <- catPlot + stat_smooth(method = "lm") +
  labs(y = "Simulated Minimum Viable Population",
       x = "Probability of Catastrophe") +
  guides() + theme_light(base_size = 16)

jpeg(filename= "5.translocation/mvp.jpg", width = 800, height = 800, quality = 200)
mvpPlot
dev.off()

minCat |>
  group_by(pcat) |>
  summarize(meanMVP = mean(mvpSize),
            sdMVP = sd(mvpSize)) |>
  print(n = 25)

# More catastrophe means higher MVP, but uncertainty is cyclical (weirdly)
# Report the model coefficient in text and use the mean MVP for simulations
