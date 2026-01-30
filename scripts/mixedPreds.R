#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' L. perennis Census Mixed Model NHFG
#' @date 2026-01-20
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Load libraries
library(tidyverse)
library(lme4)
library(broom.mixed)
library(performance)
library(modelbased)

### Transects have minimal width, only intercepting clumps are counted.

#### Load data ####
census <- read_csv("data/Combined_Lupine_Sheets.csv") |>
  # Calculate log scaled sizes
  mutate(Ellipse = Length*Width*pi/4) |>
  filter(!is.na(Ellipse)) %>%
  replace(is.na(.), 0)

# Chart size distribution by year
ggplot(census, mapping = aes(x = Ellipse)) +
  geom_histogram() +
  facet_wrap(vars(Year)) # For validation: 2010/2011, 2017/2018, and 2023/2024 need to be combined

# Number of flowers distribution by year
ggplot(census, mapping = aes(x = Flowers)) +
  geom_histogram() +
  facet_wrap(vars(Year))

census |>
  filter(Flowers > 0) |>
  ggplot(mapping = aes(x = Flowers)) +
  geom_histogram() +
  facet_wrap(vars(Year))

# Calculate probability of flowering and average number of flowers/flowering individual
census |>
  filter(Flowers > 0) |>
  summarize(nflow = mean(Flowers))

sum(census$Flowers == 0)/nrow(census)

#### Calculate BaseLine observations ####
cenSum <- census |>
  group_by(BaseLine, Year) |>
  summarize(cov = sum(Ellipse),
            count = n(),
            nflowers = sum(Flowers))

# Metadata from 2001 lupine transect establishment datasheet
meta <- read_csv(file = "data/censusMeta.csv") |>
  mutate(AREA_SURVEYED = NTRANSECTS*TRANSECT_LENGTH) |>
  rename(BaseLine = SITE)

# Merge coverage data with metadata
basePop <- merge(cenSum, meta) |>
  mutate(countSurveyed = count/AREA_SURVEYED,
         covSurveyed = cov/AREA_SURVEYED,
         totalArea = BASELINE_LENGTH*TRANSECT_LENGTH,
         countTotal = count/totalArea,
         covTotal = cov/totalArea) |>
  select(BaseLine, Year, countSurveyed, covSurveyed, countTotal, covTotal, totalArea,
         cov, count, nflowers, AREA_SURVEYED)


#### Model population density as a function of continuous Year and RE for baseline ####
# We know lupine is affected by density dependent reproduction rates so let's start here
censusMixed <- lmer(data = basePop,
                    formula = covSurveyed ~ Year + (1|BaseLine))

check_model(censusMixed) # Model fits well

# Assess coefficients
tidy(censusMixed) # Steady population density decline of ~ 1.63cm^2/m of baseline per year
summary(censusMixed)
confint(censusMixed) # Population density decline is not significantly negative
estimate_relation(censusMixed) |> plot() # Visualize year fixed effect

#Visualize random effects
estimate_relation(censusMixed,
                  include_random = TRUE) |> plot()

#Show the marginal (fixed effect) curve on top of the random effects
estimate_relation(censusMixed,
                  include_random = TRUE) |>
  plot(ribbon = list(alpha = 0)) +
  geom_line(data = estimate_relation(censusMixed),
            mapping = aes(x = Year,
                          y = Predicted),
            color = "black",
            linewidth = 2) +
  theme_classic()

#### Model population count density as a function of continuous Year and RE for baseline ####
# This should theoretically give the same result as a glm with dens ~ Year + AREA_SURVEYED
# And probably a similar result to count ~ Year + (1|BaseLine)
countsDensMixed <- lmer(data = basePop,
                     formula = countSurveyed ~ Year + (1|BaseLine))

check_model(countsDensMixed) # Model fits well

# Assess coefficients
tidy(countsDensMixed) # Population decline of ~ 0.004 individuals/m of baseline per year
summary(countsDensMixed)
confint(countsDensMixed) # Population decline is significantly negative
estimate_relation(countsDensMixed) |> plot() # Visualize year fixed effect

#Visualize random effects
estimate_relation(countsDensMixed,
                  include_random = TRUE) |> plot()

#Show the marginal (fixed effect) curve on top of the random effects
estimate_relation(countsDensMixed,
                  include_random = TRUE) |>
  plot(ribbon = list(alpha = 0)) +
  geom_line(data = estimate_relation(countsDensMixed),
            mapping = aes(x = Year,
                          y = Predicted),
            color = "black",
            linewidth = 2) +
  theme_classic()

#### Overall Population Size Trend ####
# Estimate total population size each year-- this involves combining subpops so no error bars
# Though using the previous numbers you could calculate a percent variation and scale it up to this level
censusMetrics <- basePop |>
  group_by(Year) |>
  summarize(countDensity = average(countSurveyed),
            countSD = sd(countSurveyed)
            censusFlowers = sum(nflowers),
            censusArea = sum(totalArea),
            censusCoverage = sum(cov)) |> # Summarize transects to the area of subpopulation
  mutate(countDensity = censusCount/censusArea, # Calculate overall 
         popSize = countDensity*)

popSizeMixed <- lmer(data = basePop,
                     formula = dens ~ Year + (1|BaseLine))

check_model(countsDensMixed) # Model fits well

# Assess coefficients
tidy(countsDensMixed) # Population decline of ~ 0.004 individuals/m of baseline per year
summary(countsDensMixed)
confint(countsDensMixed) # Population decline is significantly negative
estimate_relation(countsDensMixed) |> plot() # Visualize year fixed effect

#Visualize random effects
estimate_relation(countsDensMixed,
                  include_random = TRUE) |> plot()

#Show the marginal (fixed effect) curve on top of the random effects
estimate_relation(countsDensMixed,
                  include_random = TRUE) |>
  plot(ribbon = list(alpha = 0)) +
  geom_line(data = estimate_relation(countsDensMixed),
            mapping = aes(x = Year,
                          y = Predicted),
            color = "black",
            linewidth = 2) +
  theme_classic()



