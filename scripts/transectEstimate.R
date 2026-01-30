#'''''''''''''''''''''''''''''''''''''''''''''''
#' Transect Point-Intercept Population Estimate
#' @author Cooper Kimball-Rhines
#' @date 2026-01-27
#'''''''''''''''''''''''''''''''''''''''''''''''

#' This script interprets point-intercept transect data collected by New Hampshire Fish and Game
#' from 2001-2023 into population size estimates for the surveyed subpopulations of L. perennis.
#' Raw data provided by Heidi Holman, NHFG.

# Load libraries
library(tidyverse)
library(performance)
library(broom)

# Load collated data sheet
transects <- read_csv("data/combinedTransectIntercepts.csv") |>
  # Combine 2010 and 2011 which together make a complete set of transects
  mutate(YEAR = replace(YEAR, YEAR %in% c(2010, 2011), 2010.5),
         YEAR = replace(YEAR, YEAR %in% c(2017, 2018), 2017.5))

meta <- read_csv("data/transectMeta.csv") |> 
  mutate(totalLength = TRANSECT_LENGTH * 100 * N_TRANSECTS,
         siteArea = BASELINE_LENGTH * 100 * TRANSECT_LENGTH * 100) |>
  select(SITE, N_TRANSECTS, TRANSECT_LENGTH, BASELINE_LENGTH, totalLength, siteArea)

# Which sites were surveyed in which years?
transects |> group_by(YEAR) |> reframe(sites = unique(SITE)) |>
  ggplot(mapping = aes(x = as.factor(YEAR), y = sites)) +
  geom_tile()

# The raw data includes more species than just L. perennis, which should be excluded
unique(transects$SPP)

transectFilt <- transects |>
  filter(SPP %in% c("LUPE", "LP", "lp", "L")) |>
  mutate(TRANSECT = str_split_i(INTERCEPT, pattern = "-", 1), # Pull transect from intercept column
         INTcm = as.numeric(INTcm)) |> # Fix data type of intercept
  filter(!is.na(INTcm),
         INTcm > 0) |> # Filter out missing intercept data
  select(SITE, YEAR, TRANSECT, INTcm) # Only keeping columns relevant to this calculation

#' Eberhardt 1978 presents equations for population estimate from point intercept data
#' https://www.jstor.org/stable/pdf/3800685.pdf
#' Chords are not recorded, so the equation assumes all plants are perfect circles
#' There are two ways to do this calculation: Consider transects as independent observations
#' and use them to calculate independent population estimates so you can get SD and SE
#' or sum the transect values and calculate them all at once.

#### Option 1: Transects are independent observations
intLength1 <- transectFilt |>
  group_by(YEAR, SITE, TRANSECT) |>
  summarize(invINT = sum(1/INTcm), # Sum 1/interceptLength of each plant 
            plants = n())

# Merge with metadata to get transect lengths
intMeta1 <- meta |>
  merge(intLength1, by = "SITE")

# Make a per-transect density and subpopulation size estimate
estimates1 <- intMeta1 |>
  mutate(D1 = ((2/pi)*invINT)/(TRANSECT_LENGTH * 100),
         N1 = D1 * siteArea)

# Some transects do not show up in the dataset because they have no lupine
missingData <- intMeta1 |> group_by(SITE, YEAR) |>
  summarize(obsTransects = n(),
            missTransects = N_TRANSECTS - n()) |>
  distinct()

aveEst1 <- merge(missingData, estimates1) |>
  group_by(SITE, YEAR) |>
  summarize(aveD1 = sum(D1)/N_TRANSECTS,
            aveN1 = sum(N1)/N_TRANSECTS) |>
  distinct()


#### Option 2: Sum all transect intercepts by site
intLength2 <- transectFilt |>
  group_by(YEAR, SITE) |>
  summarize(invINT = sum(1/INTcm), # Sum 1/interceptLength of each plant 
            plants = n())

# Merge with metadata to get transect lengths
intMeta2 <- meta |>
  merge(intLength2, by = "SITE")

# Make a per-transect density and subpopulation size estimate
aveEst2 <- intMeta2 |>
  mutate(aveD2 = ((2/pi)*invINT)/totalLength,
         aveN2 = aveD2 * siteArea) |>
  select(SITE, YEAR, aveD2, aveN2)

# Both options result in the same answer, because the math is the same. 
# Option 2 is quicker and easier to understand, but option 1 could be used to 
# calculate standard deviations if zeros are accounted for


#### Change in Population Size
#' The problem with this is that not all sites were surveyed every year. If sites
#' were skipped because they were developed/because plants do not grow there anymore, 
#' then the years are comparable. But if they were skipped by mistake, then they are not.
#' 
yearEst <- aveEst1 |> group_by(YEAR) |>
  summarize(N = sum(aveN1))

### Assuming skipped sites had zero plants and 2014 was a fluke:
yearEst |>
  filter(YEAR != 2014) |>
  ggplot(mapping = aes(x = YEAR, y = N)) +
  geom_point() +
  stat_smooth(method = "lm") +
  theme_light(base_size = 16) +
  labs(x = "Year", y = "Transect Population Size")

# Create the model seen in the graph
estLM <- yearEst |>
  filter(YEAR != 2014) %>%
  lm(formula = N ~ YEAR, data = .)

check_model(estLM) # These look weird because there's very little data
tidy(estLM)
# We're seeing an average decrease of 2243 individuals in the transect-measured
# population per year (SE = 326, p = 0.006)

### Assuming change needs to be calculated per site:
library(lme4)
library(broom.mixed)
library(modelbased)

aveEst2 |>
  ggplot(mapping = aes(x = YEAR, y = aveN2, color = SITE)) +
  geom_point() +
  stat_smooth(method = "lm")

estRE <- lmer(formula = aveN2 ~ YEAR + (1|SITE),
              data = aveEst2)

check_model(estRE) # Fits well (but there are also very few data points, so it should)

tidy(estRE) # Each subpopulation is losing an average of 395 individuals/year
# So for 8 subpopulations this equals a loss of ~3160 plants/year
summary(estRE)
confint(estRE) # Average subpopulation decline is negative with high confidence (-585, -204)
estimate_relation(estRE) |> plot() # Visualize year fixed effect

#Visualize random effects (each site)
estimate_relation(estRE,
                  include_random = TRUE) |> plot()

#Show the marginal (fixed effect) curve on top of the random effects
# ie. Plot each site trend and the overall average population trend together
estimate_relation(estRE,
                  include_random = TRUE) |>
  plot(ribbon = list(alpha = 0)) +
  geom_line(data = estimate_relation(estRE),
            mapping = aes(x = YEAR,
                          y = Predicted),
            color = "black",
            linewidth = 2) +
  theme_classic() +
  labs(y = "Subpopulation Size", x = "Year")