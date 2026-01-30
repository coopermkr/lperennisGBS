#'''''''''''''''''''''''''''''''''''''
#' Vortex Parameters
#' @author Cooper Kimball-Rhines
#' @date 2026-01-13
#'''''''''''''''''''''''''''''''''''''

library(tidyverse)

#' List of Parameters and Definitions:
#' E: Critical population size to be considered extinct
#' LE: # lethal equivalents for inbreeding depression (optional)
#' YO: Age of first offspring
#' ML: Maximum lifespan
#' MP: Maximum possible progeny per individual per year
#' DD: Density dependence penalty (optional)
#' RP/RSD: % adult females in breeding pool and SD
#' O/OSD: Mean and SD number of offspring per adult per year
#' D/DSD: Mortality rates and SD by age until year 4
#' C: Frequency of catastrophe
#' CR/CS Reproductive and survival penalty of catastrophe
#' S: Selfing rate
#' N0: Initial census population size
#' K/KSD: Carrying capacity and SD
#' NT: Number of individuals translocated each year
#' P: Starting allele frequencies by population and locus

# Parameter Values:
E <- 5 # Chosen somewhat randomly
LE <- NA
YO <- 3 # Observations in field and USDA page
ML <- 10 # Guess based on observations and pers. comm.
MP <- 1000 # Number of seeds, not all viable
DD <- NA
RP <- 0.657 # From LMS
RSD <- 0.169 # From LMS
O <- 128
OSD <- 155
D1 <- 1-0.72*0.265 #mortality
DSD1 <- 0.05 # 20% of mean (same as variance below)
D <- 1-0.55 #mortality
DSD <- 0.1
C <- 0.02
CR <- 0.66
CS <- 0.5
S <- 0.1

## Population specific values:
# Vermont
N0 <- 12 # From census data
K <- 80  # Maximum population historically measured
KSD <- ? # Sensitivity analysis
NT <- 12 # From what Grace is doing

# Concord


#### Calculating Parameters
# Read LMS master datasheet
lms <- read_csv("data/LMS_master.csv") |>
  # Calculate size and fix object types
  mutate(size = max(clump_length, clump_width)*mean_stem_height,
         reproductive = as.numeric(total_flowering_stems > 0),
         total_stems = as.numeric(total_stems),
         total_flowering_stems = as.numeric(total_flowering_stems)) |>
  # Fix incomplete data
  filter(size > 0,
         !is.na(total_flowering_stems))

#' MP: Maximum possible progeny per individual per year
# Maximum number of progeny = maxFloweringStems*Seeds per inflorescence
# https://link.springer.com/article/10.1007/s00497-005-0250-3
MP <- max(lms$total_flowering_stems)*15.83

#' RP/RSD: % adult females in breeding pool and SD
# % of females in breeding pool by year to calculate SD
lms |> group_by(year) |>
  count(total_flowering_stems > 0) |> # Tally the number of reproductive and non-reproductive individuals by quadrat
  pivot_wider(names_from = `total_flowering_stems > 0`, values_from = n) |> # Pivot
  mutate(`FALSE` = replace_na(`FALSE`, 0), `TRUE` = replace_na(`TRUE`, 0)) |> # Replace NA with 0
  ungroup() |>
  summarize(RP = mean(`TRUE`/(`TRUE`+`FALSE`)),
            RSD = sd(`TRUE`/(`TRUE`+`FALSE`))) # Calculate overall mean and SD
 
#' O/OSD: Mean and SD number of offspring per adult per year
# Ramula (poly) says each inflorescence produces 101 seeds
lms |> filter(total_flowering_stems > 0) |>
  mutate(seeds = total_flowering_stems*101) |>
  summarize(O = mean(seeds),
            OSD = sd(seeds))

# Shi et al. say 15.83 seeds per inflorescence in perennis
lms |> filter(total_flowering_stems > 0) |>
  mutate(seeds = total_flowering_stems*15.83) |>
  summarize(O = mean(seeds),
            OSD = sd(seeds))

#' D/DSD: Mortality rates and SD by age until year 4
# FOR SEED BANK POPULATIONS:
# Shi et al. say 84% (SE = 0.04, n = 14) of seeds germinate 
# and 95% (SE = 0.02, n = 14) of those survive to 5 months
D1 <- 0.84*0.95
DSD1 <- sqrt((0.04^2)/14 + (0.02^2)/14)

# FOR WILD POPULATIONS:
# This is highly variable depending on management
# Mackay says 6.32% germination for unscarified seeds in lab
# USDA (Trudeau) says 69-75% germination in burned with liter/surface placed seeds
# For 1-year germinated seedling survival: 26.5% in NY, 60% in OH, and 20% in Ontario
# I would go with the NY value Zaremba reported bc OH manipulated with fire
D1 <- 0.72*0.265
DSD1

# From KBB Recovery Plan 2003: "survival of seedlings was about 50-60% per year"
D <- 0.55
DSD <- 0.1

#' Catastrophe:
# Number of droughts that had 2in of rain in NH in the last 100 years:
precip1 <- read_csv("data/precipConcord_2005_2014.csv") |>
  mutate(year = as.numeric(str_split_i(DATE, "-", 1)),
         month = as.numeric(str_split_i(DATE, "-", 2)),
         PRCP = replace_na(PRCP, 0))

precip2 <- read_csv("data/precipConcord_2015_2025.csv") |>
  mutate(year = as.numeric(str_split_i(DATE, "-", 1)),
         month = as.numeric(str_split_i(DATE, "-", 2)),
         PRCP = replace_na(PRCP, 0))

spring <- rbind(precip1, precip2) |>
  filter(month %in% 3:5,
         NAME == "CONCORD MUNICIPAL AIRPORT, NH US") |>
  group_by(year) |>
  summarize(spPrecip = sum(PRCP))

# Vermont experienced a D2 drought in 2024 that resulted in 50% mortality.
# D2 droughts happen 2% of the time now or once every 50 years
C <- 0.02
Cs <- 0.5

# But under future climate change, these may get worse
C <- 0.04
C <- 0.08


# Water availability affects percent germination. For instance, seeds that
# received 11 inches (28 cm, ambient) or 14 inches (35 cm, wet) of water over 3
# months exhibited 92% germination, while only 62% of seeds limited to 2 inches
# (6 cm, dry) of water germinated.
CR <- 0.62 # Reproductive penalty



#' Beneficial management: Could add in to imitate burns if data


#### Sensitivity Test
# No calculation for K or KSD, so we do a sensitivity test
kst <- read_csv("5.translocation/VOutput/kst.csv") |>
  select(PE, GeneDiv, meanNe, SV1, SV2) # Gene and Ne not done because run as pop

ggplot(kst, mapping = aes(x = SV1, y = PE, color = SV2)) +
  geom_point() + theme_light()

ggplot(kst, mapping = aes(x = SV2, y = PE)) +
  geom_point() + theme_light()

# It seems like K does not have an effect on P(E) once it is above like 75 individuals
# So we should set K to the historic population max and do another ST on KSD
# KSD
ksd <- read_csv("5.translocation/VOutput/KSDst.csv") |>
  select(PE, GeneDiv, meanNe, SV1)

ggplot(ksd, mapping = aes(x = SV1, y = PE)) +
  geom_point() + theme_light()

# Simulation is robust to this setting, little to no effect.
