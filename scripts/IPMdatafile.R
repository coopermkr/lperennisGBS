#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' L. perennis IPM Data
#' @date 2025-10-14
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Load library
library(tidyverse)

# Calculate size variable
census <- read_csv("data/Combined_Lupine_Sheets.csv") |>
  mutate(Size = Length*Width*Height)

# Make population estimate table
pop <- tibble(Year = c(0, 7, 10),
                  Size = c(178348, 79796, 36311))

# Filter and adjust years in census based on pop
censusFilt <- census |>
  filter(Year %in% c(2001, 2008, 2011)) |>
  mutate(Year = Year - 2001) |>
  select(Year, Size) |>
  arrange(Year) |>
  filter(!is.na(Size))

# Calculate per year sample sizes
censusFilt |>
  group_by(Year) |>
  count()

# Make parameter table
params <- tibble(Year = c(11,    # Time series length
                          1800,  # Maximum per-year sample size
                          200,   # Binning size
                          50),   # Size of smallest reproductive individual
                 Size = NA)

# Create data file
dat <- rbind(params, pop, censusFilt)

write.table(dat, file = "5.ipm/LPinvipmADMB.dat", row.names = FALSE, col.names = FALSE)
# Remove the NAs in the first four rows after this has saved just in case

filt <- census |>
  filter(Year == 2023)

#### LMS-Generated Data File
lms <- read_csv("data/LMS_master.csv") |>
  mutate(Size = clump_length*clump_width*mean_stem_height,
         Year = as.numeric(year) - min(as.numeric(year)),
         Size = as.numeric(Size)) |>
  filter(!is.na(Size),
         Size > 0) |>
  mutate(Size = log(Size)) |>
  filter(quadrat_ID %in% c("NORTH1-1", "NORTH1-2", "NORTH3.4-0", "NORTH3.4-2",
                           "SOUTHEAST6-2", "CENTER2.3-1", "CENTER2.3-2", "CENTER2.3-3",
                           "CENTER2.3-4", "CENTER2.3-5", "CENTER2.3-7", "CENTER2.3-8"))

lms |>
  filter(Size < 7) |>
  ggplot(mapping = aes(x = Size, y = total_flowering_stems)) +
  geom_point()

# Calculate per year sample sizes and populations
pop <- lms |>
  group_by(Year) |>
  count() |>
  rename(Size = n)

# Make parameter table
params <- tibble(Year = c(4,    # Time series length
                          162,  # Maximum per-year sample size
                          200,   # Binning size
                          6),   # Size of smallest reproductive individual
                 Size = NA)

# Format individual data
lmsFilt <- lms |> select(Year, Size)

# Create data file
dat <- rbind(params, pop, lmsFilt)

write.table(dat, file = "5.ipm/LPinvipmADMB.dat", row.names = FALSE, col.names = FALSE)
# Remove the NAs in the first four rows after this has saved just in case
# Go back and figure out which quadrats to keep and which to toss

# Go back in and manually log scale sizes and adjust size of smallest repro individual accordingly
