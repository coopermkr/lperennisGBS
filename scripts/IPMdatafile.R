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
