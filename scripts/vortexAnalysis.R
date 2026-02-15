#'''''''''''''''''''''''''''''''''''''
#' Vortex Parameters
#' @author Cooper Kimball-Rhines
#' @date 2026-01-13
#'''''''''''''''''''''''''''''''''''''


# Load libraries
library(tidyverse)

# Load in translocation vortex runs
tNe <- read_delim("5.translocation/VOutput/translocations/VOutput/collatedNe.txt",
                         col_names = FALSE) |>
  rename("Population" = X1, "Year" = X2, "Ne" = X3) |>
  na.omit()

ggplot(tNe, mapping = aes(x = Year, y = Ne, color = Population)) +
  geom_point()

