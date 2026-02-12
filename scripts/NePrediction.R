#''''''''''''''''''''''''''''''''''''''
#' Population size and Ne prediction
#' @date 2025-10-20
#' @author Cooper Kimball-Rhines
#''''''''''''''''''''''''''''''''''''''


library(tidyverse)


#### Ne Simulation ####
# Now write a function to simulate allele draws based on a steady stochastic 
# negative population decline from 2024

# Function 1:
# Set start parameters (pop size, pop change, fecundity estimates)
# Randomly draw from poisson distribution for new offspring population-wide (Nf)
# Calculate new total population size (N1) and number of number of previous generation survivors (Ns)
# Report Ns, Nf, N1

basePop |> filter(Year == 2023) |> summarize(N = sum(popSize))

## Calculate new population sizes
popSizes <- tibble(year = 1:1000,
                   N = numeric(1000))

# Set initial population size and number of iterations
popSizes$N[1] <- 26946
niter <- 5

# Write a function that will iterate the population draw for 100 years
censusSize <- function(iter) {
  for (i in 1:999) {
    popSizes$N[i+1] <- max(0, (popSizes$N[i] + rnorm(1, mean = -382, sd = 2)))
  }
  popSizes |> mutate(iteration = iter) # Assign iteration group
}

# Run iterations
censusIterated <- map_dfr(.x = 1:niter, .f = censusSize)

# Now write a function that calculates the number of surviving and flowering individuals in that population
censusSim <- function(iter) {
  iterSmall <<- censusIterated |> 
    filter(iteration == iter) |>
    mutate(Nf = 0,
           Ns = 0)
  
  for (i in 1:(nrow(iterSmall)-1)) {
    
    # Number of flowering individuals
    pflower <- rbinom(n = iterSmall$N[i], size = 1, prob = 0.3) |> sum()
    
    # Of the flowering individuals, randomly draw the number of offspring and multiply by multi-year establishment rate
    iterSmall$Nf[i] <- rpois(pflower, 280) |>
      sum()*.001
    
    
    # Calculate number of survivors
    iterSmall$Ns[i] <- max(0, round(max(0,iterSmall$N[i+1])) - round(max(0,iterSmall$Nf[i])))
  }
  return(iterSmall)
}

# Iterate and round the outputs
censusSimulated <- map_dfr(.x = 1:niter, .f = censusSim) |>
  mutate(Nf = round(Nf),
         Ns = round(Ns),
         N = round(N))

write_csv(censusSimulated, file = "6.translocation/censusSimulated.csv")

# Effective Population Size Calculations:
# Calculate Ne' using time step method
# Report Ne'
# Initialize population with age 0 cohort

# Load in allele frequencies
alleles <- read_tsv("6.translocation/concord.frq") |>
  filter(N_ALLELES == 2) |>
  rename(frequencies = `{ALLELE:FREQ}`) |>
  slice_sample(n = 10)

NeSim <- function(iter) {
  # Iterate over outputs from function 2
  smallSim <- censusSimulated |>
    filter(iteration == 1) |>
    mutate(Ne = 0)
  
  # Reformat frequencies into a wide dataframe same length as simulated census
  p <- str_split_i(alleles$frequencies, pattern = '\t', 1) |>
    str_split_i(pattern = ":", 2) |>
    str_split_i(pattern = '"', 1) |>
    as.numeric() %>%
    data.frame(p = .) |>
    t() %>%
    data.frame(p = .,
               n = 1:nrow(smallSim)) |>
    select(!n)
  
  # One for loop iterates over each SNP (columns)
  for (k in 1:ncol(p)) {
    # One loop iterates over each generation (rows)
    for (h in 2:nrow(p)) {
      # Perform binomial draws of size Nf for each allele
      pf <- rbinom(1, (smallSim$Nf[h-1]*2), p[h-1,k])
      # And calculate number of alleles from surviving plants
      ps <- (p[h-1, k] * smallSim$Ns[h-1] * 2)
      
      # Calculate new weighted allele frequencies and record
      p[h,k] <- (pf + ps)/(2*smallSim$N[h])
    }
  }
  
  # Copy the allele frequency table minus the start point
  dp <- p[-1,]
  
  # Calculate Ne from the alleles and return the census simulation table
  for (i in 2:nrow(p)) {
    for (j in 1:(ncol(p)-1)) {
      pvar <- p[i,j] - p[i-1,j]
      pbar <- (p[i-1,j] + p[i,j]) / 2
      dp[i-1,j] <- pvar/(pbar * (1-pbar))
    }
  }
  yrs = 2:nrow(p)
  
  dp |>
    rowwise() |>
    mutate(pF = mean(c_across(where(is.numeric)), na.rm = TRUE),
           varF = (2*pF^2)/ncol(p),
           iteration = iter) |>
    ungroup() |>
    mutate(year = 2:nrow(p)) |>
    select(pF, varF, year, iteration)
}

NeSimulated <- map_dfr(.x = 1:niter, .f = NeSim)

# Join back to the simulated census table and calculate Ne

NeN <- NeSimulated |>
  merge(censusSimulated, by = c("year", "iteration"), all.y = T) |>
  mutate(invNeI = 2*(pF - 1/20 + 1/N),
         invNeII = 2*(pF - 1/20),
         varInvNe = 4 * varF) |>
  arrange(iteration, year)

write_csv(NeN, "6.translocation/NeSimulated.csv")

#' Things to come back to:
#' Need to actually figure out how fast the population is declining
#' Need to get updated starting population size estimate
#' Need to adjust allele frequencies to account for transplant from AL-- should probably use a new script

# Things to consider adding to the function:
# Seed bank penalty dynamic (x% of N-1, N-2, N-3, N-4, and N-5 genes establish with their corresponding p)
# Mutation