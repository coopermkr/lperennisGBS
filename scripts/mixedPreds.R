#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' L. perennis Mixed Model Simulation
#' @date 2025-10-24
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Load libraries
library(tidyverse)
library(lme4)
library(broom.mixed)
library(performance)
library(modelbased)

#### Load data ####
census <- read_csv("data/Combined_Lupine_Sheets.csv") |>
  # Calculate log scaled sizes
  mutate(Size = Length*Width) |>
  filter(!is.na(Size)) %>%
  replace(is.na(.), 0)

# Chart size distribution by year
ggplot(census, mapping = aes(x = Size)) +
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

# Calculate BaseLine observations
cenSum <- census |>
  group_by(BaseLine, Year) |>
  summarize(cov = sum(Size),
            count = n(),
            nflowers = sum(Flowers))

#cenSum |>
#  filter(Year == 2008) |>
#  write_csv(file = "ratios.csv")

# Metadata
meta <- read_csv(file = "data/censusMeta.csv")

basePop <- merge(cenSum, meta) |>
  mutate(dens = cov/(`Area(ha)`*1000000),
         popSize = count*54)


# Model population size as a function of continuous Year and RE for baseline
censusMixed <- lmer(data = basePop,
                    formula = popSize ~ Year + (1|BaseLine))

check_model(censusMixed) # Model fits well

# Assess coefficients
tidy(censusMixed) # Steady population decline of ~ 500 individuals per subpop per year
summary(censusMixed)
confint(censusMixed) # Population decline is significant negative
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

#### Ne Simulation ####
# Now write a function to simulate allele draws based on a steady stochastic 
# negative population decline from 2024

# Function 1:
# Set start parameters (pop size, pop change, fecundity estimates)
# Randomly draw from poisson distribution for new offspring population-wide (Nf)
# Calculate new total population size (N1) and number of number of previous generation survivors (Ns)
# Report Ns, Nf, N1

basePop |> filter(Year == 2023) |> summarize(N = sum(popSize))

N0 <- 26946 # Starting population size
dN <- -547*7 # Mean change in population size per year
ysd <- 2
offLambda <- 8*35 # Mean seeds produced per individual (# flowers * seeds)

## Calculate new population sizes
popSizes <- tibble(year = 1:1000,
                   N = numeric(1000))

# Set initial population size
popSizes$N[1] <- 26946

# Write a function that will iterate the population draw for 100 years
censusSize <- function(iter) {
  for (i in 1:999) {
    popSizes$N[i+1] <- max(0, (popSizes$N[i] + rnorm(1, mean = -382, sd = 2)))
  }
  popSizes |> mutate(iteration = iter) # Assign iteration group
}

# Run iterations
censusIterated <- map_dfr(.x = 1:5, .f = censusSize)

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

# Iterate
censusSimulated <- map_dfr(.x = unique(censusIterated$iteration), .f = censusSim)


# Function 2:
# Iterate over outputs from function 1
# Load in allele frequencies
# Perform binomial draws of size Nf for each allele
# Calculate new weighted allele frequencies: (p * Ns + p' * Nf)/N'
# Calculate Ne' using time step method
# Report Ne'
# Initialize population with age 0 cohort
alleles <- read_tsv("6.translocation/concord.frq") |>
  filter(N_ALLELES == 2) |>
  rename(frequencies = `{ALLELE:FREQ}`)

p0 <- str_split_i(alleles$frequencies, pattern = '\t', 1) |>
  str_split_i(pattern = ":", 2) |>
  str_split_i(pattern = '"', 1) |>
  as.numeric()

pn <- rbinom(Nf*2, 1, p)

# Store frequency history
freq_history <- NULL

# Simulation loop
NeSim <- function(iter) {
  
  subIterated <- censusIterated |>
    filter(iteration == iter)
  
  for (i in 1:nrow(subIterated)) {
    pn <- rbinom(subIterated$Nf[i]*2, 1, p)
  }
}

for (t in 1:generations) {
  
  # Binomial draw for offspring
  num_offspring <- N * (1 - survival_rate) # Assuming constant population size
  offspring_alleles <- rbinom(num_offspring * 2, 1, p_parental)
  
  # Survival of current adults
  surviving_alleles_indices <- rbinom(N * 2, 1, survival_rate) == 1
  surviving_alleles <- population_alleles[surviving_alleles_indices]
  
  # Form new population
  population_alleles <- c(surviving_alleles, offspring_alleles)
  
  # Update and store frequency
  new_p <- sum(population_alleles) / length(population_alleles)
  freq_history <- c(freq_history, new_p)
  
  # Optional: Resize population to maintain N, e.g., random sampling
  # if (length(population_alleles) > N * 2) {
  #   population_alleles <- sample(population_alleles, N * 2)
  # }
}

# Things to consider adding to the function:
# Seed bank penalty dynamic (x% of N-1, N-2, N-3, N-4, and N-5 genes establish with their corresponding p)
# Mutation



