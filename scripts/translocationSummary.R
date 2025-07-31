#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' L. perennis Translocation Summary
#' @date 2025-07-30
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


# Load libraries
library(tidyverse)
library(cgsim)

load("4.translocation/vtDrift.RData")
names(vtDrift)
dimnames(vtDrift$gen_div) # Heterozygosity, allelic richness, and mean allelic diversity
dimnames(vtDrift$n_tracking)
#Ne_vp: Ne for variable pop size (9)
#Ne_vr: Ne for variable reproductive success (10)
#Nc: Census population size (not really reliable) (2)

# Calculate expected heterozygosity retention
exp <- gd_ret_exp(out = vtDrift, ne = "Nc")
exp <- apply(exp, 1, mean)

# Calculate predicted heterozygosity retention
pred <- gd_ret_pred(out = vtDrift, gd = "He") 
het <- cbind(pred[, 1], exp) 

# Summarize Ne for variable populations
nevp <- data.frame(vtDrift$n_tracking[,,9]) |>
  # Add year column
  mutate(year = 0:100) |>
  # Pivot replicates
  pivot_longer(names_to = "reps",
               values_to = "vp",
               cols = 1:5) |>
  # Calculate mean and standard error of Ne
  group_by(year) |>
  summarize(mean_vp = mean(vp),
            se_vp = sd(vp)/sqrt(n()))

# Summarize Ne for variable reproduction
nevr <- data.frame(vtDrift$n_tracking[,,10]) |>
  # Add year column
  mutate(year = 0:100) |>
  # Pivot replicates
  pivot_longer(names_to = "reps",
               values_to = "vr",
               cols = 1:5) |>
  # Calculate mean and standard error of Ne
  group_by(year) |>
  summarize(mean_vr = mean(vr),
            se_vr = sd(vr)/sqrt(n()))

# Plot
ggplot(data = nevp,
       mapping = aes(x = year, y = mean_vp)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_vp - se_vp, ymax = mean_vp + se_vp)) +
  theme_classic()

# Create a null table
ne <- data.frame(year = integer(),
                 sim = character(),
                 mean_vp = numeric(),
                 se_vp = numeric(),
                 mean_vr = numeric(),
                 se_vr = numeric())

# Write a function to process all of the simulations and combine Ne
necalcs <- function(simname) {
  # Load in saved simulation file
  simstring <- deparse(substitute(simname))
  
  nevp <- data.frame(simname$n_tracking[,,9]) |>
    # Add year column
    mutate(year = 0:100,
           sim = simstring) |>
    # Pivot replicates
    pivot_longer(names_to = "reps",
                 values_to = "vp",
                 cols = 1:5) |>
    # Calculate mean and standard error of Ne
    group_by(year, sim) |>
    summarize(mean_vp = mean(vp),
              se_vp = sd(vp)/sqrt(n()))
  
  # Summarize Ne for variable reproduction
  nevr <- data.frame(simname$n_tracking[,,10]) |>
    # Add year column
    mutate(year = 0:100,
           sim = simstring) |>
    # Pivot replicates
    pivot_longer(names_to = "reps",
                 values_to = "vr",
                 cols = 1:5) |>
    # Calculate mean and standard error of Ne
    group_by(year, sim) |>
    summarize(mean_vr = mean(vr),
              se_vr = sd(vr)/sqrt(n()))
  
  # Merge and store
  tmp <- merge(nevp, nevr) |>
    arrange(year)

  ne <<- rbind(tmp, ne)
}

# Load in simulations
load("4.translocation/vtDrift.RData")
load("4.translocation/vtBad.RData")
load("4.translocation/vtTrans.RData")
load("4.translocation/vtTransBad.RData")

# Run the function on each loaded simulation
necalcs(vtTransBad)

# Plot each simulation together
ggplot(data = ne,
       mapping = aes(x = year, y = mean_vp, color = sim)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_vp - se_vp, ymax = mean_vp + se_vp)) +
  theme_classic()
