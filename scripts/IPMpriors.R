#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' L. perennis IPM Modeling
#' @date 2025-10-13
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Load libraries
library(tidyverse)
library(lme4)
library(performance)
library(broom)


# Read LMS master datasheet
lms <- read_csv("5.ipm/LMS_master.csv") |>
  # Calculate size and fix object types
  mutate(size = clump_length*clump_width*mean_stem_height,
         reproductive = total_flowering_stems > 0,
         total_stems = as.numeric(total_stems),
         total_flowering_stems = as.numeric(total_flowering_stems)) |>
  # Fix incomplete data
  filter(size > 0,
         !is.na(total_flowering_stems))
  
# Chart size against reproductive variables
ggplot(data = lms, mapping = aes(x = size, y = total_flowering_stems)) +
  geom_point()

# Plot distributions against important variables
ggplot(data = lms, mapping = aes(x = total_stems)) +
  geom_histogram() # This should be non-zero

ggplot(data = lms, mapping = aes(x = total_flowering_stems)) +
  geom_histogram()

#### Modeling ####
### Flowering Probability ~ Size (binomial)
# GLM with random effect for plot approach
reproF <- glmer(data = lms,
             formula = reproductive ~ log(size) + 1|quadrat_ID,
             family = binomial)

check_model(reproF) # Looks ok, there's some wigglyness and overly-influential points
tidy(reproF)

# Priors from Ramula et al., 2020 https://link.springer.com/article/10.1007/s10530-020-02316-3/tables/1
# Bayesian approach
library(brms)
library(tidybayes)
library(broom.mixed)
library(bayesplot)

repro_log <- brm(data = lms,
                 family = bernoulli,
                 formula = bf(reproductive ~ lgsize,
                              lgsize ~ 0 + log(size),
                              nl = TRUE),
             
                 # Set priors from Ramula
                 prior = c(prior(logit(0.04, 0.024), nlpar = lgsize)),
                       #set_prior("normal(1.76, 0.106)", class = "intercept")),
                 iter = 2000, warmup = 1000, cores = 4, chains = 4, seed = 11)

saveRDS(repro_log, file = "5.IPM/repro_log.rds")

# Check model
plot(repro_log)
pp_check(repro_log)
brms::rhat(repro_log) |> mcmc_rhat() # At one-- convergence

# Fit with quadratic
repro_quad <- brm(data = lms,
                 family = bernoulli,
                 formula = bf(reproductive ~ lgsize + lgsize^2,
                              lgsize ~ 0 + log(size),
                              nl = TRUE),
                 
                 # Set priors from Ramula
                 prior = c(prior(logit(0.04, 0.024), nlpar = lgsize)),
                 #set_prior("normal(1.76, 0.106)", class = "intercept")),
                 iter = 2000, warmup = 1000, cores = 4, chains = 4, seed = 11)

saveRDS(repro_quad, file = "5.IPM/repro_log.quad")

# Check model
plot(repro_quad)
pp_check(repro_quad)
brms::rhat(repro_quad) |> mcmc_rhat() # At one-- convergence

# Compare Models
reloo_log <- loo(repro_log)
reloo_quad <- loo(repro_quad)
reloo_compare(loo_log, loo_quad)

# Assess model
summary(repro_log)
fixef(repro_log)

### Flowering Stems ~ Size
flstems_log <- brm(data = lms,
             family = poisson,
             formula = bf(total_flowering_stems ~ lgsize,
                          lgsize ~ 0 + log(size),
                          nl = TRUE),
             
             # Set priors from Ramula
             prior = c(prior(normal(0.32, 0.005), nlpar = lgsize)),
             #set_prior("normal(1.76, 0.106)", class = "intercept")),
             iter = 2000, warmup = 1000, cores = 4, chains = 4, seed = 11)

saveRDS(flstems_log, file = "5.IPM/flstems_log.rds")

plot(flstems_log)
pp_check(flstems_log) # Not a great fit, maybe quadratic?
brms::rhat(flstems_log) |> mcmc_rhat() # At one-- convergence

# Try with a quadratic term
flstems_quad <- brm(data = lms,
               family = poisson,
               formula = bf(total_flowering_stems ~ lgsize + lgsize^2,
                            lgsize ~ 0 + log(size),
                            nl = TRUE),
               
               # Set priors from Ramula
               prior = c(prior(normal(0.32, 0.005), nlpar = lgsize)),
               #set_prior("normal(1.76, 0.106)", class = "intercept")),
               iter = 6000, warmup = 1000, cores = 4, chains = 4, seed = 1989)

saveRDS(flstems_quad, file = "5.IPM/flstems_quad.rds")

# Check model
plot(flstems_quad)
pp_check(flstems_quad) # Still bad... Maybe a random effect?
brms::rhat(flstems_quad) |> mcmc_rhat() # At one-- convergence

# Compare models
flloo_log <- loo(flstems_log)
flloo_quad <- loo(flstems_quad) # Winner!
loo_compare(flloo_log, flloo_quad)


