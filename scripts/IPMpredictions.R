#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' L. perennis IPM Predictions
#' @date 2025-10-21
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Load libraries
library(tidyverse)
library(ipmr)

#### Load data ####
census <- read_csv("data/Combined_Lupine_Sheets.csv") |>
  # Calculate log scaled sizes
  mutate(Size = pmax(Length, Width) * Height,
         Size = log(Size)) |>
  filter(!is.infinite(Size),
         !is.na(Size)) # Drops 11 infinite or NA rows

# Chart size distribution by year
ggplot(census, mapping = aes(x = Size)) +
  geom_histogram() +
  facet_wrap(vars(Year)) # For validation: 2010/2011, 2017/2018, and 2023/2024 need to be combined


#### IPM ####
#Pull out 2023 size distribution
dist2023 <- census |>
  filter(Year == 2023)

# Set upper and lower integration limits
L <- min(dist2023$Size)*0.95
U <- max(dist2023$Size)*1.05

# Load parameter file
par <- read_table("5.ipm/badBounds/LPinvipmADMB.std") |>
  mutate(lo = value - std.dev,
         hi = value + std.dev)

# Or pull parameters from lit
par <- list(s_int = -9.17,    #Survival intercept (logit transformed)
            s_slope = 2.82,   #Survival slope (logit transformed)
            g_int = 4.12,     #Growth intercept
            g_slope = 0.26,   #Growth slope
            sd_g = 0.735,     #Growth normal distribution sd
            r_r_int = 1.76,   #Flowering probability intercept (logit transformed)
            r_r_slope = 0.04, #Flowering probability slope (logit transformed)
            r_s_int = -0.32,  #Number of inflorescence intercept (log transformed)
            r_s_slope = 0.32, #Number of inflorescence slope (log transformed)
            e = 0.12,         #Probability of plant establishment
            seed = 101,       #Number of seeds per inflorescence
            mu_fd = 1.92,     #Recruit size distribution mean
            sd_fd = 0.412)    #Recruit size distribution sd

# Compute kernels
kerns <- data.frame(logSize = seq(L, U, length.out = 1000)) |>
  mutate(surv = 1/(1+exp(-par$value[1]-par$value[2]*logSize-par$value[3]*logSize^2)),
         survLo = 1/(1+exp(-par$lo[1]-par$lo[2]*logSize-par$lo[3]*logSize^2)),
         survHi = 1/(1+exp(-par$hi[1]-par$hi[2]*logSize-par$hi[3]*logSize^2)))

# Plot survival kernel
ggplot(kerns, mapping = aes(x = logSize, y = surv)) +
  geom_point() +
  geom_point(mapping = aes(x = logSize, y = survLo)) +
  geom_point(mapping = aes(x = logSize, y = survHi))

growMean <- par$value[4] + par$value[5]*kerns$logSize
V <- pnorm(1000, mean = growMean, sd = par$std.dev[6])

# Plot growth kernel
ggplot(kerns, mapping = aes(x = logSize, y = grow)) +
  geom_point()

#### Deterministic IPM
detIPM <- init_ipm(sim_gen = "simple",
                     di_dd = "di",
                     det_stoch = "det")

# Define survival and growth kernel
detIPM <- define_kernel(proto_ipm = detIPM,
                        name = "Psimple",
                        family = "CC",
                        
                        # Set kernel function and core parameters
                        formula = s*g,
                        s = plogis(s_int + s_slope*dbh_1),
                        g = dnorm(dbh_2, mu_g, sd_g),
                        mu_g = g_int + g_slope*dbh_1,
                        
                        # Set parameter source
                        data_list = par,
                        states = list(c('dbh')),
                        
                        # Add eviction correction
                        evict_cor = TRUE,
                        evict_fun = truncated_distributions(fun = 'norm',
                                                            target = 'g'))

# Define reproduction kernel
detIPM <- define_kernel(proto_ipm = detIPM,
                        name = "Fsimple",
                        formula = r_r*r_s*r_d*p_est*n_seed,
                        family = "CC",
                        
                        # Set parameter values
                        r_r = plogis(r_r_int + r_r_slope*dbh_1), #Flowering probability
                        r_s = exp(r_s_int + r_s_slope*dbh_1),    #Number of flowering stems
                        n_seed = seed,                             #Number of seeds per stem
                        r_d = dnorm(dbh_2, mu_fd, sd_fd),        #Recruit size distribution
                        p_est = e,                                   #Establishment probability constant
                        
                        #Set parameter source
                        data_list = par,
                        states = list(c('dbh')),
                        
                        # Eviction correction (same as above)
                        evict_cor = TRUE,
                        evict_fun = truncated_distributions(fun = 'norm',
                                                            target = 'r_d'))

# Define implementation
detIPM <- define_impl(proto_ipm = detIPM,
                      make_impl_args_list(
                        kernel_names = c("Psimple", "Fsimple"),
                        int_rule     = rep("cdf", 2),
                        state_start  = rep("dbh", 2),
                        state_end    = rep("dbh", 2))) 

# Define integration
detIPM <- define_domains(proto_ipm = detIPM,
                         dbh = c(L, U, 500)) #Lower bound, Upper bound, number of bins (needs to equal n of population)

# Set starting population distribution
detIPM <- define_pop_state(proto_ipm = detIPM,
                           n_dbh = sample(dist2023$Size, 500))

# Run IPM
detModel <- make_ipm(proto_ipm = detIPM,
                    iterations = 200)

# Pull out vital rates
detLambda <- lambda(detModel)
w_ipmr      <- right_ev(detModel, iterations = 200)
v_ipmr      <- left_ev(detModel, iterations = 200)

# make_ipm_report works on either proto_ipms or made ipm objects

make_ipm_report(detModel, 
                render_output = TRUE, 
                title         = "my_simple_ipm_report") # Saves to working directory

