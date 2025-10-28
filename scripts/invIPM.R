### PAPER TITLE: Inverse estimation of integral projection model parameters using time series of population-level data
### DATE: 14 September 2015
### AUTHOR: Edgar J. González (edgarjgonzalez@ymail.com)

## DESCRIPTION: This code estimates unobserved structured vital rates using as input a time series of population structures and sizes.
# Known vital rates can be incorporated by placing the corresponding parameter values in the "known" column in the intervals.csv file and putting a -1 in the corresponding parameter in the PARAMETER_SECTION of the invipmADMB.tpl file.

## INPUT
## invipmADMB.dat: contains the time series (see below the section: Reading the information contained in the .dat file)
## invipmADMB.tpl: contains the code required by ADMB to calculate the composite negative log-likelihood
## invipmDREAM.cpp: contains the code required by DREAM to calculate the composite negative log-likelihood
## intervals.csv: contains the parameter values used to generate the data and the parameter intervals determining the hypercube of possible parameter values

## OUTPUT
## estimates: a warning is printed if the estimates reached a parameter bound 
## composite.ll: composite log-likelihood associated to the estimates

library(dream)
library(R2admb)
library(Rcpp)

setwd("/work/pi_brook_moyers_umb_edu/lupine")

## Set working directory (where INPUT is located)
# setwd("...") 

## Number of parameters
# this number can't be changed
n.p <- 10

## Parameter names
par.names <- paste("beta", seq(0, n.p-1), sep = "")

## Identifying the parameters for  which their values are fixed (i.e. we don't want to estimate)
# the fixed parameters have a -1 in the PARAMETER SECTION of the .tpl file
tpl <- readLines("5.ipm/LPinvipmADMB.tpl")
tpl <- tpl[grep("PARAMETER_SECTION", tpl):grep("objective_function_value", tpl)]
fixed <- c()
for (i in 1:n.p) {
  line.declare <- tpl[grep(par.names[i], tpl)[1]]
  if (grepl("-1)", line.declare, fixed = T))
    fixed <- c(fixed, i)	
}

## Hypercube limits
conf.int <- read.csv("5.ipm/LPintervals.csv", header = T)
hyper.min <- conf.int$hyper.min
hyper.max <- conf.int$hyper.max
lims <- list()
for (i in 1:n.p)
 if(all(fixed != i))
  lims <- append(lims, list(c(hyper.min[i], hyper.max[i])))
if (length(fixed) == 0)
 names(lims) <- par.names
if (length(fixed) > 0)
 names(lims) <- par.names[-fixed]

## Parameter values used to generate the data
known <- conf.int$known

## Reading the information contained in the .dat file
data <- read.table("5.ipm/LPinvipmADMB.dat", header = F, fill = T)
n.t <- data[1, 1]  # time series length
n.x <- data[2, 1]  # maximum per-year sample size (note that zeros are imputed in order to have an n.x x n.t matrix)
n.i <- data[3, 1]  # binning size
x.minf <- data[4, 1]  # size of smallest reproductive individual
t <- as.vector(data[seq(4+1, 4+n.t, 1), 1])  # observed times in the time series
d.t <- as.vector(data[seq(4+1, 4+n.t, 1), 2])  # population sizes at each observed time
D <- cbind(t, d.t)
X <- cbind(data[seq(5+n.t, dim(data)[1], 1), 1], 
  data[seq(5+n.t, dim(data)[1], 1), 2])  # individual measurements at each time
X[, 2] <- scale(X[, 2]) # if the sizes are not scaled
scale.x <- X[, 2]
X[, 1] <- X[, 1] - min(t)
t <- t-min(t)
x.p <- (2*scale.x - max(scale.x, na.rm = TRUE) - min(scale.x, na.rm = TRUE))/(max(scale.x, na.rm = TRUE) - min(scale.x, na.rm = TRUE))
a1 <- sum(x.p)/n.x
phi1 <- 2*(x.p - a1)
b1 <- sum(phi1^2)/n.x
a2 <- sum(x.p*phi1^2)/sum(phi1^2)
Mx <- max(scale.x)
mx <- min(scale.x)
min.x <- mx
max.x <- Mx

## Generating the time series of population structures from the individual measurements contained in the .dat file
n.x.t = table(X[,1])
z.i <- seq(min.x, max.x, length.out = n.i+1)
x.i <- (z.i[2:(n.i+1)] + z.i[1:(n.i)])/2
x.t <- matrix(0,nrow = n.x, ncol = n.t)
freq.i.t <- matrix(0, nrow = n.i, ncol = n.t)
for(i in 1:length(t)) {
    sub <- subset(X[, 2], X[, 1] == t[i])
    x.t[1:length(sub), i] <- sub
    freq.i.t[, i] <- hist(x.t[1:length(sub), i], breaks = z.i, plot = F)$counts    
}

## Transforming the invipmDREAM.cpp C code into an R function
sourceCpp("5.ipm/LPinvipmDREAM.cpp")

invipm1 <- function(pars) {
  # R function when all parameters are estimated
  invipm1 <- invipm(pars, n.i, x.minf, D, Mx, mx, freq.i.t)
}
invipm2 <- function(pars) {
  # R function when some parameters are fixed
  pars1 <- rep(0, n.p)
  pars1[fixed] <- known[fixed]
  pars1[-fixed] <- pars
  invipm2 <- invipm(pars1, n.i, x.minf, D, Mx, mx, freq.i.t)
}

## Run the MCMC part of the algorithm (DREAM)
# version when all parameters are estimated
n.p.f <- n.p-length(fixed)
if (length(fixed) == 0)
  results.dream <- dream(FUN = invipm1, func.type = "logposterior.density", 
   pars = lims, control = list(nseq = 3*n.p.f, burnin.length = 1e4, 
   thin.t = 1e2, REPORT = 1e4, ndraw = 1e7, Rthres = 1.01))
# version when some parameters are fixed
if (length(fixed) > 0) {
  results.dream <- dream(FUN = invipm2, func.type = "logposterior.density", 
   pars = lims, control = list(nseq = 3*n_p_f, burnin.length = 1e4, 
   thin.t = 1e2, REPORT = 1e4, ndraw = 1e7, Rthres = 1.01))
  pin <- rep(0, n.p)
  pin[fixed] <- known[fixed]
  pin[-fixed] <- results.dream$par
}
par <- results.dream$X[which(results.dream$X[, n.p.f + 1] == 
  max(results.dream$X[, n.p.f + 1])),1:n.p.f]
write.table(par, "5.ipm/LPinvipmADMB.pin", row.names = F, col.names = F)	

par <- read.table("5.ipm/LPinvipmADMB.pin", col.names = F)

## Run the gradient-based part of the algorithm (ADMB)
# admb must be installed in your computer
setup_admb("C:/Users/coope/OneDrive/Desktop/PineBarrens/RAD-Seq/lperennisGBS/5.ipm/ADMB-13.2")  # location of the admb software
compile_admb("5.ipm/LPinvipmADMB", verbose = T)
run_admb("5.ipm/LPinvipmADMB", verbose = T, profile = T, mcmc = T,
         mcmc.opts = mcmc.control(mcmcpars = "S1"))
# run with option profile = T to run likelihood profiling
# run with option mcmc = T to run mcmc

## Reading and printing output
output <- read_admb("5.ipm/LPinvipmADMB")
estimates <- output$coefficients
thresholds <- (hyper.max - hyper.min)/1000
for (i in 1:n.p)
  if (abs(estimates[i] - hyper.min[i]) < thresholds[i] | 
    abs(estimates[i] - hyper.max[i]) < thresholds[i])
    warning(paste("beta", i-1, " reached a bound", sep = ""))
estimates  # parameter estimates
(composite.ll = -output$loglik)  # composite log-likelihood associated to the estimates
(hessian = output$hes)  # if NULL the model reached a flat region
