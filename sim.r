# simulated data from he10 designed to choose parameters and a seed
# such that the 20-city simulation from 1950 to 1964 looks
# qualitatively similar to he10

# when testing on simulated data, we need to check that the simulation is
# representative of the process we want to study. For example, if simulations
# can have 2-year or annual cycles, we may want to choose a simulation
# with 2-year cycles for comparison with the data

# also, we want to be able to test on subsets in time/space while using a
# subset of one dataset

# so, we should simulate at U=20, Tmax=1964 (CHECK THIS MATCHES HE10 EXACTLY)
# we can fit to subsets, understanding that for spatial subsets they will
# not exactly be a draw from the model on the subset.

library(spatPomp)
library(tidyverse)
library(doParallel)
library(doRNG)
registerDoParallel()
registerDoRNG(4454)

source("ibpf.R")
source("param_formats.R")
load("data/twentycities.rda")
source("he10.R")

code_version <- "sim"
code_file <- paste0(code_version,".r")
out_dir <- paste0(code_version,"/")
if(!dir.exists(out_dir)) dir.create(out_dir)

U <- 20

# build a spatPomp object with all fixed parameters
# (i.e., no unit-specific or expanded parameters)
bake(file=paste0(out_dir,"he10.rds"),
  he10(U=U,simulated=FALSE,  fixedParNames=c("alpha","iota","R0",
    "cohort","amplitude","gamma",
    "sigma","mu","sigmaSE","rho","psi","g","S_0","E_0","I_0"),
    dt=1/365)
) -> spo


# N <- 20

test <- 7

if(test==1){
## parameters comparable to a typical city for he10, generally more like London
basicParams <- c(
  alpha=0.98,
  iota=0,
  R0=30,
  cohort=0.5,
  amplitude=0.4,
  gamma=52,
  sigma=52,
  mu=0.02,
  sigmaSE=0.1,
  rho=0.5,
  psi=0.15,
  g=400,
  S_0=0.035,
  E_0=0.00005,
  I_0=0.00004
)
sim_seed <- 5
}

if(test==2){
## parameters comparable to a typical city for he10, generally more like London
basicParams <- c(
  alpha=0.99,
  iota=0,
  R0=30,
  cohort=0.5,
  amplitude=0.4,
  gamma=52,
  sigma=52,
  mu=0.02,
  sigmaSE=0.1,
  rho=0.5,
  psi=0.15,
  g=500,
  S_0=0.035,
  E_0=0.00006,
  I_0=0.00005
)
sim_seed <- 22
}

if(test==3){
## parameters comparable to a typical city for he10, generally more like London
basicParams <- c(
  alpha=0.99,
  iota=0,
  R0=30,
  cohort=0.5,
  amplitude=0.5,
  gamma=80,
  sigma=52,
#  gamma = 36.5,
#  sigma = 36.5, 
  mu=0.02,
  sigmaSE=0.09,
  rho=0.5,
  psi=0.12,
  g=500,
  S_0=0.035,
  E_0=0.00006,
  I_0=0.00005
)
sim_seed <- 27
}


if(test==4){
## parameters comparable to a typical city for he10, generally more like London
basicParams <- c(
  alpha =0.99,      iota=0,          R0=30,
  cohort=0.5,  amplitude=0.4,     gamma=80,
  sigma=52,           mu=0.02,  sigmaSE=0.08,
  rho=0.5,           psi=0.1,         g=500,
  S_0=0.035,         E_0=0.00006,   I_0=0.00005
)
sim_seed <- 701
#sim_seed <- 101
sim_seed <- 22

}

if(test==5){
## parameters comparable to a typical city for he10, generally more like London
basicParams <- c(
  alpha =0.98,      iota=0,          R0=30,
  cohort=0.5,  amplitude=0.4,     gamma=80,
  sigma=52,           mu=0.02,  sigmaSE=0.08,
  rho=0.5,           psi=0.1,         g=500,
  S_0=0.035,         E_0=0.00006,   I_0=0.00005
)
#sim_seed <- 23
sim_seed <- 24
}

if(test==6){
## parameters comparable to a typical city for he10, generally more like London
basicParams <- c(
  alpha =0.98,      iota=0,          R0=30,
  cohort=0.5,  amplitude=0.4,     gamma=80,
  sigma=80,           mu=0.02,  sigmaSE=0.05,
  rho=0.5,           psi=0.1,         g=500,
  S_0=0.035,         E_0=0.00006,   I_0=0.00005
)

sim_seed <- 14
}

if(test==7){
## parameters comparable to a typical city for he10, generally more like London
basicParams <- c(
  alpha =0.99,      iota=0,          R0=30,
  cohort=0.5,  amplitude=0.3,     gamma=80,
  sigma=52,           mu=0.02,  sigmaSE=0.05,
  rho=0.5,           psi=0.1,         g=500,
  S_0=0.036,         E_0=0.00007,   I_0=0.00006
)

sim_seed <- 20
}

set.seed(sim_seed)

measles_params <- coef(spo)
for(p in names(basicParams)) measles_params[paste0(p,1)] <- basicParams[p]

sim <- simulate(spo,params=measles_params)

ggsave(filename=paste0(out_dir,"sim-",test,".pdf"),
  plot=plot(sim,log=T,type="l")
)

saveRDS(sim,file=paste0(out_dir,"sim.rds"))

saveRDS(sim_seed,file=paste0(out_dir,"sim_seed.rds"))







