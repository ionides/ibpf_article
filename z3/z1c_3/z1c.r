# z2 is a new attempt at maximizing likelihood for data & model of he10

devel_version <- "z1"
old_version <- "NA" # used to read in starting values when continuation=TRUE
run_version <- "c"
code_version <- paste0(devel_version,run_version)
run_level <- 3
code_file <- paste0(code_version,".r")
out_dir <- paste0(code_version,"_",run_level,"/")
if(!dir.exists(out_dir)) dir.create(out_dir)

continuation <- FALSE # set to be false when starting a new search
if(continuation){
  old_code_version <- paste0(old_version,run_version)
  old_out_dir <- paste0(old_code_version,"_",run_level,"/")
  old_results <- read.csv(file=paste0(old_out_dir,"global.csv"))
}

sim_data <- FALSE # fit simulated data rather than the actual data

estimated <- switch(run_version,
  a="mostly shared",
  b="all unit-specific",
  c="he10",
  d="mostly shared"
) # he10 has all unit-specific with immigration but not coupling

# estimated time on greatlakes at run level 3: 24hr with 4000 particles

continuation <- FALSE # set to be false when starting a new search

code_file <- paste0(code_version,".r")
out_dir <- paste0(code_version,"_",run_level,"/")
if(!dir.exists(out_dir)) dir.create(out_dir)
if(continuation){
  # old results used for starting the new search
  old_results <- read.csv(file=paste0(out_dir,"global.csv"))
}

rw_r_vec <- c(0.02,0.01,0.005, 0.0025, 0.00125)
rw_ivp_vec <- c(0.05,0.025,0.01, 0.005, 0.0025)
# set rw.sd smaller for alpha since it has uncertainty +/- 5 percent
# whereas +/- 50 percent is more generally the case.
rw_alpha_frac <- 0.1

# rate of autoregression of means (zero for uncoupled)
spat_reg <- switch(run_version,
  "a" = 0.1,
  "b" = 0.1,
  "c" = 0.1,
  "d" = 0.001
)

experiment <- 3
rw_r <- rw_r_vec[experiment]
rw_ivp <- rw_ivp_vec[experiment]
rw_alpha <- rw_r * rw_alpha_frac

print(paste0("code version =",code_version))
print(paste0("run_level = ",run_level))
print(paste0("experiment = ",experiment))

# save various things of possible later interest
# when doing the main computations
# set to FALSE when re-running to analyze the baked results
save_session_data <- TRUE
# save_session_data <- FALSE

# do_bake = FALSE helps debug when bake triggers recomputation unintentionally
#do_bake <- FALSE 
do_bake <- TRUE

# compute likelihood at initial guesses when we want to
# check local optimization: how well does each search do from where it starts?
compute_start_likelihood <- TRUE
# compute_start_likelihood <- experiment==1
# compute_start_likelihood <- FALSE

# standard pomp pseudocode notation: U = units, N = time points,
# J = particles, M = filter iterations.

if(run_level==1){
  U <- 20
  N <- 730
  J <- 20
  M <- 2
  Nblocks <- U
  lik_eval_Reps <- 2
  lik_eval_J <- J
  UniqueStart <- 2
  RepStart <- 1
  # Reps is UniqueStart * RepStart
} else if(run_level==2){  
  U <- 20
  N <- 730
  J <- 1000
  M <- 50
  Nblocks <- U
  lik_eval_Reps <- 5
  lik_eval_J <- 2000
  UniqueStart <- 18
  RepStart <- 2
  # Reps is UniqueStart * RepStart  
} else if(run_level==3){
  U <- 20
  N <- 730
  J <- 4000
  M <- 100
  Nblocks <- U
  lik_eval_Reps <- 10
  lik_eval_J <- 8000
  UniqueStart <- 18
  RepStart <- 2
  # Reps is UniqueStart * RepStart  
}

library(spatPomp)
library(tidyverse)
source("../ibpf.R")
source("../param_formats.R")
load("../data/twentycities.rda")
source("../he10.R")

library(doParallel)
cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset=NA))
if(is.na(cores)) cores <- detectCores()  
registerDoParallel(cores)

library(doRNG)
registerDoRNG(4454)

if(save_session_data) file.copy(code_file,
  paste0(out_dir,code_file),overwrite=TRUE)
sbat_file <- paste0(code_version,"_",run_level,".sbat")

## parameters comparable to a typical city for he10
## selected and tested in sim.r
sim_test <- readRDS(file="../sim/sim.rds")

basicParams <- coef(sim_test)
names(basicParams) <- gsub("[1-9]+[0-9]*$","",names(basicParams))
  parNames <- names(basicParams)

if(estimated == "mostly fixed"){
  sharedParNames <- c("R0","psi")
  unitParNames <- c("rho","S_0")
  estParNames <- c(sharedParNames,unitParNames)
  fixedParNames <- setdiff(parNames,estParNames)
}

if(estimated == "mostly shared"){
  sharedParNames <- c("alpha","R0","psi","g","sigma","gamma","amplitude","cohort","sigmaSE")
  unitParNames <- c("rho","S_0","E_0","I_0")
  estParNames <- c(sharedParNames,unitParNames)
  fixedParNames <- setdiff(parNames,estParNames)
}

if(estimated == "all plausible parameters shared"){
  # all the parameters estimated by he10, table 2, but shared when
  # that makes mechanistic sense
  sharedParNames <- c("alpha","R0","g","sigma","gamma","amplitude","cohort")
  unitParNames <- c("sigmaSE","S_0","E_0","I_0","rho","psi")
  estParNames <- c(sharedParNames,unitParNames)
  fixedParNames <- setdiff(parNames,estParNames)
}

if(estimated == "all unit-specific"){
  # all the parameters estimated by he10, table 2, except for immigration
  fixedParNames <- c("mu","iota")
  sharedParNames <- NULL
  unitParNames <- setdiff(parNames,fixedParNames)
  estParNames <- c(sharedParNames,unitParNames)
}

if(estimated == "he10"){
  # all the parameters estimated by he10, table 2, including immigration
  fixedParNames <- c("mu","g")
  sharedParNames <- NULL
  unitParNames <- setdiff(parNames,fixedParNames)
  estParNames <- c(sharedParNames,unitParNames)  
}

ivpParNames <- c("S_0","E_0","I_0")
ivpEstParNames <- intersect(ivpParNames,estParNames)
regEstParNames <- setdiff(estParNames,ivpParNames)

Tmax <- min(1950 + N/52, 1964)

spo <-  he10(U=U,Tmax=Tmax,
    sharedParNames=sharedParNames,
    fixedParNames=fixedParNames,
    unitParNames=unitParNames,
    basic_params=basicParams,
    dt=1/365
  )

if(sim_data){
  sim_seed <- readRDS("../sim/sim_seed.rds")
  set.seed(sim_seed)
  spo <- simulate(spo)
  saveRDS(file=paste0(out_dir,"model.rds"),spo)
  # For the case where we are matching the he10 dimensions, the result
  # should coincide with sim_test. Otherwise, replace the simulated data
  # with a subset of sim_test to obtain a nicely nested sequence of problems
  if( all(dim(obs(spo))==dim(obs(sim_test))) ){
   if(any(obs(spo)!=obs(sim_test))) stop("spo and sim_test do not match")
  } else {
   spo@data <- sim_test@data[1:nrow(spo@data),1:ncol(spo@data)]
  }
}

contracted_params <- contract_params(
  coef(spo), average=FALSE,U=U,
  fixedParNames=fixedParNames,
  sharedParNames=sharedParNames,
  unitParNames=unitParNames
)

estParNames_expanded <- unlist(lapply(estParNames,function(x)paste0(x,1:U)))
regEstParNames_expanded <- unlist(lapply(regEstParNames,function(x)paste0(x,1:U)))
ivpEstParNames_expanded <- unlist(lapply(ivpEstParNames,function(x)paste0(x,1:U)))
fixedParNames_expanded <- paste0(fixedParNames,1)

reg_rw.sd <- rep(list(rw_r),times=length(regEstParNames_expanded))
names(reg_rw.sd) <- regEstParNames_expanded
if("alpha"%in%estParNames) reg_rw.sd[paste0("alpha",1:U)] <- rw_alpha

 
ivp_rw.sd <- lapply(ivpEstParNames_expanded,function(x)expression(ivp(rw_ivp)))
names(ivp_rw.sd) <- ivpEstParNames_expanded
measles_rw.sd <- do.call(rw.sd,c(reg_rw.sd,ivp_rw.sd))

all_units = seq_len(length(unit_names(spo)))
nblocks = Nblocks
block_list = split(all_units, sort(all_units %% nblocks))
block_list <- lapply(block_list, as.integer)

if(!continuation){
  fixed_params <- coef(spo,fixedParNames_expanded)
  if(estimated == "he10") fixed_params["g1"] <- 0
  estParTransCenter <- coef(spo,estParNames_expanded,transform=TRUE)
  if(estimated == "he10") estParTransCenter[paste0("iota",1:U)] <- log(0.1)
  set.seed(53125455)
  runif_design(
    lower=estParTransCenter -0.1,
    upper=estParTransCenter +0.1,
    nseq=UniqueStart*RepStart
  ) -> guesses_transformed
  guesses <- as.data.frame(t(apply(guesses_transformed,1,function(p,po,fi) {
      partrans(object=spo,unlist(c(p,fi)),dir="fromEst")
    }, po=spo,fi=fixed_params)))
} else {
  start_params <- old_results[,1:(length(coef(spo)))]
  guesses <- matrix(NA,0,0)
  for(i in 1:UniqueStart) {
    for(j in 1:RepStart) guesses <- rbind(guesses,start_params[i,])
  }
}

set.seed(0)
bake(file=paste0(out_dir,"global.rds"),{
  registerDoRNG(1270401374)
  foreach(guess=iter(guesses,by="row"), .combine=rbind) %dopar% {
    library(pomp)
    library(tidyverse)
    ibpf(spo,
      params=guess,
      sharedParNames=sharedParNames,
      unitParNames=unitParNames,
      Nbpf=M,
      spat_regression=spat_reg,
      Np=J,
      rw.sd=measles_rw.sd,
      cooling.fraction.50=0.5,
      block_list=block_list,
      tol=1e-10 
    ) -> mf
    replicate(
      lik_eval_Reps,
      mf %>% bpfilter(Np=lik_eval_J,
        block_list=block_list) %>% logLik()
    ) %>%
      logmeanexp(se=TRUE) -> ll
    mf %>% coef() %>% dplyr::bind_rows() %>%
      dplyr::bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- getDoParWorkers()
  results
}) -> results


if(compute_start_likelihood){
  bake(file=paste0(out_dir,"initial_lik.rds"),{
    registerDoRNG(3889100)
      foreach(guess=iter(guesses,by="row"), .combine=rbind) %dopar% {
      library(pomp)
      library(tidyverse)
      replicate( lik_eval_Reps,
        logLik(bpfilter(spo,params=guess,Np=lik_eval_J,
          block_list=block_list))
      ) %>% logmeanexp(se=TRUE) -> ll
    c(start_loglik=ll[1],start_loglik.se=ll[2])
    } -> start_loglik
    attr(start_loglik,"ncpu") <- getDoParWorkers()
    start_loglik
  }) -> start_loglik
  results <- cbind(results,start_loglik)
  print(start_loglik)
}

results <-  filter(results,is.finite(loglik))
t_global <- attr(results,"system.time")
ncpu_global <- attr(results,"ncpu")

print(t_global)
print(ncpu_global)

#  read_csv("measles_params.csv") %>%
results %>%
    filter(is.finite(loglik)) %>%
    arrange(-loglik) %>%
    write_csv(file=paste0(out_dir,"global.csv"))
    

print(results[,c("loglik","loglik.se",estParNames_expanded[1:5])])

estParNames_contracted <- c(sharedParNames,as.vector(sapply(unitParNames,function(x)paste0(x,1:U))))

theta_true <- contract_params(coef(spo),
  fixedParNames=fixedParNames,
  sharedParNames=sharedParNames,
  unitParNames=unitParNames,
  U=U
)[estParNames_contracted]

results %>%
  filter(loglik>max(loglik)-5000) %>% 
  bind_rows(guesses) %>%
  mutate(type=if_else(is.na(loglik),"guess","result")) %>%
  arrange(type) -> all


pdf(file=paste0(out_dir,"pairs.pdf"),width=12,height=8)
if(estimated == "mostly fixed"){
  pairs(~loglik+rho1+psi1+R01+S_01, data=all, pch=16, cex=0.8,
    col=ifelse(all$type=="guess",grey(0.5),"red"))
}

if(estimated == "mostly shared"){
  pairs(~loglik+rho1+psi1+R01+alpha1+g1+sigma1+amplitude1+sigmaSE1+
    cohort1+gamma1+S_01+E_01+I_01, data=all, pch=16, cex=0.8,
    col=ifelse(all$type=="guess",grey(0.5),"red"))
}

dev.off()

if(save_session_data){
  if(sim_data){
    plot(spo,log=T,type="l")
    ggsave(filename=paste0(out_dir,"sim.pdf"))
  }
  cat(capture.output(sessionInfo()),
    file=paste0(out_dir,"sessionInfo.txt"),sep="\n")

  # keep a copy of the sbat file if it exists
  if(file.exists(sbat_file)) file.copy(sbat_file,
    paste0(out_dir,sbat_file),overwrite=TRUE)		   

}

