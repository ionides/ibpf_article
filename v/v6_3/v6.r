
# experimenting with values of spat_regression.
# such experiments could perhaps be carried out on a smaller example,
# but we want to get experience relevant to the "target problem" of
# data analysis on the he10 dataset.

# this experiment corresponds to a broad search (starting values +/- 1 log
# unit from the truth, rw.sd = 0.02). we could check later whether similar
# conclusions hold for a local search

# results on greatlakes, 
# reps cores time(min) memory(GB) dir
# 36   36    12:23     11.88/180    v6_3  with likelihood evaluation
# 36   36    10:48     11/48/180    v6_3 without likelihood evaluation

code_version <- "v6"
run_level <- 3

experiment <- 1
spat_regression <- c(0, 1/64, 1/32, 1/16, 1/8, 1/4, 1/2, 1) 

print(paste0("run_level = ",run_level))
print(paste0("experiment = ",experiment))

truth_loglik <- -38882.7 # from w1.r 

# save various things of possible later interest
# when doing the main computations
# set to FALSE when re-running to analyze the baked results
# save_session_data <- TRUE
save_session_data <- experiment==1
# save_session_data <- FALSE

# do_bake = FALSE helps debug when bake triggers recomputation unintentionally
#do_bake <- FALSE 
do_bake <- TRUE

# compute likelihood at initial guesses when we want to
# check local optimization: how well does each search do from where it starts?
# compute_start_likelihood <- TRUE
compute_start_likelihood <- experiment==1
# compute_start_likelihood <- FALSE

# standard pomp pseudocode notation: U = units, N = time points,
# J = particles, M = filter iterations.
# Also, Reps = replicates
# estimated is in ("mostly shared","mostly fixed",
#   "all plausible parameters shared","all unit-specific")

if(run_level==1){
  estimated <- "mostly shared"
  U <- 5
  N <- 10
  J <- 50
  M <- 2
  Nblocks <- U
  lik_eval_Reps <- 2
  lik_eval_J <- J
  Reps <- 2
} else if(run_level==2){  
  estimated <- "mostly shared"
  Reps <- 36
  U <- 20
  N <- 730
  J <- 1000
  M <- 50
  Nblocks <- U
  lik_eval_Reps <- 5
  lik_eval_J <- 2000
} else if(run_level==3){
  estimated <- "mostly shared"
  Reps <- 36
  U <- 20
  N <- 730
  J <- 2000
  M <- 100
  Nblocks <- U
  lik_eval_Reps <- 10
  lik_eval_J <- 5000
}

library(spatPomp)
library(tidyverse)
source("ibpf.R")
source("param_formats.R")
source("he10.R")

library(doParallel)
cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset=NA))
if(is.na(cores)) cores <- detectCores()  
registerDoParallel(cores)

library(doRNG)
registerDoRNG(4454)

code_file <- paste0(code_version,".r")
out_dir <- paste0(code_version,"_",run_level,"/")
sbat_file <- paste0(code_version,"_",run_level,".sbat")

if(!dir.exists(out_dir)) dir.create(out_dir)
if(save_session_data) file.copy(code_file,
  paste0(out_dir,code_file),overwrite=TRUE)


## parameters comparable to a typical city for he10, generally more like London
## selected and tested in sim.r
sim_test <- readRDS(file="sim/sim.rds")

basicParams <- coef(sim_test)
names(basicParams) <- gsub("[1-9]+[0-9]*$","",names(basicParams))

if(estimated == "mostly fixed"){
  parNames <- names(basicParams)
  sharedParNames <- c("R0","psi")
  unitParNames <- c("rho","S_0")
  estParNames <- c(sharedParNames,unitParNames)
  fixedParNames <- setdiff(parNames,estParNames)
  ivpParNames <- c("S_0","E_0","I_0")
  ivpEstParNames <- intersect(ivpParNames,estParNames)
  regEstParNames <- setdiff(estParNames,ivpParNames)
}

if(estimated == "mostly shared"){
  parNames <- names(basicParams)
  sharedParNames <- c("alpha","R0","psi","g","sigma","gamma","amplitude","cohort","sigmaSE")
  unitParNames <- c("rho","S_0","E_0","I_0")
  estParNames <- c(sharedParNames,unitParNames)
  fixedParNames <- setdiff(parNames,estParNames)
  ivpParNames <- c("S_0","E_0","I_0")
  ivpEstParNames <- intersect(ivpParNames,estParNames)
  regEstParNames <- setdiff(estParNames,ivpParNames)
}

if(estimated == "all plausible parameters shared"){
  # all the parameters estimated by he10, table 2, but shared when
  # that makes mechanistic sense
  parNames <- names(basicParams)
  sharedParNames <- c("alpha","R0","g","sigma","gamma","amplitude","cohort")
  unitParNames <- c("sigmaSE","S_0","E_0","I_0","rho","psi")
  estParNames <- c(sharedParNames,unitParNames)
  fixedParNames <- setdiff(parNames,estParNames)
  ivpParNames <- c("S_0","E_0","I_0")
  ivpEstParNames <- intersect(ivpParNames,estParNames)
  regEstParNames <- setdiff(estParNames,ivpParNames)
}

if(estimated == "all unit-specific"){
  # all the parameters estimated by he10, table 2, except for immigration
  parNames <- names(basicParams)
  fixedParNames <- c("mu","iota")
  sharedParNames <- NULL
  unitParNames <- setdiff(parNames,fixedParNames)
  estParNames <- c(sharedParNames,unitParNames)
  ivpParNames <- c("S_0","E_0","I_0")
  ivpEstParNames <- intersect(ivpParNames,estParNames)
  regEstParNames <- setdiff(estParNames,ivpParNames)
}

Tmax <- min(1950 + N/52, 1964)

sim_seed <- readRDS("sim/sim_seed.rds")

sim <- bake(file=paste0(out_dir,"model.rds"),{
  set.seed(sim_seed)
  he10(U=U,Tmax=Tmax,
    sharedParNames=sharedParNames,
    fixedParNames=fixedParNames,
    unitParNames=unitParNames,
    simulated=TRUE,
    basic_params=basicParams
  )
})

# For the case where we are matching the he10 dimensions, the result
# should coincide with sim_test. Otherwise, replace the simulated data
# with a subset of sim_test to obtain a nicely nested sequence of problems
if( all(dim(obs(sim))==dim(obs(sim_test))) ){
 if(any(obs(sim)!=obs(sim_test))) stop("sim and sim_test do not match")
} else {
 sim@data <- sim_test@data[1:nrow(sim@data),1:ncol(sim@data)]
}

contracted_params <- contract_params(
  coef(sim), average=FALSE,U=U,
  fixedParNames=fixedParNames,
  sharedParNames=sharedParNames,
  unitParNames=unitParNames
)

estParNames_expanded <- unlist(lapply(estParNames,function(x)paste0(x,1:U)))
regEstParNames_expanded <- unlist(lapply(regEstParNames,function(x)paste0(x,1:U)))
ivpEstParNames_expanded <- unlist(lapply(ivpEstParNames,function(x)paste0(x,1:U)))
fixedParNames_expanded <- paste0(fixedParNames,1)

reg_rw.sd <- lapply(regEstParNames_expanded,function(x)0.02)
names(reg_rw.sd) <- regEstParNames_expanded
ivp_rw.sd <- lapply(ivpEstParNames_expanded,function(x)expression(ivp(0.05)))
names(ivp_rw.sd) <- ivpEstParNames_expanded
measles_rw.sd <- do.call(rw.sd,c(reg_rw.sd,ivp_rw.sd))

all_units = seq_len(length(unit_names(sim)))
nblocks = Nblocks
block_list = split(all_units, sort(all_units %% nblocks))
block_list <- lapply(block_list, as.integer)

fixed_params <- coef(sim,fixedParNames_expanded)

set.seed(53125455)
  runif_design(
    lower=coef(sim,estParNames_expanded,transform=TRUE) -1,
    upper=coef(sim,estParNames_expanded,transform=TRUE) +1,
    nseq=Reps
  ) -> guesses_transformed

guesses <- as.data.frame(t(apply(guesses_transformed,1,function(p,po,fi) {
    partrans(object=sim,unlist(c(p,fi)),dir="fromEst")
  }, po=sim,fi=fixed_params)))

set.seed(0)
bake(file=paste0(out_dir,"global_e",experiment,".rds"),{
  registerDoRNG(1270401374)
  foreach(guess=iter(guesses,by="row"), .combine=rbind) %dopar% {
    library(pomp)
    library(tidyverse)
    ibpf(sim,
      params=guess,
      sharedParNames=sharedParNames,
      unitParNames=unitParNames,
      Nbpf=M,
      spat_regression=spat_regression[experiment],
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
        logLik(bpfilter(sim,params=guess,Np=lik_eval_J,
          block_list=block_list))
      ) %>% logmeanexp(se=TRUE) -> ll
    c(start_loglik=ll[1],start_loglik.se=ll[2])
    } -> start_loglik
    attr(start_loglik,"ncpu") <- getDoParWorkers()
    start_loglik
  }) -> start_loglik
  results <- cbind(results,start_loglik)
  print(start_loglik)
  pdf(file=paste0(out_dir,"lik_change.pdf"))
  plot(results[,"start_loglik"],results[,"loglik"],
    xlab="log likelihood at search start",
    ylab="log likelihood at search end (dotted line at the truth)")
  abline(h=truth_loglik,lty="dashed",col="blue")
  dev.off()
}

results <-  filter(results,is.finite(loglik))
t_global <- attr(results,"system.time")
ncpu_global <- attr(results,"ncpu")

print(t_global)
print(ncpu_global)

results %>%
    filter(is.finite(loglik)) %>%
    arrange(-loglik) %>%
    write_csv(file=paste0(out_dir,"global_e",experiment,".csv"))

print(results[,c("loglik","loglik.se",estParNames_expanded[1:5])])

estParNames_contracted <- c(sharedParNames,as.vector(sapply(unitParNames,function(x)paste0(x,1:U))))

theta_true <- contract_params(coef(sim),
  fixedParNames=fixedParNames,
  sharedParNames=sharedParNames,
  unitParNames=unitParNames,
  U=U
)[estParNames_contracted]

if(0){
theta_est <- contract_params(coef(i1),
  fixedParNames=fixedParNames,
  sharedParNames=sharedParNames,
  unitParNames=unitParNames,
  average=TRUE,
  U=U
)[estParNames_contracted]
}

results %>%
#  filter(loglik>max(loglik)-50) %>% ## show all results, at least initially
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

  plot(sim,log=T,type="l")
  ggsave(filename=paste0(out_dir,"sim.pdf"))

  cat(capture.output(sessionInfo()),
    file=paste0(out_dir,"sessionInfo.txt"),sep="\n")

  # keep a copy of the sbat file if it exists
  if(file.exists(sbat_file)) file.copy(sbat_file,
    paste0(out_dir,sbat_file),overwrite=TRUE)		   

}
