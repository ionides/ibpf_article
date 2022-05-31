
# evaluating the likelihood at the truth on simulated data
# how many particles / replicates do we need?

# conclusions: 20 evaluations with 5000 particles
# for J=6400, we obtained one evaluation with considerably higher likelihood,
# leading to a large se. this suggests that there may be more subtlety to
# the evaluation.

# estimated time for run level 3: 8hr on greatlakes

#             400           800          1600          3200          6400
#   -39016.501733 -38945.237090 -38911.400026 -3.890015e+04 -38882.711908
#se      8.883076      1.915404      1.107329  9.755192e-01      9.926747

code_version <- "w1"
run_level <- 3
save_session_data <- TRUE

print(paste0("run_level = ",run_level))

library(spatPomp)
library(tidyverse)
source("../ibpf.R")
source("../param_formats.R")
source("../he10.R")

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

# standard pomp pseudocode notation: U = units, N = time points,
# J = particles, M = filter iterations.
# Also, Reps = replicates

# the number of estimated parameters should be irrelevant for likelihood
# evaluation. it does affect the extended model, but that is
# checked elsewhere to give rise to the same simulations, so
# here we use sim_test which has all fixed parameters


## parameters comparable to a typical city for he10, generally more like London
## selected and tested in sim.r
sim <- readRDS(file="../sim/sim.rds")

# here, we're interested in how to do this for the he10 data
U <- nrow(obs(sim)) # i.e., 20

# based on ionides21, we expect to use each city as a block
block_size <- 1

if(run_level==1){
  Jvec <- c(100,50)
  Reps <- 3
} else if(run_level==2){
  Jvec <- c(800,400,200,100)
  Reps <- 5
} else if(run_level==3){
  Jvec <- c(25600,12800,6400,3200,1600)
  Reps <- 100
}

#do_bake <- FALSE
do_bake <- TRUE
# for debugging problems of bake recomputing global.rds

set.seed(12345)

if(do_bake){
bake(file=paste0(out_dir,"lik.rds"),{
  registerDoRNG(1270401374)
  foreach(J=iter(rep(Jvec,each=Reps)), .combine=rbind) %dopar% {
    library(pomp)
    library(tidyverse)
    tic <- Sys.time()
    L <- logLik(bpfilter(sim,Np=J,block_size=block_size))
    toc <- Sys.time()
    c(J,toc-tic,L)
  } -> results
  colnames(results) <- c("J","time","logLik")
  attr(results,"ncpu") <- getDoParWorkers()
  results
}) -> results
}

if(!do_bake) results <- readRDS(file=paste0(out_dir,"lik.rds"))

t_global <- attr(results,"system.time")
ncpu_global <- attr(results,"ncpu")

print(t_global)
print(ncpu_global)

print(results)


pdf(file=paste0(out_dir,"hist.pdf"),width=8,height=10)

sp <- split(results[,"logLik"],results[,"J"])
lik_range <- range(results[,"logLik"])
par(mfrow=c(length(sp),1))
if(run_level==3) lower <- 2 else lower <- 1
for(i in lower:length(sp)){
  #hist(sp[[i]],xlim=xlim,main = paste("J =",names(sp)[i]))
  if(run_level<=2) hist(sp[[i]],main = paste("J =",names(sp)[i]))
  if(run_level==3) hist(sp[[i]],main = paste("J =",names(sp)[i]),
    xlim=c(lik_range[2]-100,lik_range[2]),
    breaks=c(lik_range[1],seq(from= lik_range[2]-100,to=lik_range[2],length=12))
    )
  
}

dev.off()

sapply(sp,function(x)logmeanexp(x,se=T))

if(save_session_data){

  cat(capture.output(sessionInfo()),
    file=paste0(out_dir,"sessionInfo.txt"),sep="\n")

  # keep a copy of the code used for the results
  file.copy(code_file,paste0(out_dir,code_file),overwrite=TRUE)		   

  # keep a copy of the sbat file if it exists
  if(file.exists(sbat_file)) file.copy(sbat_file,
    paste0(out_dir,sbat_file),overwrite=TRUE)		   

  # plot the data to check it is okay
  plot(sim,log=T,type="l")
  ggsave(filename=paste0(out_dir,"sim.pdf"))

}
