
# evaluates the likelihood for the measles data
# how many particles / replicates do we need?

# we use the mle from he10, so this is a test of the block particle filter
# in a situation where there is no coupling

# estimates time for run level 3: 8 hr on greatlakes

code_version <- "w2"
run_level <- 3
#save_session_data <- TRUE
save_session_data <- FALSE

print(paste0("run_level = ",run_level))

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

code_file <- paste0(code_version,".r")
out_dir <- paste0(code_version,"_",run_level,"/")
sbat_file <- paste0(code_version,"_",run_level,".sbat")

if(!dir.exists(out_dir)) dir.create(out_dir)

# standard pomp pseudocode notation: U = units, N = time points,
# J = particles, M = filter iterations.
# Also, Reps = replicates


parNames <- c("S_0","E_0","I_0","alpha","R0","psi","g","sigma","gamma","amplitude","cohort","sigmaSE","rho","mu","iota")
U <- 20
N <- 730
sharedParNames <- NULL
fixedParNames <- c("mu","g")
unitParNames <- setdiff(parNames,c(sharedParNames,fixedParNames))
estParNames <- c(sharedParNames,unitParNames)
ivpParNames <- c("S_0","E_0","I_0")
ivpEstParNames <- intersect(ivpParNames,estParNames)
regEstParNames <- setdiff(estParNames,ivpParNames)

Tmax <- min(1950 + N/52, 1964)

spo <- he10(U=U, Tmax=Tmax,
  sharedParNames=sharedParNames,
  fixedParNames=fixedParNames,
  unitParNames=unitParNames,
  simulated=FALSE,
  dt=1/365
)
he10mle <- read.csv("../data/he10mle.csv")
theta <- rep(NA,length(coef(spo)))
names(theta) <- names(coef(spo))

for(u in 1:U){
  he_u <- which(he10mle[,1]==spo@unit_names[u])
  theta[paste0(unitParNames,u)]<- as.numeric(he10mle[he_u,unitParNames])  
}

theta["g1"]<-0
theta["mu1"] <- 0.02

coef(spo) <- theta

# based on ionides21, we expect to use each city as a block
block_size <- 1

if(run_level==1){
  Jvec <- c(100,50)
  Reps <- 3
} else if(run_level==2){
  Jvec <- c(10000,5000,1000)
  Reps <- 4
} else if(run_level==3){
  Jvec <- c(25600,12800,6400,3200,1600)
  Reps <- 100
}

#do_bake <- FALSE
do_bake <- TRUE
# for debugging problems of bake recomputing global.rds

set.seed(12875)

if(do_bake){
bake(file=paste0(out_dir,"lik.rds"),{
  registerDoRNG(1270401374)
  foreach(J=iter(rep(Jvec,each=Reps)), .combine=rbind) %dopar% {
    library(pomp)
    library(tidyverse)
    tic <- Sys.time()
    L <- logLik(bpfilter(spo,Np=J,block_size=block_size))
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

}
