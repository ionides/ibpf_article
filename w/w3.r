# cross-checking he10.R against he10 logliks
# using bpfilter

# CHECKS OKAY!

code_version <- "w3"
out_dir <- paste0(code_version,"/")
if(!dir.exists(out_dir)) dir.create(out_dir)

library(spatPomp)
library(tidyverse)
source("../ibpf.R")
source("../param_formats.R")
load("../data/twentycities.rda")
source("../he10.R")

parNames <- c("S_0","E_0","I_0","alpha","R0","psi","g","sigma","gamma","amplitude","cohort","sigmaSE","rho","mu","iota")
#for(u in 1:20){
  town_vec <- c(1:20) # note: town_vec should be an increasing sequence
  U <- length(town_vec)
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
    dt=1/365,
    town_vec=town_vec
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
  
  bpf1 <- bake(file=paste0(out_dir,"bpf1.rds"),{
    bpfilter(spo,Np=10000,block_size=1)
  })

  print(bpf1@unit_names)
  print(logLik(bpf1))
  print(sum(he10mle$loglik[he10mle$town%in%bpf1@unit_names]))

# [1] -40350.79
# [1] -40345.7


#}

data.frame(town=bpf1@unit_names,block_loglik=apply(bpf1@block.cond.loglik,1,sum),
he2010=he10mle$loglik[match(bpf1@unit_names,he10mle$town)])

#                town block_loglik  he2010
#1             London   -3803.4534 -3804.9
#2         Birmingham   -3248.5217 -3239.3
#3          Liverpool   -3406.0646 -3403.1
#4         Manchester   -3252.2560 -3250.9
#5          Sheffield   -2808.6594 -2810.7
#6              Leeds   -2916.6777 -2918.6
#7            Bristol   -2681.2205 -2681.6
#8         Nottingham   -2704.4924 -2703.5
#9           Bradford   -2585.7776 -2586.6
#10              Hull   -2728.7897 -2729.4
#11           Cardiff   -2364.0378 -2364.9
#12          Hastings   -1584.2113 -1583.7
#13           Consett   -1362.4865 -1362.9
#14         Bedwellty   -1125.0600 -1125.1
#15         Northwich   -1192.4434 -1195.1
#16          Oswestry    -696.6592  -696.1
#17 Dalton.in.Furness    -724.7791  -726.1
#18              Mold    -296.6964  -296.5
#19              Lees    -549.4389  -548.1
#20        Halesworth    -319.0659  -318.6
