
# mle from metapop/code_and_data/plot/mif_global_nv1.csv
par1 <- c(
  alpha_be1=0.093898173,
  Beta_be1=0.779439748,
  alpha_af1=0.461829062,
  Beta_af1=0.134595027,
  mu_be1=0.966823556,
  Z_be1=0.116508318,
  D_be1=0.934673622,
  mu_af1=0.686401734,
  Z_af1=0.044448782,
  D_af1=0.126463381,
  theta1=0.89626861,
  tau1=0.917180577,
  size_a1=0.550023483/6, # adaptation for new code
  size_b1=0.056420797/6, # adaptation for new code
  Td_be1=0.8,
  Td_af1=0.2,
  E_01=0.44362243,
  Iu_01=0.984098384
)

parli1 <- c(
  alpha_be1=0.14,
  Beta_be1=1.12,
  alpha_af1=0.00001,
  Beta_af1=0.000001,
  mu_be1=0.55,
  Z_be1=3.69,
  D_be1=3.47,
  mu_af1=0.0000001,
  Z_af1=0.000001,
  D_af1=0.000001,
  theta1=1.36,
  tau1=0.0000001,
  size_a1=0.550023483/6, # adaptation for new code
  size_b1=0.056420797/6, # adaptation for new code
  Td_be1=0.8,
  Td_af1=0.2,
  E_01=0.5,
  Iu_01=0.5,
)

#run_level <- 1
#Np <-switch(run_level,5, 500, 1e3)
#mif <-switch(run_level, 3, 20, 40)
#Nbpf <-switch(run_level, 3, 20, 40)
# Nrep <-switch(run_level, 3,100,100)
# U <-switch(run_level, 5,50,375)

library(metapoppkg)
library(doParallel)
# plot covid model

library(doRNG)
library(patchwork)
library(ggplot2)

m1=li20(U=375,overdispersion_dynamic = T, overdispersion_measurement = T)

## early epidemic
if(0){
  e1 <- li20(U=100,overdispersion_dynamic = T, overdispersion_measurement = T)
  e1@times <- e1@times + 9 # to match days in January
  plot(window(e1,16,23),ty="h",log=T,plot_unit=F)
  ggsave(file="early.pdf",width=18,height=12,units="cm")
}  
  
if(0){
  e2 <- li20(U=100,overdispersion_dynamic = T, overdispersion_measurement = T)
  e2 <- window(e2,1,14)
  plot(e2,ty="h",log=T,plot_unit=F)
  ggsave(file="e2.pdf",width=18,height=12,units="cm")

#  set.seed(1)
  se2 <- simulate(e2,params=par1,seed=21)
  plot(se2,ty="h",log=T,plot_unit=F)
  ggsave(file="se2.pdf",width=18,height=12,units="cm")

for(i in 20:50){
  se <- simulate(e2,params=par1,seed=i)
  plot(se,ty="h",log=T,plot_unit=F)
  ggsave(file=paste0("se-",i,".pdf"),width=18,height=12,units="cm")
}




}  


coef(m1) <- par1

s1 <- simulate(m1)
p1 <- plot(s1,ty="h",log=T,plot_unit=F)

m2=li20(U=5,overdispersion_dynamic = T, overdispersion_measurement = T)

coef(m2) <- par1

s2 <- simulate(m2)
p2 <- plot(s2,ty="h",log=T)

# ## uses old code
# setwd("~/git/metapop/code_and_data/plot")
# rm(list=ls())
# library(doParallel)
# cores <- detectCores()
# registerDoParallel(cores)
# library(doRNG)
# library(ggplot2)
# library(reshape2)
# library(tidyverse)
# registerDoRNG(3123465)
# set.seed(12312)
# source("object_mifshare_nb.R")
# para=read.csv('mif_local_nb1.csv')[1,4:19]
# sim=simulate(covid,params=para)


