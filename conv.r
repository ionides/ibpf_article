
#file <- "v18_3/global_e5.rds"
#file <- "v16_3/global_e3.rds"
#file <- "v14_3/global_e3.rds"
file <- "v15_3/global_e5.rds"

conv <- function(file,sim_dat=TRUE,U=20){
  x <- readRDS(file)
  plotPars <- c("alpha","psi","R0","gamma","sigma","sigmaSE",
    "cohort","amplitude","rho","g","S_0", "E_0", "I_0")
    
  if(sim_dat){
    library(pomp)
    readRDS("sim/sim.rds") -> sim
    params <- coef(sim)[paste0(plotPars,1)]
    names(params) <- plotPars
  }

  par(mfrow=c(length(plotPars),1))
  par(mai=c(0.05,0.8,0.01,0.1))
  par(omi=c(0.3,0,0.1,0))
  for(i in seq_along(plotPars)){
    p <- plotPars[i]
    unit_names <- rep("",U)
    if(i==length(plot)) unit_names <- 1:U
    boxplot(x[,paste0(p,1:U)],names=unit_names,
      xaxt=ifelse(i==length(plotPars),"t","n") )
    mtext(p,side=2,line=2.6)
    if(sim_dat) abline(h=params[p],col="red",lty="dashed")
  }
  mtext(1:U,line=0.7,at=seq(from=0.125,to=0.94,length=20),side=1,outer=T)
}

conv("v18_3/global_e5.rds")
#file <- "v16_3/global_e3.rds"
#file <- "v14_3/global_e3.rds"
#conv("v15_3/global_e5.rds")

