
# extension <- "_3"
# > q75
#        0      0.1      0.2      0.3     0.4       0.5     0.6       0.7 
# -38662.1 -35355.7 -35250.9 -35308.9 -35350.7 -35365.4 -35284.1 -35385.0 
#      0.8      0.9        1 
# -35252.8 -35383.1 -35388.2 

# boxplot.pdf shows that all values of spat_regression in [0.1,1] are
# more-or-less equally good. A linear model cannot pick up a linear
# or quadratic trend. However, spat_regression = 0 fails.

# interpretation: in the context of this test (a matter of approaching the
# mle, but not requiring to get within 100 log units) it does not matter
# the rate at which the units talk to each other concerning shared parameters.
# there may be a tradeoff - more diversity can lead to more information, but
# quicker collaboration can direct the search more quickly to good places.

# most runs are finding a lower-R0 regime (see R0 plot below for unit 5)
# though some are finding values close to the truth of 30.

loglik_max <- -35041 # from w1.r

library(foreach)

code_version <- "v6"

extension <- "_3"

out_dir <- paste0(code_version,extension,"/")
# spat_regression <- c(0.3,0.4,0.5,0.6,0.7)
spat_regression <- c(0.3,0.4,0.5,0.6,0.7,0,0.1,0.2,0.8,0.9,1)
# Experiment <- 5 # initial test
Experiment <- length(spat_regression)

results <- foreach(experiment=1:Experiment,.combine=rbind) %do% {
  r <- readRDS(file=paste0(out_dir,"global_",experiment,".rds"))
  cbind(r,sp_reg=spat_regression[experiment])
}

dim(results)

results <- as.data.frame(results)

result_list <- split(results$loglik,results$sp_reg)


density_list <- lapply(result_list,
  function(x) density(x,n=1000,bw=50,from=loglik_max-1000,to=loglik_max)
)
density_matrix <- sapply(density_list,function(x) x$y)


pdf(file=paste0(out_dir,"density-plot.pdf"))
matplot(x=density_list[[1]]$x, y=density_matrix,ty="l",
  lwd=as.numeric(names(result_list))*5,
  lty="solid",col="black",
  main = "loglik of 36 runs, spat_regression proportional to line width",
  ylab="kernel density",
  xlab="final log likelihood"
)
dev.off()

sp_reg <- as.numeric(names(result_list))

rho5 <- 0.5
rho5_list <- split(results$rho5,results$sp_reg)
rho5_density_list <- lapply(rho5_list,
  function(x) density(x,n=1000,from=0.25,to=0.7)
)
rho5_density_matrix <- sapply(rho5_density_list,function(x) x$y)
matplot(x=rho5_density_list[[1]]$x, y=rho5_density_matrix,ty="l",
  lwd=sp_reg*5,
  lty="solid",col=(1+as.numeric(sp_reg==0.6)),
  main = "rho5 for 36 runs, spat_regression proportional to line width",
  ylab="kernel density",
  xlab="final rho5 (0.6 is red)"
)

R05 <- 30
R05_list <- split(results$R05,results$sp_reg)
R05_density_list <- lapply(R05_list,
  function(x) density(x,n=1000,from=5,to=50)
)
R05_density_matrix <- sapply(R05_density_list,function(x) x$y)
matplot(x=R05_density_list[[1]]$x, y=R05_density_matrix,ty="l",
  lwd=sp_reg*5,
  lty="solid",col=(1+as.numeric(sp_reg==0.6)),
  main = "R05 for 36 runs, spat_regression proportional to line width",
  ylab="kernel density",
  xlab="final R05 (0.6 is red)"
)

q75 <- sapply(result_list,function(x) quantile(x,0.75))
sp_reg <- as.numeric(names(result_list))
names(q75) <- sp_reg
q75
q75lm <- lm(q75~sp_reg+I(sp_reg^2))
summary(q75lm)

plot(sp_reg,q75,ylim=c(min(q75),loglik_max))
abline(h=loglik_max)

pdf(file=paste0(out_dir,"boxplot.pdf"))
boxplot(loglik~sp_reg,data=results,ylim=c(-60000,loglik_max),
  ylab="log likelihood", xlab="spatial autoregression parameter")
abline(h=loglik_max,col="blue",lty="dotted")
dev.off()


# sp_reg=0 is much worse; no clear pattern in the remaining values.
plot(sp_reg[-1],q75[-1])

q75lm2 <- lm(q75[-1]~sp_reg[-1]+I(sp_reg[-1]^2))
summary(q75lm2)




