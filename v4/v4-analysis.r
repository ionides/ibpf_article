
library(pomp)
truth_eval <- readRDS("../w/w1_3/lik.rds")
truth_loglik<-logmeanexp(truth_eval[truth_eval[,"J"]==25600,"logLik"])

# global.rds is ordered by starting value enumeration
# global.csv is ordered by likelihood
LIK_ORDER <- TRUE 
if(LIK_ORDER){
  v1a <- read.csv(file="v1a_3/global.csv")
  v1b <- read.csv(file="v1b_3/global.csv")
  v1c <- read.csv(file="v1c_3/global.csv")
  v2a <- read.csv(file="v2a_3/global.csv")
  v2b <- read.csv(file="v2b_3/global.csv")
  v2c <- read.csv(file="v2c_3/global.csv")
  v3a <- read.csv(file="v3a_3/global.csv")
  v3b <- read.csv(file="v3b_3/global.csv")
  v3c <- read.csv(file="v3c_3/global.csv")
  v4a <- read.csv(file="v4a_3/global.csv")
  v4b <- read.csv(file="v4b_3/global.csv")
  v4c <- read.csv(file="v4c_3/global.csv")
} else {
  v1a <- readRDS(file="v1a_3/global.rds")
  v1b <- readRDS(file="v1b_3/global.rds")
  v1c <- readRDS(file="v1c_3/global.rds")
  v2a <- readRDS(file="v2a_3/global.rds")
  v2b <- readRDS(file="v2b_3/global.rds")
  v2c <- readRDS(file="v2c_3/global.rds")
}

# REMAKE=TRUE
REMAKE = FALSE
if(REMAKE){

  pdf(file="v1a-vs-initial.pdf")
  x <- v1a[,"start_loglik"]
  y <- v1a$loglik
  x_lim <- range(c(x,truth_loglik,max(y)))
  y_lim <- range(c(y,truth_loglik,max(x)))
  plot(x,y,xlim=x_lim,ylim=y_lim,
    xlab="starting log lik",
    ylab="v1 log lik"
  )
  abline(a=0,b=1,lty="dashed",col="blue")
  abline(h=truth_loglik,lty="dashed",col="blue")
  abline(v=truth_loglik,lty="dashed",col="blue")
  dev.off()

  pdf(file="v2a-vs-v1a.pdf")
  x <- rep(v1a$loglik[1:18],each=2)
  y <- v2a$loglik
  x_lim <- range(c(x,truth_loglik,max(y)))
  y_lim <- range(c(y,truth_loglik,max(x)))
  plot(x,y,xlim=x_lim,ylim=y_lim,
    xlab="v1 log lik",
    ylab="v2 log lik"
  )
  abline(a=0,b=1,lty="dashed",col="blue")
  abline(h=truth_loglik,lty="dashed",col="blue")
  abline(v=truth_loglik,lty="dashed",col="blue")
  dev.off()

  pdf(file="v2b-vs-v1b.pdf")
  x <- rep(v1b$loglik[1:18],each=2)
  y <- v2b$loglik
  x_lim <- range(c(x,truth_loglik,max(y)))
  y_lim <- range(c(y,truth_loglik,max(x)))
  plot(x,y,xlim=x_lim,ylim=y_lim,
    xlab="v1 log lik",
    ylab="v2 log lik"
  )
  abline(a=0,b=1,lty="dashed",col="blue")
  abline(h=truth_loglik,lty="dashed",col="blue")
  abline(v=truth_loglik,lty="dashed",col="blue")
  dev.off()

  pdf(file="v3a-vs-v2a.pdf")
  x <- rep(v2a$loglik[1:18],each=2)
  y <- v3a$loglik
  x_lim <- range(c(x,truth_loglik,max(y)))
  y_lim <- range(c(y,truth_loglik,max(x)))
  plot(x,y,xlim=x_lim,ylim=y_lim,
    xlab="v2a log lik",
    ylab="v3a log lik"
  )
  abline(a=0,b=1,lty="dashed",col="blue")
  abline(h=truth_loglik,lty="dashed",col="blue")
  abline(v=truth_loglik,lty="dashed",col="blue")
  dev.off()

  pdf(file="v3b-vs-v2b.pdf")
  x <- rep(v2b$loglik[1:18],each=2)
  y <- v3b$loglik
  x_lim <- range(c(x,truth_loglik,max(y)))
  y_lim <- range(c(y,truth_loglik,max(x)))
  plot(x,y,xlim=x_lim,ylim=y_lim,
    xlab="v2b log lik",
    ylab="v3b log lik"
  )
  abline(a=0,b=1,lty="dashed",col="blue")
  abline(h=truth_loglik,lty="dashed",col="blue")
  abline(v=truth_loglik,lty="dashed",col="blue")
  dev.off()

}

  pdf(file="v4a-vs-v3a.pdf")
  x <- rep(v3a$loglik[1:18],each=2)
  y <- v4a$loglik
  x_lim <- range(c(x,truth_loglik,max(y)))
  y_lim <- range(c(y,truth_loglik,max(x)))
  plot(x,y,xlim=x_lim,ylim=y_lim,
    xlab="v3a log lik",
    ylab="v4a log lik"
  )
  abline(a=0,b=1,lty="dashed",col="blue")
  abline(h=truth_loglik,lty="dashed",col="blue")
  abline(v=truth_loglik,lty="dashed",col="blue")
  dev.off()

  pdf(file="v4b-vs-v3b.pdf")
  x <- rep(v3b$loglik[1:18],each=2)
  y <- v4b$loglik
  x_lim <- range(c(x,truth_loglik,max(y)))
  y_lim <- range(c(y,truth_loglik,max(x)))
  plot(x,y,xlim=x_lim,ylim=y_lim,
    xlab="v3b log lik",
    ylab="v4b log lik"
  )
  abline(a=0,b=1,lty="dashed",col="blue")
  abline(h=truth_loglik,lty="dashed",col="blue")
  abline(v=truth_loglik,lty="dashed",col="blue")
  dev.off()

a <- cbind(
A1 = v1a$loglik,
A2 = v2a$loglik,
A3 = v3a$loglik,
A4 = v4a$loglik
)

b <- cbind(
B1 = v1b$loglik,
B2 = v2b$loglik,
B3 = v3b$loglik,
B4 = v4b$loglik
)

c <- cbind(
C1 = v1c$loglik,
C2 = v2c$loglik,
C3 = v3c$loglik,
C4 = v4c$loglik
)


loglik_max <- max(cbind(a,b,c),truth_loglik)

pdf(file="v4-plot.pdf")
boxplot(cbind(a,b,c),ylim=loglik_max-c(200,0))
abline(h=truth_loglik,lty="dashed",col="blue")
dev.off()

