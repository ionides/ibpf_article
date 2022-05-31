
library(pomp)
loglik_he <- sum(read.csv("../data/he10mle.csv")$loglik)
# -40345.7

# global.rds is ordered by starting value enumeration
# global.csv is ordered by likelihood
LIK_ORDER <- TRUE 
if(LIK_ORDER){
  z1a <- read.csv(file="z1a_3/global.csv")
  z1b <- read.csv(file="z1b_3/global.csv")
  z1c <- read.csv(file="z1c_3/global.csv")
  z1d <- read.csv(file="z1d_3/global.csv")
  z2a <- read.csv(file="z2a_3/global.csv")
  z2b <- read.csv(file="z2b_3/global.csv")
  z2c <- read.csv(file="z2c_3/global.csv")
  z2d <- read.csv(file="z2d_3/global.csv")
} else {
  z1a <- readRDS(file="z1a_3/global.rds")
  z1b <- readRDS(file="z1b_3/global.rds")
  z1c <- readRDS(file="z1c_3/global.rds")
  z1d <- readRDS(file="z1d_3/global.rds")
  z2a <- readRDS(file="z2a_3/global.rds")
  z2b <- readRDS(file="z2b_3/global.rds")
  z2c <- readRDS(file="z2c_3/global.rds")
  z2d <- readRDS(file="z2d_3/global.rds")
}

pdf(file="z1a-vs-initial.pdf")
x <- z1a[,"start_loglik"]
y <- z1a$loglik
x_lim <- range(c(x,loglik_he,max(y)))
y_lim <- range(c(y,loglik_he,max(x)))
plot(x,y,xlim=x_lim,ylim=y_lim,
  xlab="starting log lik",
  ylab="z1 log lik"
)
abline(a=0,b=1,lty="dashed",col="blue")
abline(h=loglik_he,lty="dashed",col="blue")
abline(v=loglik_he,lty="dashed",col="blue")
dev.off()

pdf(file="z2a-vs-z1a.pdf")
x <- rep(z1a$loglik[1:18],each=2)
y <- z2a$loglik
x_lim <- range(c(x,loglik_he,max(y)))
y_lim <- range(c(y,loglik_he,max(x)))
plot(x,y,xlim=x_lim,ylim=y_lim,
  xlab="z1 log lik",
  ylab="z2 log lik"
)
abline(a=0,b=1,lty="dashed",col="blue")
abline(h=loglik_he,lty="dashed",col="blue")
abline(v=loglik_he,lty="dashed",col="blue")
dev.off()

pdf(file="z2b-vs-z1b.pdf")
x <- rep(z1b$loglik[1:18],each=2)
y <- z2b$loglik
x_lim <- range(c(x,loglik_he,max(y)))
y_lim <- range(c(y,loglik_he,max(x)))
plot(x,y,xlim=x_lim,ylim=y_lim,
  xlab="z1 log lik",
  ylab="z2 log lik"
)
abline(a=0,b=1,lty="dashed",col="blue")
abline(h=loglik_he,lty="dashed",col="blue")
abline(v=loglik_he,lty="dashed",col="blue")
dev.off()


a <- cbind(
A1 = z1a$loglik,
A2 = z2a$loglik
)

b <- cbind(
B1 = z1b$loglik,
B2 = z2b$loglik
)

c <- cbind(
C1 = z1c$loglik,
C2 = z2c$loglik
)

d <- cbind(
D1 = z1d$loglik,
D2 = z2d$loglik
)

loglik_max <- max(cbind(a,b,c,d),loglik_he)

pdf(file="z2-plot.pdf")
boxplot(cbind(a,b,c,d),ylim=loglik_max-c(5000,0))
abline(h=loglik_he,lty="dashed",col="blue")
dev.off()

