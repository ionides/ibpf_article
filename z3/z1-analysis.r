
library(pomp)
loglik_he <- sum(read.csv("../data/he10mle.csv")$loglik)
# -40345.7

a <- cbind(
A1 = readRDS(file="z1a_3/global.rds")$loglik
)

b <- cbind(
B1 = readRDS(file="z1b_3/global.rds")$loglik
)

c <- cbind(
C1 = readRDS(file="z1c_3/global.rds")$loglik
)

d <- cbind(
D1 = readRDS(file="z1d_3/global.rds")$loglik
)

loglik_max <- max(cbind(a,b,c,d),loglik_he)

pdf(file="z1-plot.pdf")
boxplot(cbind(a,b,c,d),ylim=loglik_max-c(5000,0))
abline(h=loglik_he,lty="dashed",col="blue")
dev.off()

