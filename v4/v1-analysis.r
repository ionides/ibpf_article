
library(pomp)
truth_eval <- readRDS("../w/w1_3/lik.rds")
truth_loglik<-logmeanexp(truth_eval[truth_eval[,"J"]==25600,"logLik"])

v1a <- readRDS(file="v1a_3/global.rds")
v1b <- readRDS(file="v1b_3/global.rds")
v1c <- readRDS(file="v1c_3/global.rds")

a <- cbind(
A1 = v1a$loglik
)

b <- cbind(
B1 = v1b$loglik
)

c <- cbind(
C1 = v1c$loglik
)


loglik_max <- max(cbind(a,b,c),truth_loglik)

pdf(file="v1-plot.pdf")
boxplot(cbind(a,b,c),ylim=loglik_max-c(200,0))
abline(h=truth_loglik,lty="dashed",col="blue")
dev.off()

