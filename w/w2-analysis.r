
he10mle <- read.csv("../data/he10mle.csv")
loglik_he10 <- sum(he10mle$loglik)

results <- readRDS("w2_3/lik.rds")

library(pomp)
library(tidyverse)

logliks <- split(results[,"logLik"],results[,"J"]) %>% lapply(logmeanexp)

logmeanexp(results[results[,"J"]==25600,"logLik"])

unlist(c(he10=loglik_he10,logliks))
#      he10      1600      3200      6400     12800     25600 
# -40345.70 -40400.84 -40369.14 -40354.10 -40346.10 -40338.72 

