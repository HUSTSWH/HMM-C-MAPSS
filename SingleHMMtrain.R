# build HMM model for C-MAPSS
# Error occurs while training
# Author: Song Wenhao
# E-mail: 867217582@qq.com

library(ggplot2)
library(depmixS4)

fname <- "train_FD0015.csv"
setwd("D:/MachineLearning/HMM/new")

dd <- read.csv(fname)

ntimes <- diff(c(0, which(diff(dd$X1) == 1), length(dd$X1)))
nstate <- 4
trinit <- matrix(0, ncol = nstate, nrow = nstate)
for(i in 1:nstate-1){
  trinit[i,i] <- 0.98
  trinit[i, i+1] <- 0.02
}
trinit[nstate, nstate] <- 1

# picked_data <- data.frame(machnumber <- dd$X1, timeseq <- dd$X2, O2 <- dd$X7, O3 <- dd$X8, O4 <- dd$X9, O7 <- dd$X12, O8 <- dd$X13, O11 <- dd$X16, O12 <- dd$X17, O15 <- dd$X20, O17 <- dd$X22, O20 <- dd$X25, O21 <- dd$X26)

nobs <- 11
obinit <- matrix(0.12, ncol = 2*nobs, nrow = nstate)
baseobinit <- c(0.25, 0.45, 0.60, 0.90)

obinit[,1] <- baseobinit

obinit[,3] <- baseobinit

obinit[,5] <- baseobinit

obinit[,7] <- 1-baseobinit

obinit[,9] <- baseobinit

obinit[,11] <- baseobinit

obinit[,13] <- 1-baseobinit

obinit[,15] <- baseobinit

obinit[,17] <- baseobinit

obinit[,19] <- 1-baseobinit

obinit[,21] <- 1-baseobinit

# picked_data1 <- picked_data[which(picked_data$machnumber %in% 1:99), ]

train_data <- subset(dd, X1 %in% c(1), select = c(X1, X2, X7, X8, X9, X12, X13, X16, X17, X20, X22, X25, X26))
names(train_data) <- c("machnumber", "timeseq", "O2", "O3", "O4", "O7", "O8", "O11", "O12", "O15", "O17", "O20", "O21")


mod <- depmix(list(O2 ~ 1, O3 ~ 1, O4 ~ 1, O7 ~ 1, O8 ~ 1, O11 ~ 1, O12 ~ 1, O15 ~ 1, O17 ~ 1, O20 ~ 1, O21 ~ 1), 
              data = train_data, 
			  nstate = nstate, 
			  family = list(gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian()), 
			  instart = c(1, rep(0, nstate-1)), 
			  trstart = t(trinit), 
			  respstart = t(obinit),
			  ntimes = ntimes[1])
			  
fm <- fit(mod)
parm <- c(unlist(getpars(fm)))
Ai <- matrix(parm[(nstate+1):(nstate*(nstate+1)) ], ncol = nstate, byrow = T)
Bi <- matrix(parm[(nstate*(nstate+1)+1):npar(fm)], nrow = nstate, byrow = T)

test_data <- subset(dd, X1 %in% c(2), select = c(X1, X2, X7, X8, X9, X12, X13, X16, X17, X20, X22, X25, X26))
names(test_data) <- c("machnumber", "timeseq", "O2", "O3", "O4", "O7", "O8", "O11", "O12", "O15", "O17", "O20", "O21")

testmod <- depmix(list(O2 ~ 1, O3 ~ 1, O4 ~ 1, O7 ~ 1, O8 ~ 1, O11 ~ 1, O12 ~ 1, O15 ~ 1, O17 ~ 1, O20 ~ 1, O21 ~ 1), 
              data = test_data, 
              nstate = nstate, 
              family = list(gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian()), 
              ntimes = ntimes[2])
testmod <- setpars(testmod, parm)
fb <- forwardbackward(testmod)

lb <- vector(length = nstate)
ub <- vector(length = nstate)
temp <- c()
for(i in 1:nstate){
  lb[i] <- min(which(fb$alpha[,i]>0.5))
  ub[i] <- max(which(fb$alpha[,i]>0.5))
  temp <- c(temp, rep(i,ub[i]-lb[i]+1))
}
test_data$state <- factor(temp)
ggplot(data = test_data, mapping = aes(x = timeseq, y = O3, color = state)) + geom_line() + annotate("point", x = ub, y = test_data$O3[ub])





normobinit <- matrix(0, ncol = 2*nobs, nrow = nstate)
h <- floor(length(train_data$timeseq)/4)
for(i in 1:4){
  for(j in 1:11){
    normobinit[i,2*j-1] <- mean(train_data[((i-1)*h+1):(i*h),j+2])
    normobinit[i, 2*j] <- sqrt(var(train_data[((i-1)*h+1):(i*h),j+2]))
  }
}
nmod <- depmix(list(O2 ~ 1, O3 ~ 1, O4 ~ 1, O7 ~ 1, O8 ~ 1, O11 ~ 1, O12 ~ 1, O15 ~ 1, O17 ~ 1, O20 ~ 1, O21 ~ 1), 
              data = train_data, 
              nstate = nstate, 
              family = list(gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian()), 
              instart = c(1, rep(0, nstate-1)), 
              trstart = t(trinit), 
              respstart = t(normobinit),
              ntimes = ntimes[1])
nfm <- fit(nmod)
nparm <- c(unlist(getpars(nfm)))
nAi <- matrix(nparm[(nstate+1):(nstate*(nstate+1)) ], ncol = nstate, byrow = T)
nBi <- matrix(nparm[(nstate*(nstate+1)+1):npar(fm)], nrow = nstate, byrow = T)
ntestmod <- setpars(testmod, nparm)
nfb <- forwardbackward(ntestmod)
