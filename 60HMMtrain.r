# 多台发动机，普通HMM建立模型

library(ggplot2)
library(depmixS4)

fname <- "train_FD0015.csv"
setwd("D:/MachineLearning/HMM/new")

dd <- read.csv(fname)
ntimes <- diff(c(0, which(diff(dd$X1) == 1), length(dd$X1)))


randomseq <- sample(1:length(ntimes), size = length(ntimes))
traincapacity <- floor(length(ntimes)*0.6)
train_data <- subset(dd, X1 %in% randomseq[1:traincapacity], select = c(X1, X2, X7, X8, X9, X12, X13, X16, X17, X20, X22, X25, X26))
names(train_data) <- c("machnumber", "timeseq", "O2", "O3", "O4", "O7", "O8", "O11", "O12", "O15", "O17", "O20", "O21")
# test_data <- subset(dd, (X1 %in% randomseq[1:traincapacity])==F, select = c(X1, X2, X7, X8, X9, X12, X13, X16, X17, X20, X22, X25, X26))
test_data <- subset(dd, X1 == randomseq[length(ntimes)], select = c(X1, X2, X7, X8, X9, X12, X13, X16, X17, X20, X22, X25, X26))
names(test_data) <- c("machnumber", "timeseq", "O2", "O3", "O4", "O7", "O8", "O11", "O12", "O15", "O17", "O20", "O21")

train_ntimes <- diff(c(0, which(diff(train_data$machnumber) != 0), length(train_data$machnumber)))
test_ntimes <- diff(c(0, which(diff(test_data$machnumber) != 0), length(test_data$machnumber)))

nstate <- 4
testmod <- depmix(list(O2 ~ 1, O3 ~ 1, O4 ~ 1, O7 ~ 1, O8 ~ 1, O11 ~ 1, O12 ~ 1, O15 ~ 1, O17 ~ 1, O20 ~ 1, O21 ~ 1), 
                  data = test_data, 
                  nstate = nstate, 
                  family = list(gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian()), 
                  ntimes = test_ntimes)


trinit <- matrix(0, ncol = nstate, nrow = nstate)
for(i in 1:nstate-1){
  trinit[i,i] <- 0.98
  trinit[i, i+1] <- 0.02
}
trinit[nstate, nstate] <- 1
nobs <- 11
normobinit <- matrix(0, ncol = 2*nobs, nrow = nstate)

startidx <- cumsum(c(0,train_ntimes))
for(i in 1:nobs){
  for(j in 1:4){
    temp <- c()
    for(k in 1:traincapacity){
	  l <- floor(train_ntimes[k]/4)
	  temp <- c(temp, (startidx[k]+(j-1)*l+1):(startidx[k]+j*l))
	}
	normobinit[j, 2*i-1] <- mean(train_data[temp,i+2])
	normobinit[j, 2*i] <- sqrt(var(train_data[temp, i+2]))
  }
}

normmod <- depmix(list(O2 ~ 1, O3 ~ 1, O4 ~ 1, O7 ~ 1, O8 ~ 1, O11 ~ 1, O12 ~ 1, O15 ~ 1, O17 ~ 1, O20 ~ 1, O21 ~ 1), 
                  data = train_data, 
                  nstate = nstate, 
                  family = list(gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian()), 
                  instart = c(1, rep(0, nstate-1)), 
                  trstart = t(trinit), 
                  respstart = t(normobinit),
                  ntimes = train_ntimes)
nfm <- fit(normmod)
nparm <- c(unlist(getpars(nfm)))
# nAi <- matrix(nparm[(nstate+1):(nstate*(nstate+1)) ], ncol = nstate, byrow = T)
nBi <- matrix(nparm[(nstate*(nstate+1)+1):npar(nfm)], nrow = nstate, byrow = T)
ntestmod <- setpars(testmod, nparm)
nfb <- forwardbackward(ntestmod)
statei <- c()
for(i in 1:test_ntimes){
  statei <- c(statei, which.max(nfb$alpha[i,]))
}
test_data$state <- factor(statei)
ub <- rep(0, 4)
for(i in 1:nstate){
  ub[i] <- max(which(test_data$state == i))
}
ggplot(data = test_data, mapping = aes(x = timeseq, y = O15, color = state)) + geom_line() + annotate("point", x = ub, y = test_data$O15[ub])
