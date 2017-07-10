# 多台发动机，MB-HMM建立模型

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
train_ntimes <- diff(c(0, which(diff(train_data$machnumber) != 0), length(train_data$machnumber)))
train_data$O7 <- 1-train_data$O7
train_data$O12 <- 1-train_data$O12
train_data$O20 <- 1-train_data$O20
train_data$O21 <- 1-train_data$O21

test_data <- subset(dd, X1 == randomseq[length(ntimes)], select = c(X1, X2, X7, X8, X9, X12, X13, X16, X17, X20, X22, X25, X26))
names(test_data) <- c("machnumber", "timeseq", "O2", "O3", "O4", "O7", "O8", "O11", "O12", "O15", "O17", "O20", "O21")
test_ntimes <- diff(c(0, which(diff(test_data$machnumber) != 0), length(test_data$machnumber)))
test_data$O7 <- 1-test_data$O7
test_data$O12 <- 1-test_data$O12
test_data$O20 <- 1-test_data$O20
test_data$O21 <- 1-test_data$O21

getMBmodel <- function(train_data, nstate = 4){
  train_ntimes <- diff(c(0, which(diff(train_data$machnumber) != 0), length(train_data$machnumber)))
  trainseq <- unique(train_data$machnumber)
  trainnum <- length(trainseq)
  nobs <- length(train_data[1,])-2
  ininit <- rep(c(1,rep(0, nstate-1)),trainnum)
  trinit <- matrix(0, ncol = nstate*trainnum, nrow = nstate*trainnum)
  obinit <- matrix(0, ncol = 2*nobs, nrow = nstate*trainnum)
  
  ininiti <- c(1,rep(0,nstate-1))
  triniti <- matrix(0, ncol = nstate, nrow = nstate)
  for(i in 1:nstate-1){
    triniti[i,i] <- 0.8
    triniti[i, i+1] <- 0.2
  }
  triniti[nstate, nstate] <- 1
  obiniti <- matrix(10, ncol = 2*nobs, nrow = nstate)
  for(i in 1:trainnum){
    datai <- subset(train_data, machnumber == trainseq[i])
    # datai <- train_data[which(train_data$machnumber == trainseq[i]),]
    # li <- floor(train_ntimes[i]*c(0, 0.4, 0.7, 0.9, 1.0))
    
    # for(k in 1:nobs){
    #   if(which.min(datai[,k+2])<which.max(datai[,k+2])){
    #     obiniti[,2*k-1] <- seq(min(datai[,k+2]), max(datai[,k+2]), length.out = 4)
    #   }
    #   else {
    #     obiniti[,2*k-1] <- seq(max(datai[,k+2]), min(datai[,k+2]), length.out = 4)
    #   }
    # }
    baseobinit <- c(0.10, 0.30, 0.50, 0.70)
    obiniti[,1] <- baseobinit
    obiniti[,3] <- baseobinit
    obiniti[,5] <- baseobinit
    obiniti[,7] <- baseobinit
    obiniti[,9] <- baseobinit
    obiniti[,11] <- baseobinit
    obiniti[,13] <- baseobinit
    obiniti[,15] <- baseobinit
    obiniti[,17] <- baseobinit
    obiniti[,19] <- baseobinit
    obiniti[,21] <- baseobinit
# 	  for(k in 1:nobs){
#       for(j in 1:nstate){
#         temp <- (li[j]+1):(li[j+1])
#       	obiniti[j, 2*k-1] <- mean(datai[temp,k+2])
#       	obiniti[j, 2*k] <- 2*sqrt(var(datai[temp, k+2]))
#       	# obiniti[j, 2*k] <- 1
# 	    }
# 	  }
  	modi <- depmix(list(O2 ~ 1, O3 ~ 1, O4 ~ 1, O7 ~ 1, O8 ~ 1, O11 ~ 1, O12 ~ 1, O15 ~ 1, O17 ~ 1, O20 ~ 1, O21 ~ 1), 
                   data = datai, 
                   nstate = nstate, 
                   family = list(gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian()), 
                   instart = ininiti, 
                   trstart = t(triniti), 
                   respstart = t(obiniti),
                   ntimes = train_ntimes[i])
	  fmi <- fit(modi)
	  parmi <- c(unlist(getpars(fmi)))
  	Ai <- matrix(parmi[(nstate+1):(nstate*(nstate+1)) ], ncol = nstate, byrow = T)
    Bi <- matrix(parmi[(nstate*(nstate+1)+1):npar(fmi)], nrow = nstate, byrow = T)
	  trinit[((i-1)*nstate+1):(i*nstate),((i-1)*nstate+1):(i*nstate)] <- Ai
	  obinit[((i-1)*nstate+1):(i*nstate),] <- Bi
  }
  c(ininit, t(trinit), t(obinit))
}
tteemmpp <- getMBmodel(train_data)


# 之后是从其它地方复制粘贴的，无实际用处


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
# nBi <- matrix(nparm[(nstate*(nstate+1)+1):npar(fm)], nrow = nstate, byrow = T)
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



baseobinit <- c(0.25, 0.45, 0.60, 0.90)

obiniti[,1] <- baseobinit

obiniti[,3] <- baseobinit

obiniti[,5] <- baseobinit

obiniti[,7] <- 1-baseobinit

obiniti[,9] <- baseobinit

obiniti[,11] <- baseobinit

obiniti[,13] <- 1-baseobinit

obiniti[,15] <- baseobinit

obiniti[,17] <- baseobinit

obiniti[,19] <- 1-baseobinit

obiniti[,21] <- 1-baseobinit
