# 多台发动机，MB-HMM建立模型

library(ggplot2)
library(depmixS4)

fname <- "train_FD0015.csv"
setwd("D:/MachineLearning/HMM/new")

dd <- read.csv(fname)
ntimes <- diff(c(0, which(diff(dd$X1) == 1), length(dd$X1)))

train_data <- read.table("train_data.txt", header = TRUE)
test_data <- read.table("test_data.txt", header = TRUE)
train_data_divide <- read.table("train_data_divide.txt")

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
  obiniti <- matrix(10, ncol = 2*nobs, nrow = nstate)
  for(i in 1:trainnum){
    datai <- subset(train_data, machnumber == trainseq[i] & timeseq >8 )
	li <- c(0,as.integer(unlist(train_data_divide[i,1:4]) )-8 )
	
	for(k in 1:nobs){
      for(j in 1:nstate){
        temp <- (li[j]+1):(li[j+1])
        triniti[j,j] <- 1-1.0/length(temp)
        if(j<nstate){
          triniti[j,j+1] <- 1.0/length(temp)
        }
      	obiniti[j, 2*k-1] <- mean(datai[temp,k+2])
      	obiniti[j, 2*k] <- 2*sqrt(var(datai[temp, k+2]))
		    # obiniti[j, 2*k] <- 1
	  }
	}
  	modi <- depmix(list(O2 ~ 1, O3 ~ 1, O4 ~ 1, O7 ~ 1, O8 ~ 1, O11 ~ 1, O12 ~ 1, O15 ~ 1, O17 ~ 1, O20 ~ 1, O21 ~ 1), 
                   data = datai, 
                   nstate = nstate, 
                   family = list(gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian(), gaussian()), 
                   instart = ininiti, 
                   trstart = t(triniti), 
                   respstart = t(obiniti),
                   ntimes = nrow(datai))
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
