#==========================================================#
# R examples: structural modeling by using
#             overlapped group penalty for two data sets
# Copyright (c) 2018-2020 Chong Ma
#==========================================================#

library(smog)
require(coxed)

n=50;p=100
set.seed(2018)
# generate design matrix x
s=10
x=matrix(0,n,1+2*p)
x[,1]=sample(c(0,1),n,replace = T)
x[,seq(2,1+2*p,2)]=matrix(rnorm(n*p),n,p)
x[,seq(3,1+2*p,2)]=x[,seq(2,1+2*p,2)]*x[,1]

g=c(p+1,rep(1:p,rep(2,p)))  # groups
v=c(0,rep(1,2*p))           # penalization status
label=c("t",rep(c("prog","pred"),p))  # type of predictor variables

# generate beta
beta=c(rnorm(13,0,2),rep(0,ncol(x)-13))
beta[c(2,4,7,9)]=0

# generate y
data1=x%*%beta
noise1=rnorm(n)
snr1=as.numeric(sqrt(var(data1)/(s*var(noise1))))
y1=data1+snr1*noise1
lfit1=smog(x,y1,g,v,label,lambda1=10,lambda2=0,lambda3=10,family = "gaussian")
cvfit1=cv.smog(x,y1,g,v,labe,type = "GCV",family = "gaussian")

## generate binomial data
prob=exp(as.matrix(x)%*%as.matrix(beta))/(1+exp(as.matrix(x)%*%as.matrix(beta)))
y2=ifelse(prob>0.5,0,1)
lfit2=smog(x,y2,g,v,label,lambda1=0.05,lambda2=0,lambda3=0.05,family = "binomial")
cvfit2=cv.smog(x,y2,g,v,labe,type = "GCV",family = "binomial")

## generate survival data
data3=sim.survdata(N=n,T=100,X=x,beta=beta)
y3=data3$data[,c("y","failed")]
y3$failed=ifelse(y3$failed,1,0)
colnames(y3)=c("time","status")
lfit3=smog(x,y3,g,v,label,lambda1=0.2,lambda2=0,lambda3=0.2,family = "coxph")
cvfit3=cv.smog(x,y3,g,v,labe,type = "GCV",family = "coxph")















