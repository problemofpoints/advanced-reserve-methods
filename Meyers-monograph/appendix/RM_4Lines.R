#
# Script to get Risk Margin for CSR model with four independent lines
# This script calculates the cost of capital risk margins for
#  All four lines combined
#  The marginal risk margin for each line
#  The standalone risk margin for each line - the simulation is done
#    through a copula so the results will be slightly different that the other
#    standalone CSR risk margin
#  Each insurer takes about 35 minutes to run in my computer
#
# by Glenn Meyers
#
rm(list = ls())      # clear workspace")
t1=Sys.time()
#
# R packages
#
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(data.table)
library(parallel)
library(doParallel)
library(copula)
library(ChainLadder)
#
# function to get the loss data
#
ins.line.data=function(a,g.code){
  b=subset(a,a$GRCODE==g.code)
  name=b$GRNAME
  grpcode=b$GRCODE
  w=b$AccidentYear-min(b$AccidentYear)+1
  d=b$DevelopmentLag
  cal=w+d-1
  cum_incloss=b[,6]
  cum_pdloss=b[,7]
  bulk_loss=b[,8]
  dir_premium=b[,9]
  ced_premium=b[,10]
  net_premium=b[,11]
  single=b[,12]
  posted_reserve97=b[,13]
  cum_incloss=cum_incloss-bulk_loss #per CAS Loss Reserve Database
  # get incremental paid losses - assume data is sorted by ay and lag
  inc_pdloss=numeric(0)
  for (i in unique(w)){
    s=(w==i)
    pl=c(0,cum_pdloss[s])
    ndev=length(pl)-1
    il=rep(0,ndev)
    for (j in 1:ndev){
      il[j]=pl[j+1]-pl[j]
    }
    inc_pdloss=c(inc_pdloss,il)
  }
  # ad hoc adjustment to accomodate lognormal distribution
  cum_incloss=pmax(cum_incloss,1)
  cum_pdloss=pmax(cum_pdloss,1)
  #
  data.out=data.table(grpcode,w,d,cal,net_premium,dir_premium,ced_premium,
                      cum_pdloss,cum_incloss,bulk_loss,inc_pdloss,single,
                      posted_reserve97)
  data.out=data.out[order(d,w)]
  return(data.out)
}
#
# function to make simulations reproducible - from Zongwen Tan
#
extract_ordered <- function(stanfit) {
  library(dplyr)

  array <- rstan::extract(stanfit, permuted = FALSE)

  dim_fit <- dim(array)
  n_chains <- dim_fit[2]
  n_parameters <- dim_fit[3]

  mat_list <- lapply(1:n_chains, function(x) as.data.frame(array[, x, ]) )

  mat_df <- bind_rows(mat_list)

  cols_all <- colnames(mat_df)
  collapse_names <- gsub("[0-9]|\\[|\\]", "", cols_all) %>% unique()

  posterior_list <- lapply(collapse_names,
                           function(x) mat_df[cols_all[grepl(x, cols_all)]] %>% as.matrix)
  names(posterior_list) <- collapse_names
  return(posterior_list)
}
#
# Stan CSR script
#
CSRmodel_stan = "
data{
int<lower=1> len_data;
real logprem[len_data];
real logloss[len_data];
int<lower=1,upper=10> w[len_data];
int<lower=1,upper=10> d[len_data];
}
parameters{
real r_alpha[9];
real r_beta[9];
real logelr;
real <lower=0,upper=100000> a_ig[10];
real gamma;
}
transformed parameters{
real alpha[10];
real beta[10];
real speedup[10];
real sig2[10];
real sig[10];
real mu[len_data];
alpha[1] = 0;
for (i in 2:10) alpha[i] = r_alpha[i-1];
for (i in 1:9) beta[i] = r_beta[i];
beta[10] = 0;
speedup[1] = 1;
for (i in 2:10) speedup[i] = speedup[i-1]*(1-gamma);
sig2[10] = gamma_cdf(1/a_ig[10],1,1);
for (i in 1:9) sig2[10-i] = sig2[11-i]+gamma_cdf(1/a_ig[i],1,1);
for (i in 1:10) sig[i] = sqrt(sig2[i]);
for (i in 1:len_data){
mu[i] = logprem[i]+logelr+alpha[w[i]]+beta[d[i]]*speedup[w[i]];
}
}
model {
r_alpha ~ normal(0,3.162);
r_beta ~ normal(0,3.162);
for (i in 1:10) a_ig[i] ~ inv_gamma(1,1);
logelr ~ normal(-.4,3.162);
gamma ~ normal(0,0.05);
for (i in 1:len_data) logloss[i] ~ normal(mu[i],sig[d[i]]);
}
generated quantities{
vector[len_data] log_lik;
for (i in 1:len_data) log_lik[i] = normal_lpdf(logloss[i]|mu[i],sig[d[i]]);
}
"
#
# dummy data for compiling
#
data.dummy=list(len_data = 55,
                logprem  = rep(8,55),
                logloss  = rep(8,55),
                w        = c(1:10,1:9,1:8,1:7,1:6,1:5,1:4,1:3,1,2,1),
                d        = c(rep(1,10),rep(2,9),rep(3,8),rep(4,7),rep(5,6),
                             rep(6,5),rep(7,4),rep(8,3),rep(9,2),10))
#
# initialization function for CSRmodel_stan
#
init.CSR=function(chain_id){
  set.seed(12345)
  list(r_alpha=rnorm(9,0,0.2),r_beta=runif(9),a_ig=runif(10),
       logelr=runif(1,-0.75,-0.5),gamma=rnorm(1,0,0.1))
}
#
# compile the CSR model
#
fitC_CSR = stan(model_code=CSRmodel_stan,data=data.dummy,
                seed=12345,
                init=init.CSR,chains=0)
#
pars.list=c("alpha","beta","gamma","logelr","sig","log_lik")
#
#
# function to generate capital scenarios
# for one line model at ultimate time horizon
#
line_ultbe=function(CASdata,grpcode,fixed.rate){
  lossData <- ins.line.data(CASdata,grpcode)
  train_data=subset(lossData,lossData$cal<11)
  Premium=train_data$net_premium[1:10]
  logprem=log(Premium)
  loss=train_data$cum_pdloss
  #
  #  data for the model
  #
  data.CSR=list(len_data = length(loss),
                logprem  = log(train_data$net_premium),
                logloss  = log(loss),
                w        = train_data$w,
                d        = train_data$d)
  #
  # run the model
  #
  stan_thin=1
  stan_iter=5000
  Rhat_target=1.05
  max_Rhat=2
  while ((max_Rhat > Rhat_target)&(stan_thin<17)){
    fitCSR=stan(fit = fitC_CSR, data = data.CSR,init=init.CSR,
                seed = 12345,iter=stan_iter,thin=stan_thin,
                chains = 4,pars=pars.list,
                control=list(adapt_delta=.9999,max_treedepth=50),
                refresh=0)
    fitCSR_summary=
      as.matrix(summary(fitCSR)$summary)[1:32,c(1,3,10)]
    mrh=subset(fitCSR_summary,is.na(fitCSR_summary[,3])==F)
    max_Rhat=round(max(mrh[,3]),4)
    print(paste("Maximum Rhat =",max_Rhat,"Thin =",stan_thin))
    stan_thin=2*stan_thin
    stan_iter=2*stan_iter
  }
  stan_thin=stan_thin/2
  #
  # extract information from stan output to process in R
  #
  b <- extract_ordered(fitCSR)
  alpha=b$alpha
  beta=b$beta
  gamma=b$gamma
  logelr=b$logelr
  sig=b$sig
  num.mcmc=length(logelr)
  #
  # get speedup rates by accident year
  #
  speedup=matrix(1,num.mcmc,10)
  for (i in 2:10){
    speedup[,i]=speedup[,i-1]*(1-gamma)
  }
  #
  # get the best estimate
  #
  trpaid=as.triangle(lossData,origin="w",dev="d",value="cum_pdloss")
  cl <- makeCluster(4)  # replace unpaid (w,d) cells with expected paid value
  registerDoParallel(cl)
  bestest=foreach (i=1:num.mcmc,.combine=rbind) %dopar%{
    for (cy in 1:9){
      for (w in (cy+1):10){
        d=11+cy-w
        trpaid[w,d]=exp(logprem[w]+logelr[i]+alpha[i,w]+
                          beta[i,d]*speedup[i,w]+sig[i,d]^2/2)
      }
    }
    pv=0  # get the present value of payouts
    for (cy in 1:9){
      for (w in (cy+1):10){
        d=11+cy-w
        pv=pv+(trpaid[w,d]-trpaid[w,d-1])/(1+fixed.rate)^(cy-.5)
      }
    }
    best=pv
  }
  stopCluster(cl)
  best.estimate=round(mean(bestest))
  #
  # get log of simulated losses from each parameter set
  #
  mu_p=array(data=0,dim=c(num.mcmc,10,10))
  logloss_p=array(data=0,dim=c(num.mcmc,10,10))
  #
  # get the mus generated from the lower triangle

  for (d in 2:10){
    for (w in (12-d):10){
      mu_p[,w,d]=logprem[w]+logelr+alpha[,w]+beta[,d]
      logloss_p[,w,d]=rnorm(num.mcmc,mu_p[,w,d],sig[,d])
    }
  }
  #
  # unconditional ultimate loss estimates by scenario
  #
  mean.ult=matrix(lossData$cum_pdloss[91],num.mcmc,10)
  # = paid loss for w=1,d=10
  for (w in 2:10){
    mean.ult[,w]=Premium[w]*exp(logelr+alpha[,w]+sig[,10]^2/2)
  }
  ultall=rowSums(mean.ult)
  #
  # that's all for this function
  #
  ultbe=list(best.estimate = best.estimate,
             ultall        = ultall,
             logloss.lowtri= logloss_p,
             mu.lowtri     = mu_p,
             sig           = sig)
  return(ultbe)
}
#
# calculate the posterior post_TVaR for the needed lines combinations-
#
post_assets=function(cop,post1,ult1,post2,ult2,post3,ult3,post4,ult4,k){
  set.seed(k)
  x1=sample(ult1,10000,replace=T,post1)
  x2=sample(ult2,10000,replace=T,post2)
  x3=sample(ult3,10000,replace=T,post3)
  x4=sample(ult4,10000,replace=T,post4)
  x1=quantile(x1,cop[,1])
  x2=quantile(x2,cop[,2])
  x3=quantile(x3,cop[,3])
  x4=quantile(x4,cop[,4])
  pa=rep(0,18)
  pa[1]=mean(x2+x3+x4)
  pa[2]=mean(sort(x2+x3+x4)[TVaR.Range])
  pa[3]=mean(x1+x3+x4)
  pa[4]=mean(sort(x1+x3+x4)[TVaR.Range])
  pa[5]=mean(x1+x2+x4)
  pa[6]=mean(sort(x1+x2+x4)[TVaR.Range])
  pa[7]=mean(x1+x2+x3)
  pa[8]=mean(sort(x1+x2+x3)[TVaR.Range])
  pa[9]=mean(x1+x2+x3+x4)
  pa[10]=mean(sort(x1+x2+x3+x4)[TVaR.Range])
  pa[11]=mean(x1)
  pa[12]=mean(sort(x1)[TVaR.Range])
  pa[13]=mean(x2)
  pa[14]=mean(sort(x2)[TVaR.Range])
  pa[15]=mean(x3)
  pa[16]=mean(sort(x3)[TVaR.Range])
  pa[17]=mean(x4)
  pa[18]=mean(sort(x4)[TVaR.Range])
  return(pa)
}
#
# get the list of four-line insurers
#
CA.list=read.csv("~/Dropbox/CAS Loss Reserve Database/Selected Insurers CA.csv")
PA.list=read.csv("~/Dropbox/CAS Loss Reserve Database/Selected Insurers PA.csv")
WC.list=read.csv("~/Dropbox/CAS Loss Reserve Database/Selected Insurers WC.csv")
OL.list=read.csv("~/Dropbox/CAS Loss Reserve Database/Selected Insurers OL.csv")
inall=intersect(CA.list$x,PA.list$x)
inall=intersect(inall,WC.list$x)
inall=intersect(inall,OL.list$x)
#
fixed.rate=0.04
risky.rate=0.10
rho=0.0
cop=rCopula(10000, normalCopula(rho, dim = 4))
TVaR.Range=9701:10000
#
outframe=NULL
outfile=paste("Fourline_Risk_Margin",round(rho,2),".csv",sep="_")
#
for (g in inall[3]){   # the [4} is only one of the seven available]
grpcode=g
#
CASdata =
  read.csv("http://www.casact.org/research/reserve_data/comauto_pos.csv")
ultbe=line_ultbe(CASdata,grpcode,fixed.rate)
ult1=ultbe$ultall
logloss.lowtri1=ultbe$logloss.lowtri
best.estimate1=ultbe$best.estimate
mu.lowtri1=ultbe$mu.lowtri
sig1=ultbe$sig
#
CASdata =
  read.csv("http://www.casact.org/research/reserve_data/ppauto_pos.csv")
ultbe=line_ultbe(CASdata,grpcode,fixed.rate)
ult2=ultbe$ultall
logloss.lowtri2=ultbe$logloss.lowtri
best.estimate2=ultbe$best.estimate
mu.lowtri2=ultbe$mu.lowtri
sig2=ultbe$sig
#
CASdata =
  read.csv("http://www.casact.org/research/reserve_data/wkcomp_pos.csv")
ultbe=line_ultbe(CASdata,grpcode,fixed.rate)
ult3=ultbe$ultall
logloss.lowtri3=ultbe$logloss.lowtri
best.estimate3=ultbe$best.estimate
mu.lowtri3=ultbe$mu.lowtri
sig3=ultbe$sig
#
CASdata =
  read.csv("http://www.casact.org/research/reserve_data/othliab_pos.csv")
ultbe=line_ultbe(CASdata,grpcode,fixed.rate)
ult4=ultbe$ultall
logloss.lowtri4=ultbe$logloss.lowtri
best.estimate4=ultbe$best.estimate
mu.lowtri4=ultbe$mu.lowtri
sig4=ultbe$sig
#
# likelihood function of the observed new cy
#
llike=function(x_p,mu_p,sig,cy,sz){
  ll=rep(0,sz)
  for (w in (1+cy):10){
    ll=ll+dnorm(x_p[w,11+cy-w],mu_p[,w,11+cy-w],sig[,11+cy-w],log=T)
  }
  return(ll)
}
#
# get conditional estimates
#
p.mean.all=rep(0,10)
p.assets.all=rep(0,10)
p.mean.1=rep(0,10)
p.assets.1=rep(0,10)
p.mean.2=rep(0,10)
p.assets.2=rep(0,10)
p.mean.3=rep(0,10)
p.assets.3=rep(0,10)
p.mean.4=rep(0,10)
p.assets.4=rep(0,10)
p.mean.u1=rep(0,10)
p.assets.u1=rep(0,10)
p.mean.u2=rep(0,10)
p.assets.u2=rep(0,10)
p.mean.u3=rep(0,10)
p.assets.u3=rep(0,10)
p.mean.u4=rep(0,10)
p.assets.u4=rep(0,10)
#

cl <- makePSOCKcluster(4)
registerDoParallel(cl)
pred.mean.assets=foreach (i=1:10000,.combine=rbind) %dopar%{
  x1=logloss.lowtri1[i,,]
  x2=logloss.lowtri2[i,,]
  x3=logloss.lowtri3[i,,]
  x4=logloss.lowtri4[i,,]
  #
  loglike.1=rep(0,10000)
  loglike.2=rep(0,10000)
  loglike.3=rep(0,10000)
  loglike.4=rep(0,10000)
  p0=rep(.0001,10000)
  call.pa=post_assets(cop,p0,ult1,p0,ult2,p0,ult3,p0,ult4,i)
      p.mean.1[1]=call.pa[1]
    p.assets.1[1]=call.pa[2]
      p.mean.2[1]=call.pa[3]
    p.assets.2[1]=call.pa[4]
      p.mean.3[1]=call.pa[5]
    p.assets.3[1]=call.pa[6]
      p.mean.4[1]=call.pa[7]
    p.assets.4[1]=call.pa[8]
    p.mean.all[1]=call.pa[9]
  p.assets.all[1]=call.pa[10]
     p.mean.u1[1]=call.pa[11]
   p.assets.u1[1]=call.pa[12]
     p.mean.u2[1]=call.pa[13]
   p.assets.u2[1]=call.pa[14]
     p.mean.u3[1]=call.pa[15]
   p.assets.u3[1]=call.pa[16]
     p.mean.u4[1]=call.pa[17]
   p.assets.u4[1]=call.pa[18]

  #
  loglike.1=loglike.1+llike(x1,mu.lowtri1,sig1,1,10000)
  loglike=loglike.1-max(loglike.1)
  postint=sum(exp(loglike))
  posterior.1=exp(loglike)/postint
  #
  loglike.2=loglike.2+llike(x2,mu.lowtri2,sig2,1,10000)
  loglike=loglike.2-max(loglike.2)
  postint=sum(exp(loglike))
  posterior.2=exp(loglike)/postint
  #
  loglike.3=loglike.3+llike(x3,mu.lowtri3,sig3,1,10000)
  loglike=loglike.3-max(loglike.3)
  postint=sum(exp(loglike))
  posterior.3=exp(loglike)/postint
  #
  loglike.4=loglike.4+llike(x4,mu.lowtri4,sig4,1,10000)
  loglike=loglike.4-max(loglike.4)
  postint=sum(exp(loglike))
  posterior.4=exp(loglike)/postint
  #
  call.pa=post_assets(cop,
                      posterior.1,ult1,
                      posterior.2,ult2,
                      posterior.3,ult3,
                      posterior.4,ult4,i)
      p.mean.1[2]=call.pa[1]
    p.assets.1[2]=call.pa[2]
      p.mean.2[2]=call.pa[3]
    p.assets.2[2]=call.pa[4]
      p.mean.3[2]=call.pa[5]
    p.assets.3[2]=call.pa[6]
      p.mean.4[2]=call.pa[7]
    p.assets.4[2]=call.pa[8]
    p.mean.all[2]=call.pa[9]
  p.assets.all[2]=call.pa[10]
     p.mean.u1[2]=call.pa[11]
   p.assets.u1[2]=call.pa[12]
     p.mean.u2[2]=call.pa[13]
   p.assets.u2[2]=call.pa[14]
     p.mean.u3[2]=call.pa[15]
   p.assets.u3[2]=call.pa[16]
     p.mean.u4[2]=call.pa[17]
   p.assets.u4[2]=call.pa[18]
  #
  #
   loglike.1=loglike.1+llike(x1,mu.lowtri1,sig1,2,10000)
   loglike=loglike.1-max(loglike.1)
   postint=sum(exp(loglike))
   posterior.1=exp(loglike)/postint
   #
   loglike.2=loglike.2+llike(x2,mu.lowtri2,sig2,2,10000)
   loglike=loglike.2-max(loglike.2)
   postint=sum(exp(loglike))
   posterior.2=exp(loglike)/postint
   #
   loglike.3=loglike.3+llike(x3,mu.lowtri3,sig3,2,10000)
   loglike=loglike.3-max(loglike.3)
   postint=sum(exp(loglike))
   posterior.3=exp(loglike)/postint
   #
   loglike.4=loglike.4+llike(x4,mu.lowtri4,sig4,2,10000)
   loglike=loglike.4-max(loglike.4)
   postint=sum(exp(loglike))
   posterior.4=exp(loglike)/postint
   #
  call.pa=post_assets(cop,
                      posterior.1,ult1,
                      posterior.2,ult2,
                      posterior.3,ult3,
                      posterior.4,ult4,i)
      p.mean.1[3]=call.pa[1]
    p.assets.1[3]=call.pa[2]
      p.mean.2[3]=call.pa[3]
    p.assets.2[3]=call.pa[4]
      p.mean.3[3]=call.pa[5]
    p.assets.3[3]=call.pa[6]
      p.mean.4[3]=call.pa[7]
    p.assets.4[3]=call.pa[8]
    p.mean.all[3]=call.pa[9]
  p.assets.all[3]=call.pa[10]
     p.mean.u1[3]=call.pa[11]
   p.assets.u1[3]=call.pa[12]
     p.mean.u2[3]=call.pa[13]
   p.assets.u2[3]=call.pa[14]
     p.mean.u3[3]=call.pa[15]
   p.assets.u3[3]=call.pa[16]
     p.mean.u4[3]=call.pa[17]
   p.assets.u4[3]=call.pa[18]
  #
  #
   loglike.1=loglike.1+llike(x1,mu.lowtri1,sig1,3,10000)
   loglike=loglike.1-max(loglike.1)
   postint=sum(exp(loglike))
   posterior.1=exp(loglike)/postint
   #
   loglike.2=loglike.2+llike(x2,mu.lowtri2,sig2,3,10000)
   loglike=loglike.2-max(loglike.2)
   postint=sum(exp(loglike))
   posterior.2=exp(loglike)/postint
   #
   loglike.3=loglike.3+llike(x3,mu.lowtri3,sig3,3,10000)
   loglike=loglike.3-max(loglike.3)
   postint=sum(exp(loglike))
   posterior.3=exp(loglike)/postint
   #
   loglike.4=loglike.4+llike(x4,mu.lowtri4,sig4,3,10000)
   loglike=loglike.4-max(loglike.4)
   postint=sum(exp(loglike))
   posterior.4=exp(loglike)/postint
   #
  call.pa=post_assets(cop,
                      posterior.1,ult1,
                      posterior.2,ult2,
                      posterior.3,ult3,
                      posterior.4,ult4,i)
      p.mean.1[4]=call.pa[1]
    p.assets.1[4]=call.pa[2]
      p.mean.2[4]=call.pa[3]
    p.assets.2[4]=call.pa[4]
      p.mean.3[4]=call.pa[5]
    p.assets.3[4]=call.pa[6]
      p.mean.4[4]=call.pa[7]
    p.assets.4[4]=call.pa[8]
    p.mean.all[4]=call.pa[9]
  p.assets.all[4]=call.pa[10]
     p.mean.u1[4]=call.pa[11]
   p.assets.u1[4]=call.pa[12]
     p.mean.u2[4]=call.pa[13]
   p.assets.u2[4]=call.pa[14]
     p.mean.u3[4]=call.pa[15]
   p.assets.u3[4]=call.pa[16]
     p.mean.u4[4]=call.pa[17]
   p.assets.u4[4]=call.pa[18]
  #
  #
  #
   loglike.1=loglike.1+llike(x1,mu.lowtri1,sig1,4,10000)
   loglike=loglike.1-max(loglike.1)
   postint=sum(exp(loglike))
   posterior.1=exp(loglike)/postint
   #
   loglike.2=loglike.2+llike(x2,mu.lowtri2,sig2,4,10000)
   loglike=loglike.2-max(loglike.2)
   postint=sum(exp(loglike))
   posterior.2=exp(loglike)/postint
   #
   loglike.3=loglike.3+llike(x3,mu.lowtri3,sig3,4,10000)
   loglike=loglike.3-max(loglike.3)
   postint=sum(exp(loglike))
   posterior.3=exp(loglike)/postint
   #
   loglike.4=loglike.4+llike(x4,mu.lowtri4,sig4,4,10000)
   loglike=loglike.4-max(loglike.4)
   postint=sum(exp(loglike))
   posterior.4=exp(loglike)/postint
   #
  call.pa=post_assets(cop,
                      posterior.1,ult1,
                      posterior.2,ult2,
                      posterior.3,ult3,
                      posterior.4,ult4,i)
      p.mean.1[5]=call.pa[1]
    p.assets.1[5]=call.pa[2]
      p.mean.2[5]=call.pa[3]
    p.assets.2[5]=call.pa[4]
      p.mean.3[5]=call.pa[5]
    p.assets.3[5]=call.pa[6]
      p.mean.4[5]=call.pa[7]
    p.assets.4[5]=call.pa[8]
    p.mean.all[5]=call.pa[9]
  p.assets.all[5]=call.pa[10]
     p.mean.u1[5]=call.pa[11]
   p.assets.u1[5]=call.pa[12]
     p.mean.u2[5]=call.pa[13]
   p.assets.u2[5]=call.pa[14]
     p.mean.u3[5]=call.pa[15]
   p.assets.u3[5]=call.pa[16]
     p.mean.u4[5]=call.pa[17]
   p.assets.u4[5]=call.pa[18]
  #
  #
   loglike.1=loglike.1+llike(x1,mu.lowtri1,sig1,5,10000)
   loglike=loglike.1-max(loglike.1)
   postint=sum(exp(loglike))
   posterior.1=exp(loglike)/postint
   #
   loglike.2=loglike.2+llike(x2,mu.lowtri2,sig2,5,10000)
   loglike=loglike.2-max(loglike.2)
   postint=sum(exp(loglike))
   posterior.2=exp(loglike)/postint
   #
   loglike.3=loglike.3+llike(x3,mu.lowtri3,sig3,5,10000)
   loglike=loglike.3-max(loglike.3)
   postint=sum(exp(loglike))
   posterior.3=exp(loglike)/postint
   #
   loglike.4=loglike.4+llike(x4,mu.lowtri4,sig4,5,10000)
   loglike=loglike.4-max(loglike.4)
   postint=sum(exp(loglike))
   posterior.4=exp(loglike)/postint
   #
   call.pa=post_assets(cop,
                      posterior.1,ult1,
                      posterior.2,ult2,
                      posterior.3,ult3,
                      posterior.4,ult4,i)
      p.mean.1[6]=call.pa[1]
    p.assets.1[6]=call.pa[2]
      p.mean.2[6]=call.pa[3]
    p.assets.2[6]=call.pa[4]
      p.mean.3[6]=call.pa[5]
    p.assets.3[6]=call.pa[6]
      p.mean.4[6]=call.pa[7]
    p.assets.4[6]=call.pa[8]
    p.mean.all[6]=call.pa[9]
  p.assets.all[6]=call.pa[10]
     p.mean.u1[6]=call.pa[11]
   p.assets.u1[6]=call.pa[12]
     p.mean.u2[6]=call.pa[13]
   p.assets.u2[6]=call.pa[14]
     p.mean.u3[6]=call.pa[15]
   p.assets.u3[6]=call.pa[16]
     p.mean.u4[6]=call.pa[17]
   p.assets.u4[6]=call.pa[18]
  #
  #
   loglike.1=loglike.1+llike(x1,mu.lowtri1,sig1,6,10000)
   loglike=loglike.1-max(loglike.1)
   postint=sum(exp(loglike))
   posterior.1=exp(loglike)/postint
   #
   loglike.2=loglike.2+llike(x2,mu.lowtri2,sig2,6,10000)
   loglike=loglike.2-max(loglike.2)
   postint=sum(exp(loglike))
   posterior.2=exp(loglike)/postint
   #
   loglike.3=loglike.3+llike(x3,mu.lowtri3,sig3,6,10000)
   loglike=loglike.3-max(loglike.3)
   postint=sum(exp(loglike))
   posterior.3=exp(loglike)/postint
   #
   loglike.4=loglike.4+llike(x4,mu.lowtri4,sig4,6,10000)
   loglike=loglike.4-max(loglike.4)
   postint=sum(exp(loglike))
   posterior.4=exp(loglike)/postint
   #
  call.pa=post_assets(cop,
                      posterior.1,ult1,
                      posterior.2,ult2,
                      posterior.3,ult3,
                      posterior.4,ult4,i)
      p.mean.1[7]=call.pa[1]
    p.assets.1[7]=call.pa[2]
      p.mean.2[7]=call.pa[3]
    p.assets.2[7]=call.pa[4]
      p.mean.3[7]=call.pa[5]
    p.assets.3[7]=call.pa[6]
      p.mean.4[7]=call.pa[7]
    p.assets.4[7]=call.pa[8]
    p.mean.all[7]=call.pa[9]
  p.assets.all[7]=call.pa[10]
     p.mean.u1[7]=call.pa[11]
   p.assets.u1[7]=call.pa[12]
     p.mean.u2[7]=call.pa[13]
   p.assets.u2[7]=call.pa[14]
     p.mean.u3[7]=call.pa[15]
   p.assets.u3[7]=call.pa[16]
     p.mean.u4[7]=call.pa[17]
   p.assets.u4[7]=call.pa[18]
  #
  #
   loglike.1=loglike.1+llike(x1,mu.lowtri1,sig1,7,10000)
   loglike=loglike.1-max(loglike.1)
   postint=sum(exp(loglike))
   posterior.1=exp(loglike)/postint
   #
   loglike.2=loglike.2+llike(x2,mu.lowtri2,sig2,7,10000)
   loglike=loglike.2-max(loglike.2)
   postint=sum(exp(loglike))
   posterior.2=exp(loglike)/postint
   #
   loglike.3=loglike.3+llike(x3,mu.lowtri3,sig3,7,10000)
   loglike=loglike.3-max(loglike.3)
   postint=sum(exp(loglike))
   posterior.3=exp(loglike)/postint
   #
   loglike.4=loglike.4+llike(x4,mu.lowtri4,sig4,7,10000)
   loglike=loglike.4-max(loglike.4)
   postint=sum(exp(loglike))
   posterior.4=exp(loglike)/postint
   #
  call.pa=post_assets(cop,
                      posterior.1,ult1,
                      posterior.2,ult2,
                      posterior.3,ult3,
                      posterior.4,ult4,i)
      p.mean.1[8]=call.pa[1]
    p.assets.1[8]=call.pa[2]
      p.mean.2[8]=call.pa[3]
    p.assets.2[8]=call.pa[4]
      p.mean.3[8]=call.pa[5]
    p.assets.3[8]=call.pa[6]
      p.mean.4[8]=call.pa[7]
    p.assets.4[8]=call.pa[8]
    p.mean.all[8]=call.pa[9]
  p.assets.all[8]=call.pa[10]
     p.mean.u1[8]=call.pa[11]
   p.assets.u1[8]=call.pa[12]
     p.mean.u2[8]=call.pa[13]
   p.assets.u2[8]=call.pa[14]
     p.mean.u3[8]=call.pa[15]
   p.assets.u3[8]=call.pa[16]
     p.mean.u4[8]=call.pa[17]
   p.assets.u4[8]=call.pa[18]
  #
  #
   loglike.1=loglike.1+llike(x1,mu.lowtri1,sig1,8,10000)
   loglike=loglike.1-max(loglike.1)
   postint=sum(exp(loglike))
   posterior.1=exp(loglike)/postint
   #
   loglike.2=loglike.2+llike(x2,mu.lowtri2,sig2,8,10000)
   loglike=loglike.2-max(loglike.2)
   postint=sum(exp(loglike))
   posterior.2=exp(loglike)/postint
   #
   loglike.3=loglike.3+llike(x3,mu.lowtri3,sig3,8,10000)
   loglike=loglike.3-max(loglike.3)
   postint=sum(exp(loglike))
   posterior.3=exp(loglike)/postint
   #
   loglike.4=loglike.4+llike(x4,mu.lowtri4,sig4,8,10000)
   loglike=loglike.4-max(loglike.4)
   postint=sum(exp(loglike))
   posterior.4=exp(loglike)/postint
   #
  call.pa=post_assets(cop,
                      posterior.1,ult1,
                      posterior.2,ult2,
                      posterior.3,ult3,
                      posterior.4,ult4,i)
      p.mean.1[9]=call.pa[1]
    p.assets.1[9]=call.pa[2]
      p.mean.2[9]=call.pa[3]
    p.assets.2[9]=call.pa[4]
      p.mean.3[9]=call.pa[5]
    p.assets.3[9]=call.pa[6]
      p.mean.4[9]=call.pa[7]
    p.assets.4[9]=call.pa[8]
    p.mean.all[9]=call.pa[9]
  p.assets.all[9]=call.pa[10]
     p.mean.u1[9]=call.pa[11]
   p.assets.u1[9]=call.pa[12]
     p.mean.u2[9]=call.pa[13]
   p.assets.u2[9]=call.pa[14]
     p.mean.u3[9]=call.pa[15]
   p.assets.u3[9]=call.pa[16]
     p.mean.u4[9]=call.pa[17]
   p.assets.u4[9]=call.pa[18]
  #
  #
   loglike.1=loglike.1+llike(x1,mu.lowtri1,sig1,9,10000)
   loglike=loglike.1-max(loglike.1)
   postint=sum(exp(loglike))
   posterior.1=exp(loglike)/postint
   #
   loglike.2=loglike.2+llike(x2,mu.lowtri2,sig2,9,10000)
   loglike=loglike.2-max(loglike.2)
   postint=sum(exp(loglike))
   posterior.2=exp(loglike)/postint
   #
   loglike.3=loglike.3+llike(x3,mu.lowtri3,sig3,9,10000)
   loglike=loglike.3-max(loglike.3)
   postint=sum(exp(loglike))
   posterior.3=exp(loglike)/postint
   #
   loglike.4=loglike.4+llike(x4,mu.lowtri4,sig4,9,10000)
   loglike=loglike.4-max(loglike.4)
   postint=sum(exp(loglike))
   posterior.4=exp(loglike)/postint
   #
  call.pa=post_assets(cop,
                      posterior.1,ult1,
                      posterior.2,ult2,
                      posterior.3,ult3,
                      posterior.4,ult4,i)
      p.mean.1[10]=call.pa[1]
    p.assets.1[10]=call.pa[2]
      p.mean.2[10]=call.pa[3]
    p.assets.2[10]=call.pa[4]
      p.mean.3[10]=call.pa[5]
    p.assets.3[10]=call.pa[6]
      p.mean.4[10]=call.pa[7]
    p.assets.4[10]=call.pa[8]
    p.mean.all[10]=call.pa[9]
  p.assets.all[10]=call.pa[10]
     p.mean.u1[10]=call.pa[11]
   p.assets.u1[10]=call.pa[12]
     p.mean.u2[10]=call.pa[13]
   p.assets.u2[10]=call.pa[14]
     p.mean.u3[10]=call.pa[15]
   p.assets.u3[10]=call.pa[16]
     p.mean.u4[10]=call.pa[17]
   p.assets.u4[10]=call.pa[18]
  #
  result=c(p.mean.1,p.assets.1,p.mean.2,p.assets.2,p.mean.3,p.assets.3,
           p.mean.4,p.assets.4,p.mean.all,p.assets.all,p.mean.u1,p.assets.u1,
           p.mean.u2,p.assets.u2,p.mean.u3,p.assets.u3,p.mean.u4,p.assets.u4)
}
stopCluster(cl)
#
# calculate marginal capital
#
pred.E.all=as.matrix(pred.mean.assets[,81:90])
pred.A.all=as.matrix(pred.mean.assets[,91:100])
pred.C.all=pred.A.all-pred.E.all
release=pred.C.all[,1:9]*(1+fixed.rate)-pred.C.all[,2:10]
risk.margin.all=pred.C.all[,1]
for (i in 1:9){
  risk.margin.all=risk.margin.all-release[,i]/(1+risky.rate)^i
}

#
pred.E.m1=as.matrix(pred.mean.assets[,1:10])
pred.A.m1=as.matrix(pred.mean.assets[,11:20])
pred.C.m1=pred.A.m1-pred.E.m1
release=pred.C.m1[,1:9]*(1+fixed.rate)-pred.C.m1[,2:10]
risk.margin.m1=pred.C.m1[,1]
for (i in 1:9){
  risk.margin.m1=risk.margin.m1-release[,i]/(1+risky.rate)^i
}
#
pred.E.m2=as.matrix(pred.mean.assets[,21:30])
pred.A.m2=as.matrix(pred.mean.assets[,31:40])
pred.C.m2=pred.A.m2-pred.E.m2
release=pred.C.m2[,1:9]*(1+fixed.rate)-pred.C.m2[,2:10]
risk.margin.m2=pred.C.m2[,1]
for (i in 1:9){
  risk.margin.m2=risk.margin.m2-release[,i]/(1+risky.rate)^i
}
#
pred.E.m3=as.matrix(pred.mean.assets[,41:50])
pred.A.m3=as.matrix(pred.mean.assets[,51:60])
pred.C.m3=pred.A.m3-pred.E.m3
release=pred.C.m3[,1:9]*(1+fixed.rate)-pred.C.m3[,2:10]
risk.margin.m3=pred.C.m3[,1]
for (i in 1:9){
  risk.margin.m3=risk.margin.m3-release[,i]/(1+risky.rate)^i
}
#
pred.E.m4=as.matrix(pred.mean.assets[,61:70])
pred.A.m4=as.matrix(pred.mean.assets[,71:80])
pred.C.m4=pred.A.m4-pred.E.m4
release=pred.C.m4[,1:9]*(1+fixed.rate)-pred.C.m4[,2:10]
risk.margin.m4=pred.C.m4[,1]
for (i in 1:9){
  risk.margin.m4=risk.margin.m4-release[,i]/(1+risky.rate)^i
}
#
pred.E.1=as.matrix(pred.mean.assets[,101:110])
pred.A.1=as.matrix(pred.mean.assets[,111:120])
pred.C.1=pred.A.1-pred.E.1
release=pred.C.1[,1:9]*(1+fixed.rate)-pred.C.1[,2:10]
risk.margin.1=pred.C.1[,1]
for (i in 1:9){
  risk.margin.1=risk.margin.1-release[,i]/(1+risky.rate)^i
}
#
pred.E.2=as.matrix(pred.mean.assets[,121:130])
pred.A.2=as.matrix(pred.mean.assets[,131:140])
pred.C.2=pred.A.2-pred.E.2
release=pred.C.2[,1:9]*(1+fixed.rate)-pred.C.2[,2:10]
risk.margin.2=pred.C.2[,1]
for (i in 1:9){
  risk.margin.2=risk.margin.2-release[,i]/(1+risky.rate)^i
}
#
pred.E.3=as.matrix(pred.mean.assets[,141:150])
pred.A.3=as.matrix(pred.mean.assets[,151:160])
pred.C.3=pred.A.3-pred.E.3
release=pred.C.3[,1:9]*(1+fixed.rate)-pred.C.3[,2:10]
risk.margin.3=pred.C.3[,1]
for (i in 1:9){
  risk.margin.3=risk.margin.3-release[,i]/(1+risky.rate)^i
}
#
pred.E.4=as.matrix(pred.mean.assets[,161:170])
pred.A.4=as.matrix(pred.mean.assets[,171:180])
pred.C.4=pred.A.4-pred.E.4
release=pred.C.4[,1:9]*(1+fixed.rate)-pred.C.4[,2:10]
risk.margin.4=pred.C.4[,1]
for (i in 1:9){
  risk.margin.4=risk.margin.4-release[,i]/(1+risky.rate)^i
}
#

mrm1=mean(risk.margin.all)-mean(risk.margin.m1)
mrm2=mean(risk.margin.all)-mean(risk.margin.m2)
mrm3=mean(risk.margin.all)-mean(risk.margin.m3)
mrm4=mean(risk.margin.all)-mean(risk.margin.m4)
rmall=mean(risk.margin.all)
rm1=mean(risk.margin.1)
rm2=mean(risk.margin.2)
rm3=mean(risk.margin.3)
rm4=mean(risk.margin.4)
mrm.sum=mrm1+mrm2+mrm3+mrm4
mrm1.alloc=rmall*mrm1/mrm.sum
mrm2.alloc=rmall*mrm2/mrm.sum
mrm3.alloc=rmall*mrm3/mrm.sum
mrm4.alloc=rmall*mrm4/mrm.sum
#
Ultimate_Loss=round(c(mean(ult1),mean(ult2),
                      mean(ult3),mean(ult4),
                      mean(ult1+ult2+ult3+ult4)))
Best_Estimate=c(best.estimate1,best.estimate2,
                best.estimate3,best.estimate4,
                best.estimate1+best.estimate2+
                best.estimate3+best.estimate4)
Marg_Risk_Margin=round(c(mrm1,mrm2,mrm3,mrm4,sum(mrm1,mrm2,mrm3,mrm4)))
Risk_Margin_Alloc=round(c(mrm1.alloc,mrm2.alloc,mrm3.alloc,mrm4.alloc,rmall))
Univ_Risk_Margin=round(c(rm1,rm2,rm3,rm4))
Univ_Risk_Margin=c(Univ_Risk_Margin,sum(Univ_Risk_Margin))
Group=rep(g,5)
outframe=rbind(outframe,data.frame(Group,Ultimate_Loss,Best_Estimate,Marg_Risk_Margin,
                    Risk_Margin_Alloc,Univ_Risk_Margin))
#write.csv(outframe,outfile)
print(outframe)
#
t2=Sys.time()
print(t2-t1)
}
