#
# Script to run the CSR risk margin model for the ultimate time horizon
# on data from the CAS Loss Reserve Database
# Uses Stan for the MCMC run
# by Glenn Meyers
# Run time on my computer - about 3 minutes
#
rm(list = ls())      # clear workspace")
t0=Sys.time()
#
# R packages
#
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(data.table)
library(parallel)
library(doParallel)
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
"
#
# initialization function for CSRModel
#
init.CSR=function(chain_id){
  set.seed(12345)
  list(r_alpha=rnorm(9,0,0.2),r_beta=runif(9),a_ig=runif(10),
       logelr=runif(1,-0.75,-0.5),gamma=rnorm(1,0,0.1))
}
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
# compile the CSR model
#
fitC_CSR = stan(model_code=CSRmodel_stan,data=data.dummy,
                seed=12345,
                init=init.CSR,chains=0)
#
pars.list=c("logelr","alpha","beta","gamma","sig")
#
# function to generate capital scenarios
# for one line model at ultimate time horizon
#
risk_margin_ult_horizon=function(CASdata,grpcode,fixed.rate,TVaR.Range){
  lossData <- ins.line.data(CASdata,grpcode)
  train_data=subset(lossData,lossData$cal<11)
  test_data=subset(lossData,lossData$cal>10)
  Premium=train_data$net_premium[1:10]
  logprem=log(Premium)
  loss=train_data$cum_pdloss
  test=test_data$cum_pdloss
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
  #
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
  # likelihood function of the observed new cy
  #
  llike=function(x_p,cy,sz){
    ll=rep(0,sz)
    for (w in (1+cy):10){
      ll=ll+dnorm(x_p[w,11+cy-w],mu_p[,w,11+cy-w],sig[,11+cy-w],log=T)
    }
    return(ll)
  }
  #
  # calculate the posterior mean and TVaR -
  #
  post_assets=function(post,ult){
    x=sample(ult,10000,replace=T,post)
    pa=rep(0,2)
    pa[1]=mean(x)
    pa[2]=mean(sort(x)[TVaR.Range])
    return(pa)
  }
  #
  # get  estimates by scenario of posterior required capital
  #
  p.mean=rep(0,10)
  p.assets=rep(0,10)
  cl <- makePSOCKcluster(4)
  registerDoParallel(cl)
  pred.mean.assets=foreach (i=1:num.mcmc,.combine=rbind) %dopar%{
    x_p=logloss_p[i,,]
    p.mean[1]=mean(ultall)
    p.assets[1]=mean(sort(ultall)[TVaR.Range])
    loglike=rep(0,num.mcmc)
    loglike=loglike+llike(x_p,1,num.mcmc)
    loglike2=loglike-max(loglike)
    postint=sum(exp(loglike2))
    posterior=exp(loglike2)/postint
    call.pa=post_assets(posterior,ultall)
    p.mean[2]=call.pa[1]
    p.assets[2]=call.pa[2]
    #
    loglike=loglike+llike(x_p,2,num.mcmc)
    loglike2=loglike-max(loglike)
    postint=sum(exp(loglike2))
    posterior=exp(loglike2)/postint
    call.pa=post_assets(posterior,ultall)
    p.mean[3]=call.pa[1]
    p.assets[3]=call.pa[2]
    #
    loglike=loglike+llike(x_p,3,num.mcmc)
    loglike2=loglike-max(loglike)
    postint=sum(exp(loglike2))
    posterior=exp(loglike2)/postint
    call.pa=post_assets(posterior,ultall)
    p.mean[4]=call.pa[1]
    p.assets[4]=call.pa[2]
    #
    loglike=loglike+llike(x_p,4,num.mcmc)
    loglike2=loglike-max(loglike)
    postint=sum(exp(loglike2))
    posterior=exp(loglike2)/postint
    call.pa=post_assets(posterior,ultall)
    p.mean[5]=call.pa[1]
    p.assets[5]=call.pa[2]
    #
    loglike=loglike+llike(x_p,5,num.mcmc)
    loglike2=loglike-max(loglike)
    postint=sum(exp(loglike2))
    posterior=exp(loglike2)/postint
    call.pa=post_assets(posterior,ultall)
    p.mean[6]=call.pa[1]
    p.assets[6]=call.pa[2]
    #
    loglike=loglike+llike(x_p,6,num.mcmc)
    loglike2=loglike-max(loglike)
    postint=sum(exp(loglike2))
    posterior=exp(loglike2)/postint
    call.pa=post_assets(posterior,ultall)
    p.mean[7]=call.pa[1]
    p.assets[7]=call.pa[2]
    #
    loglike=loglike+llike(x_p,7,num.mcmc)
    loglike2=loglike-max(loglike)
    postint=sum(exp(loglike2))
    posterior=exp(loglike2)/postint
    call.pa=post_assets(posterior,ultall)
    p.mean[8]=call.pa[1]
    p.assets[8]=call.pa[2]  #
    loglike=loglike+llike(x_p,8,num.mcmc)
    loglike2=loglike-max(loglike)
    postint=sum(exp(loglike2))
    posterior=exp(loglike2)/postint
    call.pa=post_assets(posterior,ultall)
    p.mean[9]=call.pa[1]
    p.assets[9]=call.pa[2]
    #
    loglike=loglike+llike(x_p,9,num.mcmc)
    loglike2=loglike-max(loglike)
    postint=sum(exp(loglike2))
    posterior=exp(loglike2)/postint
    call.pa=post_assets(posterior,ultall)
    p.mean[10]=call.pa[1]
    p.assets[10]=call.pa[2]
    result=c(p.mean,p.assets)
    }
  stopCluster(cl)
  #
  pred.E=as.matrix(pred.mean.assets[,1:10])
  pred.A=as.matrix(pred.mean.assets[,11:20])
  pred.C=pred.A-pred.E
  #
  output=list(num.mcmc      = num.mcmc,
              ultall        = round(mean(ultall)),
              best.estimate = round(best.estimate),
              pred.E        = pred.E,
              pred.A        = pred.A,
              pred.C        = pred.C)
  return(output)
}
#
# user inputs
#
fixed.rate=0.04
risky.rate=0.10
TVaR.Range=9701:10000
#
# single insurer with exhibits
CASdata = read.csv("~/Dropbox/CAS Loss Reserve Database/comauto_pos.csv")
# Location of files in the CAS Loss Reserve Database
#    http://www.casact.org/research/reserve_data/comauto_pos.csv
#    http://www.casact.org/research/reserve_data/ppauto_pos.csv
#    http://www.casact.org/research/reserve_data/wkcomp_pos.csv
#    http://www.casact.org/research/reserve_data/othliab_pos.csv
lineobus="CA"
grpcode=353
#
# calling output for single insurer for ultimate time horizon
# using the CSR standalone model
#
co_model=risk_margin_ult_horizon(CASdata,grpcode,fixed.rate,TVaR.Range)
num.mcmc=co_model$num.mcmc
ultall=co_model$ultall
best.estimate=co_model$best.estimate
pred.E=co_model$pred.E
pred.A=co_model$pred.A
pred.C=co_model$pred.C
#
# plot sample ultimate loss estimates for future calendar years
#
cyprobs=seq(0.01,.99,.02)
cyindex=trunc(quantile(1:num.mcmc,probs=cyprobs,names=F))
cynums=order(pred.E[,10])[cyindex]
#par(mfrow=c(1,1))
sub1=paste("E[Ultimate Loss] =",ultall,
           "- Best Estimate of Liability =",best.estimate)
plot(0:9,pred.E[cynums[1],],
     main="Paths of Ultimate Loss Estimates",
     xlab="Future Calendar Year",ylab="Ultimate Loss Estimate",
     ylim=range(pred.E[cynums,]),type="l",
     sub=sub1)
for(i in 2:length(cynums)){
  par(new=T)
  plot(0:9,pred.E[cynums[i],],main="",xlab="",ylab="",
       ylim=range(pred.E[cynums,]),type="l",sub="")
}
#
# plot sample capital requirements
#
plot(0:9,pred.C[cynums[1],],
     main="Required Capital by Calendar Year",
     xlab="Future Calendar Year",ylab="Required Capital",
     ylim=range(pred.C[cynums,]),type="l",
     sub=paste("Initial Capital =",round(pred.C[1,1])))
for(i in 2:length(cynums)){
  par(new=T)
  plot(0:9,pred.C[cynums[i],],main="",xlab="",ylab="",
       ylim=range(pred.C[cynums,]),type="l",sub="")
}
#
# plot some capital release paths
#
release=pred.C[,1:10]*(1+fixed.rate)-cbind(pred.C[,2:10],rep(0,num.mcmc))
duration=rep(0,num.mcmc)
par(mfrow=c(1,1))
plot(1:10,release[cynums[1],],main="Path of Released Capital",
     xlab="Future Calendar Year",ylab="Capital Released",
     ylim=range(release[cynums,]),type="l")
for(i in 2:length(cynums)){
  par(new=T)
  plot(1:10,release[cynums[i],],main="",xlab="",ylab="",
       ylim=range(release[cynums,]),type="l",sub="")
}
abline(0,0,col="red",lwd=3)
#
# calculate risk margins
#
risk.margin=pred.C[,1]
for (i in 1:9){
  risk.margin=risk.margin-release[,i]/(1+risky.rate)^i
}
par(mfrow=c(1,1))
hist(risk.margin,xlab="Risk Margin",main="Risk Margin",
     sub=paste("Mean Risk Margin =",round(mean(risk.margin))))
#
# plot risk margin as a % of the initial capital
#
risk.margin.pct=100*risk.margin/pred.C[1,1]
hist(risk.margin.pct,main="Risk Margin as a % of Initial Capital",
     xlab="Risk Margin %",breaks=25,
     sub=paste("Mean =",round(mean(risk.margin.pct),1),"%",
              "Std. Dev. =",round(sd(risk.margin.pct),1),"%"))
t1=Sys.time()
print(t1-t0)
