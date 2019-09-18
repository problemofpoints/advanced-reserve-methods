#
# Script to run the CSR Dependency Model on data from the CAS Loss Reserve Database
# Uses Stan for the MCMC run
# by Glenn Meyers
# Run time on my computer - < 13 minutes
#
####### Most of the script is in functions - They are called at the end.
#
rm(list = ls())  		# clear workspace
t0=Sys.time()
#
# get packages
#
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(loo)
library(data.table)
library(mvtnorm)
library(parallel)
library(doParallel)
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
# residual plot function
#
bivariate_resid_plot=function(co_model,grpcode,linobus1,linobus2){
  lossData1 <- ins.line.data(CASdata1,grpcode)
  train_data1=subset(lossData1,lossData1$cal<11)
  lossData2 <- ins.line.data(CASdata2,grpcode)
  train_data2=subset(lossData2,lossData2$cal<11)
  logloss1=log(train_data1$cum_pdloss)
  logloss2=log(train_data2$cum_pdloss)
  w=train_data1$w
  d=train_data1$d
  #
  # create standardized residuals
  #
  std_resid1=NULL
  std_resid2=NULL
  ay=NULL
  dy=NULL
  samp=sample(1:length(co_model$rho),100)
  for (s in samp){
    for (i in 1:length(logloss1)){
      std_resid1=c(std_resid1,
                   (logloss1[i]-co_model$mu1[s,i])/co_model$sig1[s,d[i]])
      std_resid2=c(std_resid2,
                   (logloss2[i]-co_model$mu2[s,i])/co_model$sig2[s,d[i]])
      ay=c(ay,w[i])
      dy=c(dy,d[i])
    }
  }
  par(mfrow=c(2,1))
  #
  # scatter plot of standardized residuals
  #
  # plot(std_resid1,std_resid2,
  #      main=paste("Standardized Residuals - Group",grpcode),
  #      xlab=linobus1,ylab=linobus2)

  #
  # histogram of covariances by scenario
  #
  hist(co_model$rho,xlab=expression(rho),
       main=paste("CSR Model Correlations for",
                  "Group",grpcode,linobus1,"and",linobus2),
       sub=paste("Mean =",co_model$SumStats$mean_rho,
                 "Percentile of 0 =",co_model$SumStats$pct0))
  abline(v=0,lwd=3,col="blue")

  #
  # boxplot of standardized residuals by accident year
  #
  BoxPlotData=list(AY11=std_resid1[ay==1],
                   AY21=std_resid2[ay==1],
                   AY12=std_resid1[ay==2],
                   AY22=std_resid2[ay==2],
                   AY13=std_resid1[ay==3],
                   AY23=std_resid2[ay==3],
                   AY14=std_resid1[ay==4],
                   AY24=std_resid2[ay==4],
                   AY15=std_resid1[ay==5],
                   AY25=std_resid2[ay==5],
                   AY16=std_resid1[ay==6],
                   AY26=std_resid2[ay==6],
                   AY17=std_resid1[ay==7],
                   AY27=std_resid2[ay==7],
                   AY18=std_resid1[ay==8],
                   AY28=std_resid2[ay==8],
                   AY19=std_resid1[ay==9],
                   AY29=std_resid2[ay==9],
                   AY110=std_resid1[ay==10],
                   AY210=std_resid2[ay==10])
  #
  AY=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10)
  COL=rep(c("gray","black"),10)
  boxplot(BoxPlotData,notch=T,col=COL,names=AY,
          main="Standardized Residual Box Plots",
          sub=paste(linobus1,"in Gray ",linobus2,"in Black"),
          xlab="Accident Year")
  abline(h=0,lwd=5)
  abline(h=qnorm(.25))
  abline(h=qnorm(.75))
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
# initialization function for CSRModel
#
init.CSR=function(chain_id){
  set.seed(12345)
  list(r_alpha=rnorm(9,0,0.2),r_beta=runif(9),a_ig=runif(10),
       logelr=runif(1,-0.75,-0.5),gamma=rnorm(1,0,0.1))
}
#
# dummy data for compiling CSR model
#
CSR.data.dummy=list(len_data = 55,
                    logprem  = rep(8,55),
                    logloss  = rep(8,55),
                    w        = c(1:10,1:9,1:8,1:7,1:6,1:5,1:4,1:3,1,2,1),
                    d        = c(rep(1,10),rep(2,9),rep(3,8),rep(4,7),rep(5,6),
                                 rep(6,5),rep(7,4),rep(8,3),rep(9,2),10))
#
# compile the CSR model
#
fitC_CSR = stan(model_code=CSRmodel_stan,data=CSR.data.dummy,
                seed=12345,init=init.CSR,chains=0)
#
CSR_pars.list=c("logelr","alpha","beta","gamma","sig","mu","log_lik")
#
# stan rho script
#
rhomodel_stan="
data {
real<lower=0> sig_1[10];
real<lower=0> sig_2[10];
real mu_1[55];
real mu_2[55];
int<lower=1,upper=10> d[55];
vector[2] X[55];
}
transformed data{
vector[2] Mu[55];
for (i in 1:55){
  Mu[i][1]=mu_1[i];
  Mu[i][2]=mu_2[i];
}
}
parameters {
real<lower=0,upper=1> r_rho;
}
transformed parameters {
matrix[2,2] Sigma[55];
real rho;
rho=2*r_rho-1;
for (i in 1:55){
  Sigma[i][1,1]=sig_1[d[i]]^2;
  Sigma[i][1,2]=rho*sig_1[d[i]]*sig_2[d[i]];
  Sigma[i][2,1]=Sigma[i][1,2];
  Sigma[i][2,2]=sig_2[d[i]]^2;
}
}
model {
r_rho ~ beta(2,2);
for (i in 1:55){
  X[i] ~ multi_normal(Mu[i],Sigma[i]);
}
}
"
#
# set up function to run stan model and create output
#
model_function=function(CASdata1,CASdata2,grpcode){
  lossData1 <- ins.line.data(CASdata1,grpcode)
  train_data1=subset(lossData1,lossData1$cal<11)
  test_data1=subset(lossData1,lossData1$cal>10)
  loss1=train_data1$cum_pdloss
  logloss1=log(loss1)
  Premium1=train_data1$net_premium[1:10]
  w=train_data1$w
  d=train_data1$d
  #
  #  data for the model
  #
  data.CSR1=list(len_data = length(logloss1),
                 logprem  = log(train_data1$net_premium),
                 logloss  = logloss1,
                 w        = w,
                 d        = d)
  #
  # run the model
  #
  stan_thin=1
  stan_iter=5000
  Rhat_target=1.05
  max_Rhat=2
  while ((max_Rhat > Rhat_target)&(stan_thin<17)){
    fitCSR1=stan(fit = fitC_CSR, data = data.CSR1,init=init.CSR,
                 seed = 12345,iter=stan_iter,thin=stan_thin,
                 chains = 4,pars=CSR_pars.list,
                 control=list(adapt_delta=.9999,max_treedepth=50),
                 refresh=0)
    fitCSR1_summary=
      as.matrix(summary(fitCSR1)$summary)[1:32,c(1,3,10)]
    mrh=subset(fitCSR1_summary,is.na(fitCSR1_summary[,3])==F)
    max_Rhat=round(max(mrh[,3]),4)
    print(paste("Maximum Rhat =",max_Rhat,"Thin =",stan_thin))
    stan_thin=2*stan_thin
    stan_iter=2*stan_iter
  }
  stan_thin1=stan_thin/2
  #
  # extract information from stan output to process in R
  #
  b <- extract_ordered(fitCSR1)
  mu1=b$mu
  sig1=b$sig
  alpha1=b$alpha
  beta1=b$beta
  logelr1=b$logelr
  gamma1=b$gamma
  #
  # get mu and sig parameters for line 2
  #
  lossData2 <- ins.line.data(CASdata2,grpcode)
  train_data2=subset(lossData2,lossData2$cal<11)
  test_data2=subset(lossData2,lossData2$cal>10)
  Premium2=train_data2$net_premium[1:10]
  loss2=train_data2$cum_pdloss
  logloss2=log(loss2)
  #
  #  data for the model
  #
  data.CSR2=list(len_data = length(logloss2),
                 logprem  = log(train_data2$net_premium),
                 logloss  = logloss2,
                 w        = w,
                 d        = d)
  #
  # run the model
  #
  stan_thin=1
  stan_iter=5000
  Rhat_target=1.05
  max_Rhat=2
  while ((max_Rhat > Rhat_target)&(stan_thin<17)){
    fitCSR2=stan(fit = fitC_CSR, data = data.CSR2,init=init.CSR,
                 seed = 12345,iter=stan_iter,thin=stan_thin,
                 chains = 4,pars=CSR_pars.list,
                 control=list(adapt_delta=.9999,max_treedepth=50),
                 refresh=0)
    fitCSR2_summary=
      as.matrix(summary(fitCSR2)$summary)[1:32,c(1,3,10)]
    mrh=subset(fitCSR2_summary,is.na(fitCSR2_summary[,3])==F)
    max_Rhat=round(max(mrh[,3]),4)
    print(paste("Maximum Rhat =",max_Rhat,"Thin =",stan_thin))
    stan_thin=2*stan_thin
    stan_iter=2*stan_iter
  }
  stan_thin2=stan_thin/2
  #
  # extract information from stan output to process in R
  #
  b <- extract_ordered(fitCSR2)
  mu2=b$mu
  sig2=b$si
  alpha2=b$alpha
  beta2=b$beta
  logelr2=b$logelr
  gamma2=b$gamma
  num.mcmc=dim(mu2)[1]
  #
  # compile the rho model
  #
  init.rho=list(list(r_rho=0.5))
  rho.pars.list="rho"
  rho.dummy.data=list(sig_1    = rep(0.1,10),
                      sig_2    = rep(0.1,10),
                      mu_1     = rep(8,55),
                      mu_2     = rep(8,55),
                      d        = CSR.data.dummy$d,
                      X        = cbind(rep(8,55),rep(8,55)))
  fitC_rho = stan(model_code=rhomodel_stan,data=rho.dummy.data,
                  seed=12345,init=init.rho,chains=0)
  #
  # set up parallel processing to get the distribution of rho
  #
  cl <- makePSOCKcluster(4)
  registerDoParallel(cl)
  rho.mcmc=foreach (k=1:num.mcmc,.combine=rbind) %dopar%{
    library(rstan)
    rho.data=list(sig_1    = sig1[k,],
                  sig_2    = sig2[k,],
                  mu_1     = mu1[k,],
                  mu_2     = mu2[k,],
                  d        = d,
                  X        = cbind(logloss1,logloss2))
    #
    fit_rho = stan(fit=fitC_rho,data=rho.data,chains=1,iter=101,warmup=100,
                seed=k,init=init.rho,verbose=F,pars=rho.pars.list,
                control=list(adapt_delta=.99,max_treedepth=50),
                refresh=0)
    #
    # extract information from MCMC output to process in R
    #
    rho=extract(fit_rho)$rho[1]
    #
    # Rhat is commmented out after experimenting as it doubles the run time
    # Rhat=as.matrix(summary(fit1)$summary)[1,10]
    # rho.mcmc.out=c(rho,Rhat)
  }
  stopCluster(cl)
  #
  # calculate psis-loo model selection statistics using parallel processing
  #
  cl <- makePSOCKcluster(4)
  registerDoParallel(cl)
  loglik.biv.ind=foreach (k=1:num.mcmc,.combine=rbind) %dopar%{
    library(mvtnorm)
    lp=rep(0,55)
    lp0=rep(0,55)
    Mu=rep(0,2)
    Sigma=matrix(0,2,2)
    Sigma0=matrix(0,2,2)
    for (i in 1:55){
      Mu[1]=mu1[k,i]
      Mu[2]=mu2[k,i]
      Sigma[1,1]=sig1[k,d[i]]^2
      Sigma[2,1]=rho.mcmc[k]*sig1[k,d[i]]*sig2[k,d[i]]
      Sigma[1,2]=Sigma[2,1]
      Sigma[2,2]=sig2[k,d[i]]^2
      lp[i]=dmvnorm(c(logloss1[i],logloss2[i]),Mu,Sigma,log=T)
      Sigma0=Sigma
      Sigma0[1,2]=0
      Sigma0[2,1]=0
      lp0[i]=dmvnorm(c(logloss1[i],logloss2[i]),Mu,Sigma0,log=T)
    }
    lp1=c(lp,lp0)
    return(lp1)
  }
  stopCluster(cl)
  #
  loglik.biv=loglik.biv.ind[,1:55]
  loglik.ind=loglik.biv.ind[,56:110]
  loo.biv<-loo(loglik.biv)
  loo.ind <- loo(loglik.ind)
  #
  # calculate holdout/lower triangle model selection statistics
  #
  hw=test_data1$w
  hd=test_data1$d
  hlogloss1=log(test_data1$cum_pdloss)
  hlogloss2=log(test_data2$cum_pdloss)
  #
  cl <- makePSOCKcluster(4)
  registerDoParallel(cl)
  loglik.h.biv.ind=foreach (k=1:num.mcmc,.combine=rbind) %dopar%{
    library(mvtnorm)
    lp=rep(0,45)
    lp0=rep(0,45)
    Mu=rep(0,2)
    Sigma=matrix(0,2,2)
    Sigma0=matrix(0,2,2)
    for (i in 1:45){
      Mu[1]=log(Premium1[hw[i]])+logelr1[k]+alpha1[k,hw[i]]+
            beta1[k,hd[i]]*(1-gamma1[k])^(hw[i]-1)
      Mu[2]=log(Premium2[hw[i]])+logelr2[k]+alpha2[k,hw[i]]+
        beta2[k,hd[i]]*(1-gamma2[k])^(hw[i]-1)
      Sigma[1,1]=sig1[k,hd[i]]^2
      Sigma[2,1]=rho.mcmc[k]*sig1[k,hd[i]]*sig2[k,hd[i]]
      Sigma[1,2]=Sigma[2,1]
      Sigma[2,2]=sig2[k,hd[i]]^2
      lp[i]=dmvnorm(c(hlogloss1[i],hlogloss2[i]),Mu,Sigma,log=T)
      Sigma0=Sigma
      Sigma0[1,2]=0
      Sigma0[2,1]=0
      lp0[i]=dmvnorm(c(hlogloss1[i],hlogloss2[i]),Mu,Sigma0,log=T)
    }
    lp1=c(lp,lp0)
    return(lp1)
  }
  stopCluster(cl)
  biv_elpd_test=sum(log(colMeans(exp(loglik.h.biv.ind[,1:45]))))
  ind_elpd_test=sum(log(colMeans(exp(loglik.h.biv.ind[,46:90]))))
  #
  # simulate loss statistics for the independent case
  #
  # line 1
  #
  set.seed(12345)
  mu1.wd10=rep(0,num.mcmc)
  at1.wd10.Ind=matrix(loss1[55],num.mcmc,10)
  for (w in 2:10){
    mu1.wd10=log(Premium1[w])+alpha1[,w]+logelr1+beta1[,10]*(1-gamma1)^(w-1)
    at1.wd10.Ind[,w]=ceiling(rlnorm(num.mcmc,mu1.wd10,sig1[,10]))
  }
  #
  # line 2
  #
  set.seed(12345)
  mu2.wd10=rep(0,num.mcmc)
  at2.wd10.Ind=matrix(loss2[55],num.mcmc,10)
  for (w in 2:10){
    mu2.wd10=log(Premium2[w])+alpha2[,w]+logelr2+beta2[,10]*(1-gamma2)^(w-1)
    at2.wd10.Ind[,w]=ceiling(rlnorm(num.mcmc,mu2.wd10,sig2[,10]))
  }
  #
  ss.wd10.Ind=rep(0,10)
  ms.wd10.Ind=rep(0,10)
  #
  ms.wd10.Ind[1]=mean(at1.wd10.Ind[,1]+at2.wd10.Ind[,1])
  for (w in 2:10){
    ms.wd10.Ind[w]=mean(at1.wd10.Ind[,w]+at2.wd10.Ind[,w])
    ss.wd10.Ind[w]=sd(at1.wd10.Ind[,w]+at2.wd10.Ind[,w])
  }
  Pred.CSR.Ind=rowSums(at1.wd10.Ind+at2.wd10.Ind)
  ms.td10.Ind=mean(Pred.CSR.Ind)
  ss.td10.Ind=sd(Pred.CSR.Ind)
  CSR.Estimate.Ind=round(ms.wd10.Ind)
  CSR.SE.Ind=round(ss.wd10.Ind)
  CSR.CV.Ind=round(CSR.SE.Ind/CSR.Estimate.Ind,4)
  act=subset(lossData1$cum_pdloss,lossData1$d==10)
  act=act+subset(lossData2$cum_pdloss,lossData2$d==10)
  pct.CSR.Ind=sum(Pred.CSR.Ind<=sum(act))/length(Pred.CSR.Ind)*100
  # put CSR accident year statistics into a data frame
  #
  CSR.Estimate.Ind=c(CSR.Estimate.Ind,round(ms.td10.Ind))
  CSR.SE.Ind=c(CSR.SE.Ind,round(ss.td10.Ind))
  CSR.CV.Ind=c(CSR.CV.Ind,round(ss.td10.Ind/ms.td10.Ind,4))
  Premium=c(Premium1+Premium2,sum(Premium1+Premium2))
  Outcome=act
  Outcome=c(Outcome,sum(Outcome))
  Group=rep(grpcode,11)
  CSR.Pct.Ind=c(rep(NA,10),pct.CSR.Ind)
  W=c(1:10,"Total")
  risk.Ind=data.frame(Group,W,Premium,CSR.Estimate.Ind,CSR.SE.Ind,CSR.CV.Ind,
                      Outcome,CSR.Pct.Ind)
  #
  # simulate loss statistics for the dependent case
  #
  cl <- makePSOCKcluster(4)
  registerDoParallel(cl)
  at.wd10.Dep=foreach (k=1:num.mcmc,.combine=rbind) %dopar%{
    set.seed(12345+k)
    library(mvtnorm)
    atb1=rep(0,10)
    atb2=rep(0,10)
    Mu=rep(0,2)
    Sigma=matrix(0,2,2)
    Sigma[1,1]=sig1[k,10]^2
    Sigma[2,1]=rho.mcmc[k]*sig1[k,10]*sig2[k,10]
    Sigma[1,2]=Sigma[2,1]
    Sigma[2,2]=sig2[k,10]^2
    atb1[1]=loss1[55]
    atb2[1]=loss2[55]
    for (w in 2:10){
      Mu=c(log(Premium1[w])+logelr1[k]+alpha1[k,w],
           log(Premium2[w])+logelr2[k]+alpha2[k,w])
      lloss=rmvnorm(1,Mu,Sigma)
      atb1[w]=exp(lloss[1])
      atb2[w]=exp(lloss[2])
    }
    at=c(atb1,atb2)
  }
  stopCluster(cl)
  #
  #
  ss.wd10.Dep=rep(0,10)
  ms.wd10.Dep=rep(0,10)
  #
  ms.wd10.Dep[1]=mean(at.wd10.Dep[,1]+at.wd10.Dep[,11])
  for (w in 2:10){
    ms.wd10.Dep[w]=mean(at.wd10.Dep[,w]+at.wd10.Dep[,10+w])
    ss.wd10.Dep[w]=sd(at.wd10.Dep[,w]+at.wd10.Dep[,10+w])
  }
  Pred.CSR.Dep=rowSums(at.wd10.Dep)
  ms.td10.Dep=mean(Pred.CSR.Dep)
  ss.td10.Dep=sd(Pred.CSR.Dep)
  CSR.Estimate.Dep=round(ms.wd10.Dep)
  CSR.SE.Dep=round(ss.wd10.Dep)
  CSR.CV.Dep=round(CSR.SE.Dep/CSR.Estimate.Dep,4)
  act=subset(lossData1$cum_pdloss,lossData1$d==10)
  act=act+subset(lossData2$cum_pdloss,lossData2$d==10)
  pct.CSR.Dep=sum(Pred.CSR.Dep<=sum(act))/length(Pred.CSR.Dep)*100 #
  # put CSR accident year statistics into a data frame
  #
  CSR.Estimate.Dep=c(CSR.Estimate.Dep,round(ms.td10.Dep))
  CSR.SE.Dep=c(CSR.SE.Dep,round(ss.td10.Dep))
  CSR.CV.Dep=c(CSR.CV.Dep,round(ss.td10.Dep/ms.td10.Dep,4))
  Premium=c(Premium1+Premium2,sum(Premium1+Premium2))
  Outcome=act
  Outcome=c(Outcome,sum(Outcome))
  Group=rep(grpcode,11)
  CSR.Pct.Dep=c(rep(NA,10),pct.CSR.Dep)
  W=c(1:10,"Total")
  risk.Dep=data.frame(Group,W,Premium,CSR.Estimate.Dep,CSR.SE.Dep,CSR.CV.Dep,
                      Outcome,CSR.Pct.Dep)
  #
  mean_rho=round(mean(rho.mcmc),3)
  pct0=sum(rho.mcmc<=0)/num.mcmc*100
  Group=grpcode
  Premium=Premium[11]
  CSR.Estimate.Ind=CSR.Estimate.Ind[11]
  CSR.SE.Ind=CSR.SE.Ind[11]
  Outcome=Outcome[11]
  CSR.Pct.Ind=CSR.Pct.Ind[11]
  CSR.Estimate.Dep=CSR.Estimate.Dep[11]
  CSR.SE.Dep=CSR.SE.Dep[11]
  CSR.Pct.Dep=CSR.Pct.Dep[11]
  #
  biv_elpd_test=round(biv_elpd_test,3)
  ind_elpd_test=round(ind_elpd_test,3)
  #
  elpd_loo.Dep=round(loo.biv$estimates[1,1],3)
  p_loo.Dep=round(loo.biv$estimates[2,1],3)
  #
  elpd_loo.Ind=round(loo.ind$estimates[1,1],3)
  p_loo.Ind=round(loo.ind$estimates[2,1],3)
  #

  SumStats=data.frame(Group,Premium,
                      CSR.Estimate.Ind,CSR.SE.Ind,
                      CSR.Estimate.Dep,CSR.SE.Dep,
                      Outcome,CSR.Pct.Ind,CSR.Pct.Dep,
                      mean_rho,pct0,
                      biv_elpd_test,ind_elpd_test,
                      elpd_loo.Dep,p_loo.Dep,
                      elpd_loo.Ind,p_loo.Ind,
                      stan_thin1,stan_thin2)
  output=list(risk.Ind=risk.Ind,
              risk.Dep=risk.Dep,
              mu1=mu1,
              mu2=mu2,
              sig1=sig1,
              sig2=sig2,
              rho=rho.mcmc,
              SumStats=SumStats)
  return(output)
}
#
# Single insurer
#
CASdata1 = read.csv("~/Dropbox/CAS Loss Reserve Database/comauto_pos.csv")
CASdata2 = read.csv("~/Dropbox/CAS Loss Reserve Database/ppauto_pos.csv")
# Location of files in the CAS Loss Reserve Database
#    http://www.casact.org/research/reserve_data/comauto_pos.csv
#    http://www.casact.org/research/reserve_data/ppauto_pos.csv
#    http://www.casact.org/research/reserve_data/wkcomp_pos.csv
#    http://www.casact.org/research/reserve_data/othliab_pos.csv
linobus1="CA"
linobus2="PA"
grpcode=353
co_model=model_function(CASdata1,CASdata2,grpcode)
print(co_model$risk.Ind)
print(co_model$risk.Dep)
print(co_model$SumStats)
bivariate_resid_plot(co_model,grpcode,linobus1,linobus2)
t1=Sys.time()
print(t1-t0)
