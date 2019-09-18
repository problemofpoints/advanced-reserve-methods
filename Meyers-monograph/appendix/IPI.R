#
# Script to run the IPI (CSR_CAY (CAY formerly known as CCL)) model
#       on data from the CAS Loss Reserve Database
# Uses Stan for the MCMC run
# by Glenn Meyers
# Run time on my computer < 2 minutes
#
####### Most of the script is in functions - They are called at the end.
#
rm(list = ls())  		# clear workspace
t0=Sys.time()
set.seed(12345)
#
# get packages
#
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(loo)
library(data.table)
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
# residual plot function for paid losses
#
Comb_CSR_resid_plot=function(co_model,grpcode,lineobus){
  lossData <- ins.line.data(CASdata,grpcode)
  train_data=subset(lossData,lossData$cal<11)
  w=train_data$w
  d=train_data$d
  Premium=train_data$net_premium[1:10]
  logloss=log(train_data$cum_pdloss)
  #
  std_resid=NULL
  ay=NULL
  dy=NULL
  samp=sample(1:length(co_model$logelr),100)
  for (s in samp){
    speedup=rep(1,10)
    for (i in 2:10){
      speedup[i]=speedup[i-1]*(1-co_model$gamma[s])
    }
    for (i in 1:length(logloss)){
      mu=log(Premium[w[i]])+co_model$logelr[s]+co_model$alpha[s,w[i]]+
        co_model$beta_p[s,d[i]]*speedup[w[i]]
      std_resid=c(std_resid,(logloss[i]-mu)/co_model$sig_p[s,d[i]])
      ay=c(ay,w[i])
      dy=c(dy,d[i])
    }
  }
  par(mfrow=c(1,2))
  AY=list(AY1=std_resid[ay==1],
          AY2=std_resid[ay==2],
          AY3=std_resid[ay==3],
          AY4=std_resid[ay==4],
          AY5=std_resid[ay==5],
          AY6=std_resid[ay==6],
          AY7=std_resid[ay==7],
          AY8=std_resid[ay==8],
          AY9=std_resid[ay==9],
          AY10=std_resid[ay==10])
  boxplot(AY,notch=T,col="gray",names=1:10,
          main="",xlab="Accident Year")
  abline(h=0,lwd=5)
  abline(h=qnorm(.25))
  abline(h=qnorm(.75))
  #
  DY=list(DY1=std_resid[dy==1],
          DY2=std_resid[dy==2],
          DY3=std_resid[dy==3],
          DY4=std_resid[dy==4],
          DY5=std_resid[dy==5],
          DY6=std_resid[dy==6],
          DY7=std_resid[dy==7],
          DY8=std_resid[dy==8],
          DY9=std_resid[dy==9],
          DY10=std_resid[dy==10])
  boxplot(DY,notch=T,col="gray",names=1:10,
          main="",xlab="Development Year")
  abline(h=0,lwd=5)
  abline(h=qnorm(.25))
  abline(h=qnorm(.75))
 mtext((paste("Combined CSR Model Standardized Residual BoxPlots",
              "\n",lineobus," Group",grpcode)),
              side = 3, line = -2.75, outer = TRUE)
}
#
# residual plot function for incurred losses
#

Comb_CAY_resid_plot=function(co_model,grpcode,lineobus){
  lossData <- ins.line.data(CASdata,grpcode)
  train_data=subset(lossData,lossData$cal<11)
  wne1=ifelse(train_data$w>1,1,0)
  w=train_data$w
  d=train_data$d
  Premium=train_data$net_premium[1:10]
  logloss=log(train_data$cum_incloss)
  #
  std_resid=NULL
  ay=NULL
  dy=NULL
  samp=sample(1:length(co_model$logelr),200)
  for(s in samp){
    mu=rep(0,length(logloss))
    mu[1]=log(Premium[1])+co_model$logelr[s]+co_model$alpha[s,w[1]]+
      co_model$beta_i[s,d[1]]
    resid1=(logloss[1]-mu[1])/co_model$sig_i[s,d[1]]
    ay=c(ay,w[1])
    dy=c(dy,d[1])
    std_resid=c(std_resid,resid1)
    for(i in 2:length(logloss)){
      mu[i]=log(Premium[w[i]])+co_model$logelr[s]+co_model$alpha[s,w[i]]+
        co_model$beta_i[s,d[i]]+
        co_model$rho[s]*(logloss[i-1]-mu[i-1])*wne1[i]
      resid1=(logloss[i]-mu[i])/co_model$sig_i[s,d[i]]
      ay=c(ay,w[i])
      dy=c(dy,d[i])
      std_resid=c(std_resid,resid1)
    }
  }

  par(mfrow=c(1,2))
  AY=list(AY1=std_resid[ay==1],
          AY2=std_resid[ay==2],
          AY3=std_resid[ay==3],
          AY4=std_resid[ay==4],
          AY5=std_resid[ay==5],
          AY6=std_resid[ay==6],
          AY7=std_resid[ay==7],
          AY8=std_resid[ay==8],
          AY9=std_resid[ay==9],
          AY10=std_resid[ay==10])
  boxplot(AY,notch=T,col="gray",names=1:10,
          main="",xlab="Accident Year")
  abline(h=0,lwd=5)
  abline(h=qnorm(.25))
  abline(h=qnorm(.75))
  #
  DY=list(DY1=std_resid[dy==1],
          DY2=std_resid[dy==2],
          DY3=std_resid[dy==3],
          DY4=std_resid[dy==4],
          DY5=std_resid[dy==5],
          DY6=std_resid[dy==6],
          DY7=std_resid[dy==7],
          DY8=std_resid[dy==8],
          DY9=std_resid[dy==9],
          DY10=std_resid[dy==10])
  boxplot(DY,notch=T,col="gray",names=1:10,
          main="",xlab="Development Year")
  abline(h=0,lwd=5)
  abline(h=qnorm(.25))
  abline(h=qnorm(.75))
  mtext((paste("Combined CAY Model Standardized Residual BoxPlots",
               "\n",lineobus," Group",grpcode)),
               side = 3, line = -2.75, outer = TRUE)
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
# stan script
#
CSR_CAYmodel_stan = "
data{
int<lower=1> len_data;
int<lower=0,upper=1> wne1[len_data];
real logprem[len_data];
real logcpd[len_data];
real loginc[len_data];
int<lower=1,upper=10> w[len_data];
int<lower=1,upper=10> d[len_data];
}
parameters{
real r_alpha[9];
real beta_p[10];
real r_beta_i[9];
real logelr;
real <lower=0,upper=100000> a_p[10];
real <lower=0,upper=100000> a_i[10];
real <lower=0,upper=1> r_rho;
real gamma;
}
transformed parameters{
real alpha[10];
real beta_i[10];
real speedup[10];
real rho;
real sig2_p[10];
real sig2_i[10];
real sig_p[10];
real sig_i[10];
real mu_p[len_data];
real mu_i[len_data];
alpha[1] = 0;
for (i in 2:10) alpha[i] = r_alpha[i-1];
for (i in 1:9) beta_i[i] = r_beta_i[i];
beta_i[10] = 0;
speedup[1] = 1;
for (i in 2:10) speedup[i] = speedup[i-1]*(1-gamma);
sig2_p[10] = gamma_cdf(1/a_p[10],1,1);
for (i in 1:9) sig2_p[10-i] = sig2_p[11-i]+gamma_cdf(1/a_p[i],1,1);
for (i in 1:10) sig_p[i] = sqrt(sig2_p[i]);
sig2_i[10] = gamma_cdf(1/a_i[10],1,1);
for (i in 1:9) sig2_i[10-i] = sig2_i[11-i]+gamma_cdf(1/a_i[i],1,1);
for (i in 1:10) sig_i[i] = sqrt(sig2_i[i]);
mu_p[1] = logprem[1]+logelr+beta_p[d[1]]*speedup[w[1]];
mu_i[1] = logprem[1]+logelr+beta_i[d[1]];
rho = -2*r_rho+1;
for (i in 2:len_data){
mu_p[i] = logprem[i]+logelr+alpha[w[i]]+beta_p[d[i]]*speedup[w[i]];
mu_i[i] = logprem[i]+logelr+alpha[w[i]]+beta_i[d[i]]+
rho*(loginc[i-1]-mu_i[i-1])*wne1[i];
}
}
model {
r_alpha ~ normal(0,3.162);
beta_p ~ normal(0,3.162);
r_beta_i ~ normal(0,3.162);
a_p ~ inv_gamma(1,1);
a_i ~ inv_gamma(1,1);
logelr ~ normal(-.4,3.162);
gamma ~ normal(0,0.05);
r_rho ~ beta(2,2);
for (i in 1:len_data) {
logcpd[i] ~ normal(mu_p[i],sig_p[d[i]]);
loginc[i] ~ normal(mu_i[i],sig_i[d[i]]);
}
}
generated quantities{
vector[len_data] log_lik_cpd;
vector[len_data] log_lik_inc;
for (i in 1:len_data){
log_lik_cpd[i] = normal_lpdf(logcpd[i]|mu_p[i],sig_p[d[i]]);
log_lik_inc[i] = normal_lpdf(loginc[i]|mu_i[i],sig_i[d[i]]);
}
}
"
#
# initialization function for CSR_CAYmodel_stan
#
init.CSR_CAY=function(chain_id){
  set.seed(12345)
  list(r_alpha=rnorm(9,0,0.2),r_beta_p=runif(9),
       r_beta_i=runif(9),a_i=runif(10),a_p=runif(10),
       logelr=runif(1,-0.75,-0.25),gamma=rnorm(1,0,0.1),r_rho=runif(1))
}
#
# dummy data for compiling
#
data.dummy=list(len_data = 55,
                logprem  = rep(8,55),
                logcpd   = rep(8,55),
                loginc   = rep(8,55),
                w        = c(1:10,1:9,1:8,1:7,1:6,1:5,1:4,1:3,1,2,1),
                d        = c(rep(1,10),rep(2,9),rep(3,8),rep(4,7),rep(5,6),
                             rep(6,5),rep(7,4),rep(8,3),rep(9,2),10),
                wne1     =rep(1,55))

#
# compile the univariate model
#
fitC_CSR_CAY = stan(model_code=CSR_CAYmodel_stan,data=data.dummy,
                    seed=12345,init=init.CSR_CAY,chains=0)
#
pars.list=c("logelr","alpha","beta_p","beta_i","gamma",
            "rho","sig_p","sig_i","log_lik_cpd","log_lik_inc")
#
# set up function to run stan model and create output
#
model_function=function(CASdata,grpcode){
  lossData <- ins.line.data(CASdata,grpcode)
  train_data=subset(lossData,lossData$cal<11)
  test_data=subset(lossData,lossData$cal>10)
  wne1=ifelse(train_data$w>1,1,0)
  Premium=train_data$net_premium[1:10]
  logprem=log(Premium)
  cpdloss=train_data$cum_pdloss
  incloss=train_data$cum_incloss
  cpdtest=test_data$cum_pdloss
  inctest=test_data$cum_incloss

  #
  #  data for the model
  #
  data.CSR_CAY=list(len_data = length(train_data$w),
                    logprem  = log(train_data$net_premium),
                    logcpd   = log(cpdloss),
                    loginc   = log(incloss),
                    w        = train_data$w,
                    d        = train_data$d,
                    wne1     = wne1)
  #
  # run the model
  #
  stan_thin=1
  stan_iter=5000
  Rhat_target=1.05
  max_Rhat=2
  while ((max_Rhat > Rhat_target)&(stan_thin<17)){
    fitCSR_CAY=stan(model_code=CSR_CAYmodel_stan,data=data.CSR_CAY,
                    init=init.CSR_CAY,
                    seed = 12345,iter=stan_iter,thin=stan_thin,
                    chains = 4,pars=pars.list,
                    control=list(adapt_delta=.9999,max_treedepth=50),
                    refresh=0)
    fitCSR_CAY_summary=as.matrix(summary(fitCSR_CAY)$summary)[,c(1,3,10)]
    mrh=subset(fitCSR_CAY_summary,is.na(fitCSR_CAY_summary[,3])==F)
    max_Rhat=round(max(mrh[,3]),4)
    print(paste("Maximum Rhat =",max_Rhat,"Thin =",stan_thin))
    stan_thin=2*stan_thin
    stan_iter=2*stan_iter
  }
  stan_thin=stan_thin/2
  #
  # extract information from stan output to process in R
  #
  b=extract_ordered(fitCSR_CAY)
  alpha=b$alpha
  beta_p=b$beta_p
  beta_i=b$beta_i
  gamma=b$gamma
  logelr=b$logelr
  rho=b$rho
  sig_p=b$sig_p
  sig_i=b$sig_i
  num.mcmc=length(logelr)
  #
  # goodness of fit statistics for comparing models
  #
  elpd_cpd=loo(b$log_lik_cpd)
  elpd_inc=loo(b$log_lik_inc)

  #
  # calculate test_elpd for test data
  #
  mu.test=matrix(0,num.mcmc,10)
  test_elpd_I=0
  trloss=as.triangle(lossData,origin="w",dev="d",value="cum_incloss")
  for (dev in 2:10){
    mu.test=matrix(0,num.mcmc,10)
    mu.test[,1]=log(Premium[1])+logelr+beta_i[,dev]
    for (origin in 2:10){
      mu.test[,origin]=log(Premium[origin])+logelr+alpha[,origin]+beta_i[,dev]+
        rho*(log(trloss[origin-1,dev])-mu.test[,origin-1])
      if (dev+origin>11){
        test_elpd_I=test_elpd_I+
          log(mean(dnorm(log(trloss[origin,dev]),mu.test[,origin],
                         sig_i[,dev],log=F)))
      }
    }
  }
  test_elpd_I=round(test_elpd_I,3)
  test_elpd_P=0
  for (i in 1:dim(test_data)[1]){
    mu.test=log(Premium[test_data$w[i]])+logelr+alpha[,test_data$w[i]]+
      beta_p[,test_data$d[i]]*(1-gamma)^(test_data$w[i]-1)
    test_elpd_P=test_elpd_P+
      log(mean(dnorm(log(cpdtest[i]),mu.test,sig_p[,test_data$d[i]],log=F)))
  }
  test_elpd_P=round(test_elpd_P,3)
  #
  # simulate outcomes for d=10
  #
  at_p.wd10=matrix(cpdloss[55],num.mcmc,10)
  for (w in 2:10){
    at_p.wd10[,w]=rlnorm(num.mcmc,log(Premium[w])+logelr+
                         alpha[,w]+beta_p[,10]*(1-gamma)^(w-1),
                         sig_p[,10])
  }
  at_i.wd10=matrix(incloss[55],num.mcmc,10)
  mu_i=matrix(0,num.mcmc,10)
  mu_i[,1]=log(Premium[1])+logelr+beta_i[,10]
  for (w in 2:10){
    mu_i[,w]=log(Premium[w])+logelr+alpha[,w]+
             rho*(log(at_i.wd10[,w-1])-mu_i[,w-1])
    at_i.wd10[,w]=rlnorm(num.mcmc,mu_i[,w],sig_i[,10])
  }
  ss_p.wd10=rep(0,10)
  ms_p.wd10=rep(0,10)
  ss_i.wd10=rep(0,10)
  ms_i.wd10=rep(0,10)
  #
  ms_p.wd10[1]=mean(at_p.wd10[,1])
  ms_i.wd10[1]=mean(at_p.wd10[,1])
  for (w in 2:10){
    ms_p.wd10[w]=mean(at_p.wd10[,w])
    ss_p.wd10[w]=sd(at_p.wd10[,w])
    ms_i.wd10[w]=mean(at_i.wd10[,w])
    ss_i.wd10[w]=sd(at_i.wd10[,w])
  }
  Pred.IP_CSR=rowSums(at_p.wd10)
  ms_p.td10=mean(Pred.IP_CSR)
  ss_p.td10=sd(Pred.IP_CSR)
  IP_CSR.Estimate=round(ms_p.wd10)
  IP_CSR.SE=round(ss_p.wd10)
  IP_CSR.CV=round(IP_CSR.SE/IP_CSR.Estimate,4)
  act_i=subset(lossData$cum_incloss,lossData$d==10)
  pct.IP_CSR=sum(Pred.IP_CSR<=act_i)/length(Pred.IP_CSR)*100
  #
  # put IP_CSR accident year statistics into a data frame
  #
  Pred.IP_CSR=rowSums(at_p.wd10[,1:10])
  ms_p.td10=mean(Pred.IP_CSR)
  ss_p.td10=sd(Pred.IP_CSR)
  IP_CSR.Estimate=round(ms_p.wd10)
  IP_CSR.SE=round(ss_p.wd10)
  IP_CSR.CV=round(IP_CSR.SE/IP_CSR.Estimate,4)
  act_p=sum(subset(lossData$cum_pdloss,lossData$d==10)[1:10])
  pct.IP_CSR=sum(Pred.IP_CSR<=act_p)/length(Pred.IP_CSR)*100
  IP_CSR.Estimate=c(IP_CSR.Estimate,round(ms_p.td10))
  IP_CSR.SE=c(IP_CSR.SE,round(ss_p.td10))
  IP_CSR.CV=c(IP_CSR.CV,round(ss_p.td10/ms_p.td10,4))
  Premium=c(Premium,sum(Premium))
  Outcome_P=subset(lossData$cum_pdloss,lossData$d==10)
  Outcome_P=c(Outcome_P,sum(Outcome_P))
  Group=rep(grpcode,11)
  IP_CSR.Pct=c(rep(NA,10),pct.IP_CSR)
  W=c(1:10,"Total")
  #
  # put IP_CAY accident year statistics into a data frame
  #
  Pred.IP_CAY=rowSums(at_i.wd10)
  ms_i.td10=mean(Pred.IP_CAY)
  ss_i.td10=sd(Pred.IP_CAY)
  IP_CAY.Estimate=round(ms_i.wd10)
  IP_CAY.SE=round(ss_i.wd10)
  IP_CAY.CV=round(IP_CAY.SE/IP_CAY.Estimate,4)
  act_i=sum(subset(lossData$cum_incloss,lossData$d==10)[1:10])
  pct.IP_CAY=sum(Pred.IP_CAY<=act_i)/length(Pred.IP_CAY)*100
  IP_CAY.Estimate=c(IP_CAY.Estimate,round(ms_i.td10))
  IP_CAY.SE=c(IP_CAY.SE,round(ss_i.td10))
  IP_CAY.CV=c(IP_CAY.CV,round(ss_i.td10/ms_i.td10,4))
  Outcome_I=subset(lossData$cum_incloss,lossData$d==10)
  Outcome_I=c(Outcome_I,sum(Outcome_I))
  Group=rep(grpcode,11)
  IP_CAY.Pct=c(rep(NA,10),pct.IP_CAY)
  W=c(1:10,"Total")
  risk_P=data.frame(Group,W,Premium,IP_CSR.Estimate,IP_CSR.SE,IP_CSR.CV,
                    Outcome_P,IP_CSR.Pct)
  risk_I=data.frame(Group,W,Premium,IP_CAY.Estimate,IP_CAY.SE,IP_CAY.CV,
                    Outcome_I,IP_CAY.Pct)
  risk_IP=data.frame(Group,W,Premium,
                     IP_CAY.Estimate,IP_CAY.SE,IP_CAY.CV,Outcome_I,
                     IP_CAY.Pct,
                     IP_CSR.Estimate,IP_CSR.SE,IP_CSR.CV,Outcome_P,
                     IP_CSR.Pct)
  #
  # set up output
  #
  Group=grpcode
  Premium=Premium[11]
  IP_CSR.Estimate=IP_CSR.Estimate[11]
  IP_CSR.SE=IP_CSR.SE[11]
  IP_CSR.CV=IP_CSR.CV[11]
  Outcome_I=Outcome_I[11]
  IP_CSR.Pct=IP_CSR.Pct[11]
  #
  IP_CAY.Estimate=IP_CAY.Estimate[11]
  IP_CAY.SE=IP_CAY.SE[11]
  IP_CAY.CV=IP_CAY.CV[11]
  Outcome_P=Outcome_P[11]
  IP_CAY.Pct=IP_CAY.Pct[11]
  mean_rho=round(mean(rho),4)
  sd_rho=round(sd(rho),4)
  mean_gamma=round(mean(gamma),4)
  sd_gamma=round(sd(gamma),4)
  mean_logelr=round(mean(logelr),4)
  elpd_P=round(elpd_cpd$estimates[1,1],3)
  elpd_I=round(elpd_inc$estimates[1,1],3)
  test_elpd_P=round(test_elpd_P,3)
  test_elpd_I=round(test_elpd_I,3)
  SumStats=data.frame(Group,Premium,
                      IP_CAY.Estimate,IP_CAY.SE,IP_CAY.CV,Outcome_I,
                      IP_CAY.Pct,mean_rho,sd_rho,
                      IP_CSR.Estimate,IP_CSR.SE,IP_CSR.CV,Outcome_P,
                      IP_CSR.Pct,mean_gamma,sd_gamma,mean_logelr,
                      elpd_I,elpd_P,test_elpd_I,test_elpd_P,
                      stan_thin)
  output=list(ParmSummary=fitCSR_CAY_summary[1:53,],
              alpha      =alpha,
              beta_i     =beta_i,
              beta_p     =beta_p,
              logelr     =logelr,
              rho        =rho,
              gamma      =gamma,
              sig_i      =sig_i,
              sig_p      =sig_p,
              risk_P     =risk_P,
              risk_I     =risk_I,
              risk_IP    =risk_IP,
              Pred_IP_CAY=Pred.IP_CAY,
              Pred_IP_CSR=Pred.IP_CSR,
              SumStats =SumStats,
              b        =b)
  return(output)
}
#
# Single triangle
#
CASdata = read.csv("~/Dropbox/CAS Loss Reserve Database/comauto_pos.csv")
# Location of files in the CAS Loss Reserve Database
#    http://www.casact.org/research/reserve_data/comauto_pos.csv
#    http://www.casact.org/research/reserve_data/ppauto_pos.csv
#    http://www.casact.org/research/reserve_data/wkcomp_pos.csv
#    http://www.casact.org/research/reserve_data/othliab_pos.csv
lineobus="CA - IP"
grpcode=353
co_model=model_function(CASdata,grpcode)
print(co_model$ParmSummary)
print(co_model$risk_P)
print(co_model$risk_I)
print(co_model$SumStats)
Comb_CSR_resid_plot(co_model,grpcode,lineobus)
Comb_CAY_resid_plot(co_model,grpcode,lineobus)

par(mfrow=c(2,1))
rng=range(co_model$Pred_IP_CAY,co_model$Pred_IP_CSR)
hist(co_model$Pred_IP_CAY,xlab="Simulated Outcomes",xlim=rng,
     main="Predictive Distribution of Outcomes for Incurred Losses",
     sub=paste("Mean =",co_model$SumStats$IP_CAY.Estimate,
               " SE =",co_model$SumStats$IP_CAY.SE))
hist(co_model$Pred_IP_CSR,xlab="Simulated Outcomes",xlim=rng,
     main="Predictive Distribution of Outcomes for Paid Losses",
     sub=paste("Mean =",co_model$SumStats$IP_CSR.Estimate,
               " SE =",co_model$SumStats$IP_CSR.SE))
t1=Sys.time()
print(t1-t0)