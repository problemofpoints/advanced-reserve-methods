#
# Script to run the CSR Model on data from the CAS Loss Reserve Database
# Uses Stan for the MCMC run
# by Glenn Meyers
# Run time on my computer - < 1 minute
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
library(actuar)
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
# residual plot function
#
CSR_resid_plot=function(co_model,grpcode,lineobus){
  lossData <- ins.line.data(CASdata,grpcode)
  train_data=subset(lossData,lossData$cal<11)
  wne1=ifelse(train_data$w>1,1,0)
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
        co_model$beta[s,d[i]]*speedup[w[i]]
      std_resid=c(std_resid,(logloss[i]-mu)/co_model$sig[s,d[i]])
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
 mtext((paste("CSR Model Standardized Residual Box Plots",
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
pars.list=c("logelr","alpha","beta","gamma","sig","log_lik")
#
# set up function to run stan model and create output
#
model_function=function(CASdata,grpcode){
  lossData <- ins.line.data(CASdata,grpcode)
  train_data=subset(lossData,lossData$cal<11)
  test_data=subset(lossData,lossData$cal>10)
  Premium=train_data$net_premium[1:10]
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
  # goodness of fit statistics for comparing models
  #
  loglik1=extract_log_lik(fitCSR)
  elpd_stats=loo(loglik1)
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
  # calculate test_elpd for test data
  #
  test_elpd=0
  for (i in 1:dim(test_data)[1]){
    mu.test=log(Premium[test_data$w[i]])+logelr+alpha[,test_data$w[i]]+
      beta[,test_data$d[i]]*(1-gamma)^(test_data$w[i]-1)
    test_elpd=test_elpd+
      log(mean(dnorm(log(test[i]),mu.test,sig[,test_data$d[i]],log=F)))
  }
  #
  # simulate loss statistics by accident year
  #
  set.seed(12345)
  at.wd10=matrix(loss[55],num.mcmc,10)
  mu.wd10=rep(0,num.mcmc)
  for (w in 2:10){
    mu.wd10=log(Premium[w])+alpha[,w]+logelr
    at.wd10[,w]=ceiling(rlnorm(num.mcmc,mu.wd10,sig[,10]))
  }
  #
  # calculate loss statistics and output to data frame
  #
  ss.wd10=rep(0,10)
  ms.wd10=rep(0,10)
  #
  ms.wd10[1]=mean(at.wd10[,1])
  for (w in 2:10){
    ms.wd10[w]=mean(at.wd10[,w])
    ss.wd10[w]=sd(at.wd10[,w])
  }
  Pred.CSR=rowSums(at.wd10)
  ms.td10=mean(Pred.CSR)
  ss.td10=sd(Pred.CSR)
  CSR.Estimate=round(ms.wd10)
  CSR.SE=round(ss.wd10)
  CSR.CV=round(CSR.SE/CSR.Estimate,4)
  act=subset(lossData$cum_pdloss,lossData$d==10)
  sumact=sum(act)
  pct.CSR=sum(Pred.CSR<=sumact)/length(Pred.CSR)*100 #
  # put CSR accident year statistics into a data frame
  #
  CSR.Estimate=c(CSR.Estimate,round(ms.td10))
  CSR.SE=c(CSR.SE,round(ss.td10))
  CSR.CV=c(CSR.CV,round(ss.td10/ms.td10,4))
  Premium=c(Premium,sum(Premium))
  Outcome=act
  Outcome=c(Outcome,sum(Outcome))
  Group=rep(grpcode,11)
  CSR.Pct=c(rep(NA,10),pct.CSR)
  W=c(1:10,"Total")
  risk=data.frame(Group,W,Premium,CSR.Estimate,CSR.SE,CSR.CV,Outcome,CSR.Pct)
  Group=grpcode
  Premium=Premium[11]
  CSR.Estimate=CSR.Estimate[11]
  CSR.SE=CSR.SE[11]
  Outcome=Outcome[11]
  CSR.Pct=CSR.Pct[11]
  mean_gamma=round(mean(gamma),4)
  mean_logelr=round(mean(logelr),4)
  elpd_loo=round(elpd_stats$estimates[1,1],3)
  p_loo=round(elpd_stats$estimates[2,1],3)
  test_elpd=round(test_elpd,3)
  SumStats=data.frame(Group,Premium,CSR.Estimate,CSR.SE,Outcome,CSR.Pct,
                      mean_gamma,mean_logelr,stan_thin,
                      elpd_loo,p_loo,test_elpd)

  output=list(risk=risk,
              Predictive_Outcome=Pred.CSR,
              ParmSummary=fitCSR_summary,
              logelr=logelr,
              alpha=alpha,
              beta=beta,
              gamma=gamma,
              sig=sig,
              SumStats=SumStats)
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
lineobus="CA"
grpcode=353
co_model=model_function(CASdata,grpcode)
print(co_model$ParmSummary)
print(" ")
print(co_model$risk)
print("")
print(co_model$SumStats)
CSR_resid_plot(co_model,grpcode,lineobus)
par(mfrow=c(1,1))
hist(co_model$gamma,main="",#paste("Insurer",grpcode,lineobus),
     xlab=expression(gamma))
abline(v=0,col="blue",lwd=3)
hist(co_model$Predictive_Outcome,xlab="Simulated Outcomes",
     main="Predictive Distribution of Outcomes",
     sub=paste("Mean =",co_model$SumStats$CSR.Estimate,
               " SE =",co_model$SumStats$CSR.SE))
abline(v=co_model$SumStats$Outcome,lwd=3,col="red")
t1=Sys.time()
print(t1-t0)
