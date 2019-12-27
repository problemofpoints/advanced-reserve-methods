#
# R/Stan Script that fits the CSR-w model to data underling the
# Learning Lounge session given at the 2018 CAS Annual Meeting
# given by Bob Wolf and Mary Frances Miller
#
# Script written by Glenn Meyers
#
# Run time varies but should be less than 2 minutes with Thin=1
#
rm(list = ls())  		# clear workspace
t0=Sys.time()
#
# user input
#
run_id="P-7:16"
# First letter "P" for paid loss triangle
# First letter "I" for incurred loss triangle
keep_cal=7:16 # specifies calendar years for fitting
loss_type=substr(run_id,1,1)
#
# designate output files
#
setwd("~/Dropbox/LL CAS Annual Meeting/E-Forum Paper/Case_Study")
res.stats=NULL
#res.stats="LL_Reserve_Statistics.csv"
earlier_runs=NULL
#earlier_runs=read.csv(file=res.stats)
AY.Exhibit=NULL
#AY.Exhibit=paste(run_id,"AY_Exhibit_CSR-w.csv")
parm.summary=NULL
#parm.summary=paste(run_id,"Parm_Summary_CSR-w.csv")

#
# read in data
#
raw_paid=read.csv(file="LL_paid_triangle.csv")
raw_incurred=read.csv(file="LL_incurred_triangle.csv")
a_Premium=round(raw_paid[,2])
adata=NULL
rec=0
for (d in 1:dim(raw_paid)[1]){
  for (w in 1:(dim(raw_paid)[1]+1-d)){
    rec=rec+1
    ay=raw_paid[w,1]
    net_premium=round(raw_paid[w,2])
    ploss=round(raw_paid[w,d+2])
    iloss=round(raw_incurred[w,d+2])
    cal=w+d-1
    wne1=ifelse(w==1,1,0)
    adata=rbind(adata,c(rec,ay,net_premium,w,d,ploss,iloss,cal,wne1))
  }
}
colnames(adata)=c("rec","ay","net_premium","w","d","ploss","iloss","cal","wne1")
adata=as.data.frame(adata)
fit_rec=subset(adata$rec,adata$cal %in% keep_cal)
ip_loss=NULL
for (i in 1:dim(adata)[1]){
  ip_loss=c(ip_loss,ifelse(loss_type=="P",adata$ploss[i],adata$iloss[i]))
}
rdata=subset(adata,adata$cal %in% keep_cal)
W=max(rdata$w)
D=max(rdata$d)
C=max(rdata$cal)
Premium=rep(0,W)
for (w in adata$w){
  Premium[w]=min(subset(adata$net_premium,adata$w==w))
}
#
# get packages
#
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(loo)
library(data.table)
#
# Stan input list
#
data.CSR=list(len_data = length(ip_loss),
              len_fit  = length(fit_rec),
              fit_rec  = fit_rec,
              W        = W,
              D        = D,
              logprem  = log(adata$net_premium),
              logloss  = log(ip_loss),
              w        = adata$w,
              d        = adata$d)#
# Stan output parameters
#
pars.list=c("logelr","alpha","beta","sig","gamma","speedup","mu","log_lik")
#
# initialization function for CSRModel
#
init.CSR=function(chain_id){
  set.seed(12345)
  list(r_alpha=rnorm(W-1,0,0.2),r_beta=runif(D-1),a_ig=runif(D),
       logelr=runif(1,-0.75,-0.5),gamma=rnorm(1,0,0.1))
}
#
# Stan CSR script
#
CSRmodel_stan = "
data{
int<lower=1> len_data;
int<lower=1> len_fit;
int<lower=1> fit_rec[len_fit];
int<lower=1> W;
int<lower=1> D;
real logprem[len_data];
real logloss[len_data];
int<lower=1,upper=W> w[len_data];
int<lower=1,upper=D> d[len_data];
}
parameters{
real r_alpha[W-1];
real r_beta[D-1];
real logelr;
real <lower=0,upper=100000> a_ig[W];
real gamma;
}
transformed parameters{
real alpha[W];
real beta[D];
real speedup[W];
real sig2[D];
real sig[D];
real mu[len_data];
alpha[1] = 0;
for (i in 2:W) alpha[i] = r_alpha[i-1];
for (i in 1:(D-1)) beta[i] = r_beta[i];
beta[D] = 0;
speedup[1] = 1;
for (i in 2:W) speedup[i] = speedup[i-1]*(1-gamma);
sig2[D] = gamma_cdf(1/a_ig[D],1,1);
for (i in 1:(D-1)) sig2[D-i] = sig2[D+1-i]+gamma_cdf(1/a_ig[i],1,1);
for (i in 1:D) sig[i] = sqrt(sig2[i]);
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
for (i in 1:len_fit){
  logloss[fit_rec[i]] ~ normal(mu[fit_rec[i]],sig[d[fit_rec[i]]]);
}
}
generated quantities{
vector[len_fit] log_lik;
for (i in 1:len_fit) log_lik[i] =
    normal_lpdf(logloss[fit_rec[i]]|mu[fit_rec[i]],sig[d[fit_rec[i]]]);
}
"
#
# run the model
#
stan_thin=1
stan_iter=5000
Rhat_target=1.05
max_Rhat=2
while ((max_Rhat > Rhat_target)&(stan_thin<17)){
  fitCSR=stan(model_code=CSRmodel_stan, data = data.CSR,init=init.CSR,
              seed = 12345,iter=stan_iter,thin=stan_thin,
              chains = 4,pars=pars.list,
              control=list(adapt_delta=.9999,max_treedepth=50),
              refresh=0)
  fitCSR_summary=
    as.matrix(summary(fitCSR)$summary)[1:(2*W+2*D+2),c(1,3,10)]
  mrh=subset(fitCSR_summary,is.na(fitCSR_summary[,3])==F)
  max_Rhat=round(max(mrh[,3]),4)
  print(paste("Maximum Rhat =",max_Rhat,"Thin =",stan_thin))
  stan_thin=2*stan_thin
  stan_iter=2*stan_iter
}
stan_thin=stan_thin/2
#
# goodness of fit statistics for comparing models
# warning -- the loo package output depends upon the version of the loo package
#
loglik1=extract_log_lik(fitCSR)
CSR_loo=(loo(loglik1))
elpd_loo=round(CSR_loo$elpd_loo,3)
p_loo=round(CSR_loo$p_loo,3)
#
# extract information from stan output to process in R
#
b <- extract(fitCSR)
alpha=b$alpha
beta=b$beta
gamma=b$gamma
speedup=b$speedup
logelr=b$logelr
sig=b$sig
mu=b$mu
num.mcmc=length(logelr)
#
# calculate loss statistice by accident year for simulated outcomes
#
set.seed(12345)
at.wdD=matrix(0,num.mcmc,W)
mu.wdD=rep(0,num.mcmc)
for (w in 1:W){
  mu.wdD=log(Premium[w])+alpha[,w]+logelr
  at.wdD[,w]=exp(mu.wdD+sig[,D]^2/2) # this the estimate given mu and sigma
}
#
# calculate loss statistics and output to data frame
#
ss.wdD=rep(0,W)
ms.wdD=rep(0,W)
#
for (w in 1:W){
  ms.wdD[w]=mean(at.wdD[,w])
  ss.wdD[w]=sd(at.wdD[,w])
}
Pred.CSR=rowSums(at.wdD)
ms.tdD=mean(Pred.CSR)
ss.tdD=sd(Pred.CSR)
CSR.Estimate=round(ms.wdD)
CSR.SE=round(ss.wdD)
CSR.CV=round(CSR.SE/CSR.Estimate,4)
# append totals
CSR.Estimate=c(CSR.Estimate,round(ms.tdD))
CSR.SE=c(CSR.SE,round(ss.tdD))
CSR.CV=c(CSR.CV,round(ss.tdD/ms.tdD,4))
Premium=c(Premium,sum(Premium))
AY=c(sort(unique(adata$ay)),"Total")
AY_Exhibit=data.frame(AY,Premium,CSR.Estimate,CSR.SE,CSR.CV)
print(AY_Exhibit)
Premium=Premium[W+1]
CSR.Estimate=CSR.Estimate[W+1]
CSR.SE=CSR.SE[W+1]
mean_elr=round(mean(exp(logelr)),3)
SumStats=data.frame(Premium,CSR.Estimate,CSR.SE,
                    mean_elr,elpd_loo,p_loo,stan_thin)
print(fitCSR_summary)
print(SumStats)
#
# create histogram of reserve predictive distribution
#
par(mfrow=c(1,1))
reserve=sort(Pred.CSR-sum(rdata$ploss[rdata$cal==max(keep_cal)]))
CSR.ResMean=round(mean(reserve))
CSR.ResLow=round(reserve[250])
CSR.ResHigh=round(reserve[9750])
hist(reserve,
  main=paste("CSR-w",run_id, "Predictive Distribution of the Loss Reserve"),
  xlab=paste("Mean =",CSR.ResMean),
  sub=paste("Low End (2.5%) =",CSR.ResLow," High End (97.5%) =",CSR.ResHigh))
abline(v=CSR.ResLow,lwd=3,col="red")
abline(v=CSR.ResMean,lwd=3,col="blue")
abline(v=CSR.ResHigh,lwd=3,col="red")
#
# standardized residual plots
#
w=adata$w
d=adata$d
Premium=Premium[1:W]
logloss=log(ip_loss)
#
std_resid=NULL
ay=NULL
dy=NULL
cy=NULL
s_size=200
#
# function to equate sample sizes for all ay,dy and cy
#
eq_s=function(std_r,yr,size){
  eq_r=NULL
  eq_y=NULL
  for (y in unique(yr)){
    std_resid_y=subset(std_r,yr==y)
    yr_y=subset(yr,yr==y)
    ss=sample(1:length(yr_y),size)
    eq_r=c(eq_r,std_resid_y[ss])
    eq_y=c(eq_y,yr_y[ss])
  }
  return(data.frame(eq_r,eq_y))
}
samp=sample(1:length(logelr),s_size)
for (s in samp){
  for (i in 1:length(logloss)){
    if(i %in% fit_rec){
      std_resid=c(std_resid,(logloss[i]-mu[s,i])/sig[s,d[i]])
      ay=c(ay,w[i])
      dy=c(dy,d[i])
      cy=c(cy,w[i]+d[i]-1)
    }
  }
}
#rng=range(std_resid)
rng=c(-4.5,4.5)  #for the paper
par(mfrow=c(1,3))
temp=eq_s(std_resid,ay,s_size)
boxplot(split(temp$eq_r,temp$eq_y),notch=T,ylim=rng,
        col="gray",main="",xlab="Accident Year")
abline(h=0,lwd=5)
abline(h=qnorm(.25))
abline(h=qnorm(.75))
#
temp=eq_s(std_resid,dy,s_size)
boxplot(split(temp$eq_r,temp$eq_y),notch=T,ylim=rng,
        col="gray",main="",xlab="Development Year")
abline(h=0,lwd=5)
abline(h=qnorm(.25))
abline(h=qnorm(.75))
#
temp=eq_s(std_resid,cy,s_size)
boxplot(split(temp$eq_r,temp$eq_y),notch=T,ylim=rng,
        col="gray",main="",xlab="Calendar Year")
abline(h=0,lwd=5)
abline(h=qnorm(.25))
abline(h=qnorm(.75))
#
mtext(paste("CSR-w",run_id, "Standardized Residual Boxplots"),
      side = 3, line = -2.75, outer = TRUE)
#
# write output data files
#
Runid=paste("CSR-w",run_id,"")
Ultimate.Loss=CSR.Estimate
Ultimate.SE=CSR.SE
Reserve.Low=CSR.ResLow
Reserve.Mean=CSR.ResMean
Reserve.High=CSR.ResHigh
ELPD=elpd_loo
res_stats=rbind(earlier_runs[2:8],data.frame(Runid,Ultimate.Loss,Ultimate.SE,
                                             Reserve.Low,Reserve.Mean,Reserve.High,ELPD))
parm_summary=round(fitCSR_summary[,1:2],4)
if (is.null(res.stats)==F){
  res_stats=rbind(earlier_runs[2:8],data.frame(Runid,Ultimate.Loss,Ultimate.SE,
                                               Reserve.Low,Reserve.Mean,Reserve.High,ELPD))
  write.csv(res_stats,file=res.stats)
}
if (is.null(AY.Exhibit)==F) write.csv(AY_Exhibit,file=AY.Exhibit)
if (is.null(parm.summary)==F) write.csv(parm_summary,file=parm.summary)
t1=Sys.time()
print(t1-t0)
