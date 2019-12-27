#
# R/Stan Script that fits the POS model to data underling the
# Learning Lounge session given at the 2018 CAS Annual Meeting
# given by Bob Wolf and Mary Frances Miller
#
# Script written by Glenn Meyers
#
# Run time varies but should be less than 3 minutes with Thin=1
#
rm(list = ls())  		# clear workspace
t0=Sys.time()
#
# user input
#
run_id="PI-7:16"
# First letter "P" for paid loss triangle
# First letter "I" for incurred loss triangle
keep_cal=7:16 # specifies calendar years for fitting
#
# designate output files
#
setwd("~/Dropbox/LL CAS Annual Meeting/E-Forum Paper/Case_Study")
res.stats=NULL
#res.stats="LL_Reserve_Statistics_POS.csv"
earlier_runs=NULL
#earlier_runs=read.csv(file=res.stats)
AY.Exhibit=NULL
#AY.Exhibit=paste(run_id,"AY_Exhibit_POS_vc.csv")
parm.summary=NULL
parm.summary=paste(run_id,"Parm_Summary_POS_vc.csv")
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
    ay=raw_paid[w,1]
    net_premium=round(raw_paid[w,2])
    ploss=round(raw_paid[w,d+2])
    iloss=round(raw_incurred[w,d+2])
    cal=w+d-1
    wne1=ifelse(w==1,1,0)
    adata=rbind(adata,c(ay,net_premium,w,d,ploss,iloss,cal,wne1))
  }
}
colnames(adata)=c("ay","net_premium","w","d","ploss","iloss","cal","wne1")
adata=as.data.frame(adata)
adata=subset(adata,adata$cal<=max(keep_cal))
rec=1:dim(adata)[1]
adata=cbind(rec,adata)
fit_rec=subset(adata$rec,adata$cal %in% keep_cal)
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
data.POS=list(len_data = length(adata$w),
              len_fit  = length(fit_rec),
              fit_rec  = fit_rec,
              W        = W,
              D        = D,
              C        = C,
              c        = min(keep_cal),
              logprem  = log(adata$net_premium),
              logiloss = log(adata$iloss),
              logploss = log(adata$ploss),
              w        = adata$w,
              d        = adata$d,
              cal      = adata$cal)
#
# Stan output parameters
#
pars.list=c("logelr","alpha","beta_p","beta_os","sig_p","sig_i",
            "gamma_pd","gamma_os",
            "speedup_p","speedup_os","mu_i","mu_p",
            "log_lik_cpd","log_lik_inc")
#
# initialization function for POSModel
#
init.POS=function(chain_id){
  set.seed(12345)
  list(r_alpha=rnorm(W-1,0,0.2),r_beta_p=runif(D-1),r_beta_os=runif(D),
       a_ig_p=runif(D),a_ig_i=runif(D),
       logelr=runif(1,-0.75,-0.5),
       r_gamma_pd=rnorm(C-1,0,0.025),r_gamma_os=rnorm(C-1,0,0.025))
}
#
# Stan POS script
#
POSmodel_stan = "
data{
int<lower=1> len_data;
int<lower=1> len_fit;
int<lower=1> fit_rec[len_fit];
int<lower=1> W;
int<lower=1> D;
int<lower=1> C;
int<lower=1> c;
real logprem[len_data];
real logiloss[len_data];
real logploss[len_data];
int<lower=1,upper=W> w[len_data];
int<lower=1,upper=D> d[len_data];
int<lower=1,upper=C> cal[len_data];
}
parameters{
real r_alpha[W-1];
real r_beta_p[D-1];
real beta_os[D];
real r_gamma_pd[C-1];
real r_gamma_os[C-1];
real logelr;
real <lower=0,upper=100000> a_ig_p[D];
real <lower=0,upper=100000> a_ig_i[D];
}
transformed parameters{
real alpha[W];
real beta_p[D];
real gamma_pd[C-1];  //changed naming conventions as Stan had a prior nane
real gamma_os[C-1];
real speedup_p[C];
real speedup_os[C];
real sig2_p[D];
real sig2_i[D];
real sig_p[D];
real sig_i[D];
real mu_p[len_data];
real mu_i[len_data];
alpha[1] = 0;
for (i in 2:W) alpha[i] = r_alpha[i-1];
//
//paid loss parameters
//
for (i in 1:(D-1)) beta_p[i] = r_beta_p[i];
beta_p[D] = 0;
for (i in 1:(C-1))
 if (i<c) gamma_pd[i]=0;
 else gamma_pd[i]=r_gamma_pd[i-c+1];
speedup_p[C]=1;
for (i in 1:(C-1)) speedup_p[C-i]=speedup_p[C-i+1]*(1+gamma_pd[C-i]);
sig2_p[D] = gamma_cdf(1/a_ig_p[D],1,1);
for (i in 1:(D-1)) sig2_p[D-i] = sig2_p[D+1-i]+gamma_cdf(1/a_ig_p[i],1,1);
for (i in 1:D) sig_p[i] = sqrt(sig2_p[i]);
for (i in 1:len_data){
mu_p[i] = logprem[i]+logelr+alpha[w[i]]+beta_p[d[i]]*speedup_p[cal[i]];
}
//
// outstanding loss parameters
//
for (i in 1:(C-1))
 if (i<c) gamma_os[i]=0;
 else gamma_os[i]=r_gamma_os[i-c+1];
speedup_os[C]=1;
for (i in 1:(C-1)) speedup_os[C-i]=speedup_os[C-i+1]*(1+gamma_os[C-i]);
sig2_i[D] = gamma_cdf(1/a_ig_i[D],1,1);
for (i in 1:(D-1)) sig2_i[D-i] = sig2_i[D+1-i]+gamma_cdf(1/a_ig_i[i],1,1);
for (i in 1:D) sig_i[i] = sqrt(sig2_i[i]);
for (i in 1:len_data){
mu_i[i] = mu_p[i]+beta_os[d[i]]*speedup_os[cal[i]];
}
}
model {
r_alpha ~ normal(0,3.162);
r_beta_p ~ normal(0,3.162);
beta_os ~ normal(0,3.162);
a_ig_p ~ inv_gamma(1,1);
a_ig_i ~ inv_gamma(1,1);
logelr ~ normal(-.4,3.162);
r_gamma_pd ~ normal(0,0.05);
r_gamma_os ~ normal(0,0.05);
for (i in 1:len_fit){
logploss[fit_rec[i]] ~ normal(mu_p[fit_rec[i]],sig_p[d[fit_rec[i]]]);
logiloss[fit_rec[i]] ~ normal(mu_i[fit_rec[i]],sig_i[d[fit_rec[i]]]);
}
}
generated quantities{
vector[len_fit] log_lik_cpd;
vector[len_fit] log_lik_inc;
for (i in 1:len_fit){
log_lik_cpd[i] =
  normal_lpdf(logploss[fit_rec[i]]|mu_p[fit_rec[i]],sig_p[d[fit_rec[i]]]);
log_lik_inc[i] =
  normal_lpdf(logiloss[fit_rec[i]]|mu_i[fit_rec[i]],sig_i[d[fit_rec[i]]]);
}
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
  fitPOS=stan(model_code=POSmodel_stan, data = data.POS,init=init.POS,
              seed = 12345,iter=stan_iter,thin=stan_thin,
              chains = 4,pars=pars.list,
              control=list(adapt_delta=.9999,max_treedepth=50),
              refresh=0)
  fitPOS_summary=
    as.matrix(summary(fitPOS)$summary)[1:(W+4*D+4*C-1),c(1,3,10)]
  mrh=subset(fitPOS_summary,is.na(fitPOS_summary[,3])==F)
  max_Rhat=round(max(mrh[,3]),4)
  print(paste("Maximum Rhat =",max_Rhat,"Thin =",stan_thin))
  stan_thin=2*stan_thin
  stan_iter=2*stan_iter
}
stan_thin=stan_thin/2
print(fitPOS_summary)
#
# extract information from stan output to process in R
#
b <- extract(fitPOS)
logelr=b$logelr
alpha=b$alpha
beta_p=b$beta_p
beta_os=b$beta_os
gamma_p=b$gamma_pd
gamma_os=b$gamma_os
speedup_p=b$speedup_p
speedup_os=b$speedup_os
sig_p=b$sig_p
sig_i=b$sig_i
mu_p=b$mu_p
mu_i=b$mu_i
num.mcmc=length(logelr)
#
# goodness of fit statistics for comparing models
#
elpd_cpd=loo(b$log_lik_cpd)
elpd_inc=loo(b$log_lik_inc)
#
# calculate loss statistice by accident year for simulated outcomes
#
at.wdD_p=matrix(0,num.mcmc,W)
at.wdD_i=matrix(0,num.mcmc,W)
for (w in 1:W){
  mu.wdD=log(Premium[w])+alpha[,w]+logelr
  at.wdD_p[,w]=exp(mu.wdD+sig_p[,D]^2/2) # this the estimate given mu and sigma
  at.wdD_i[,w]=exp(mu.wdD+sig_i[,D]^2/2) # this the estimate given mu and sigma
}
#
# calculate loss statistics and output to data frame
#
ss.wdD_p=rep(0,W)
ms.wdD_p=rep(0,W)
ss.wdD_i=rep(0,W)
ms.wdD_i=rep(0,W)
#
for (w in 1:W){
  ms.wdD_p[w]=mean(at.wdD_p[,w])
  ss.wdD_p[w]=sd(at.wdD_p[,w])
  ms.wdD_i[w]=mean(at.wdD_i[,w])
  ss.wdD_i[w]=sd(at.wdD_i[,w])
}
Pred.POS_p=rowSums(at.wdD_p)
ms.tdD_p=mean(Pred.POS_p)
ss.tdD_p=sd(Pred.POS_p)
POS.Estimate_p=round(ms.wdD_p)
POS.SE_p=round(ss.wdD_p)
POS.CV_p=round(POS.SE_p/POS.Estimate_p,4)
POS.RetroRes_p=cumsum(POS.Estimate_p)
for (w in 1:W){
  POS.RetroRes_p[w]=POS.RetroRes_p[w]-sum(adata$ploss[adata$cal==w])
}
Pred.POS_i=rowSums(at.wdD_i)
ms.tdD_i=mean(Pred.POS_i)
ss.tdD_i=sd(Pred.POS_i)
POS.Estimate_i=round(ms.wdD_i)
POS.SE_i=round(ss.wdD_i)
POS.CV_i=round(POS.SE_i/POS.Estimate_i,4)
POS.RetroRes_i=cumsum(POS.Estimate_i)
for (w in 1:W){
  POS.RetroRes_i[w]=POS.RetroRes_i[w]-sum(adata$ploss[adata$cal==w])
}
# append totals
POS.Estimate_p=c(POS.Estimate_p,round(ms.tdD_p))
POS.SE_p=c(POS.SE_p,round(ss.tdD_p))
POS.CV_p=c(POS.CV_p,round(ss.tdD_p/ms.tdD_p,4))
POS.Estimate_i=c(POS.Estimate_i,round(ms.tdD_i))
POS.SE_i=c(POS.SE_i,round(ss.tdD_i))
POS.CV_i=c(POS.CV_i,round(ss.tdD_i/ms.tdD_i,4))
Premium=c(Premium,sum(Premium))
POS.RetroRes_p=c(POS.RetroRes_p,POS.RetroRes_p[W])
POS.RetroRes_i=c(POS.RetroRes_i,POS.RetroRes_i[W])
AY=c(sort(unique(adata$ay)),"Total")
AY_Exhibit=data.frame(AY,Premium,POS.Estimate_p,POS.Estimate_i,
                      POS.SE_p,POS.SE_i,POS.CV_p,POS.CV_i,
                      POS.RetroRes_p,POS.RetroRes_i)
print(AY_Exhibit)
Premium=Premium[W+1]
POS.Estimate_p=POS.Estimate_p[W+1]
POS.Estimate_i=POS.Estimate_i[W+1]
POS.SE_p=POS.SE_p[W+1]
POS.SE_i=POS.SE_i[W+1]
mean_elr=round(mean(exp(logelr)),3)
# warning -- the loo package output depends upon the version of the loo package
elpd_p=round(elpd_cpd$elpd_loo,2)
elpd_i=round(elpd_inc$elpd_loo,2)
SumStats=data.frame(Premium,POS.Estimate_p,POS.Estimate_i,POS.SE_p,POS.SE_i,
                    mean_elr,elpd_p,elpd_i,stan_thin)
print(SumStats)
#
# create histogram of reserve predictive distribution
#
par(mfrow=c(2,1))
reserve=sort(Pred.POS_p-sum(adata$ploss[adata$cal==max(keep_cal)]))
POS.ResMean_p=round(mean(reserve))
POS.ResLow_p=round(reserve[250])
POS.ResHigh_p=round(reserve[9750])
hist(reserve,
 main=paste("POS-vcp",run_id, "Predictive Distribution of the Loss Reserve"),
#     xlim=c(75000,1250000),breaks=49, #specific for paper
 xlab=paste("Mean =",POS.ResMean_p),
 sub=paste("Low End (2.5%) =",POS.ResLow_p," High End (97.5%) =",POS.ResHigh_p))
abline(v=POS.ResLow_p,lwd=3,col="red")
abline(v=POS.ResMean_p,lwd=3,col="blue")
abline(v=POS.ResHigh_p,lwd=3,col="red")
#
reserve=sort(Pred.POS_i-sum(adata$ploss[adata$cal==max(keep_cal)]))
POS.ResMean_i=round(mean(reserve))
POS.ResLow_i=round(reserve[250])
POS.ResHigh_i=round(reserve[9750])
hist(reserve,
     main=paste("POS-vci",run_id, "Predictive Distribution of the Loss Reserve"),
#     xlim=c(75000,1250000),breaks=49, #specific for paper
     xlab=paste("Mean =",POS.ResMean_i),
 sub=paste("Low End (2.5%) =",POS.ResLow_i," High End (97.5%) =",POS.ResHigh_i))
abline(v=POS.ResLow_i,lwd=3,col="red")
abline(v=POS.ResMean_i,lwd=3,col="blue")
abline(v=POS.ResHigh_i,lwd=3,col="red")
#
# standardized residual plots
#
w=adata$w
d=adata$d
Premium=Premium[1:W]
logloss_p=log(adata$ploss)
logloss_i=log(adata$iloss)
#
std_resid_p=NULL
std_resid_i=NULL
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
  for (i in 1:length(logloss_p)){
    if(i %in% fit_rec){
      std_resid_p=c(std_resid_p,(logloss_p[i]-mu_p[s,i])/sig_p[s,d[i]])
      std_resid_i=c(std_resid_i,(logloss_i[i]-mu_i[s,i])/sig_i[s,d[i]])
      ay=c(ay,w[i])
      dy=c(dy,d[i])
      cy=c(cy,w[i]+d[i]-1)
    }
  }
}
#rng=range(std_resid)
rng=c(-4.5,4.5)  #for the paper
par(mfrow=c(1,3))
#
# paid boxplots
#
temp=eq_s(std_resid_p,ay,s_size)
boxplot(split(temp$eq_r,temp$eq_y),notch=T,ylim=rng,
        col="gray",main="",xlab="Accident Year")
abline(h=0,lwd=5)
abline(h=qnorm(.25))
abline(h=qnorm(.75))
#
temp=eq_s(std_resid_p,dy,s_size)
boxplot(split(temp$eq_r,temp$eq_y),notch=T,ylim=rng,
        col="gray",main="",xlab="Development Year")
abline(h=0,lwd=5)
abline(h=qnorm(.25))
abline(h=qnorm(.75))
#
temp=eq_s(std_resid_p,cy,s_size)
boxplot(split(temp$eq_r,temp$eq_y),notch=T,ylim=rng,
        col="gray",main="",xlab="Calendar Year")
abline(h=0,lwd=5)
abline(h=qnorm(.25))
abline(h=qnorm(.75))
#
mtext(paste("POS-vc",run_id, "Standardized Residual Boxplots - Paid"),
      side = 3, line = -2.75, outer = TRUE)
#
# incurred boxplots
#
temp=eq_s(std_resid_i,ay,s_size)
boxplot(split(temp$eq_r,temp$eq_y),notch=T,ylim=rng,
        col="gray",main="",xlab="Accident Year")
abline(h=0,lwd=5)
abline(h=qnorm(.25))
abline(h=qnorm(.75))
#
temp=eq_s(std_resid_i,dy,s_size)
boxplot(split(temp$eq_r,temp$eq_y),notch=T,ylim=rng,
        col="gray",main="",xlab="Development Year")
abline(h=0,lwd=5)
abline(h=qnorm(.25))
abline(h=qnorm(.75))
#
temp=eq_s(std_resid_i,cy,s_size)
boxplot(split(temp$eq_r,temp$eq_y),notch=T,ylim=rng,
        col="gray",main="",xlab="Calendar Year")
abline(h=0,lwd=5)
abline(h=qnorm(.25))
abline(h=qnorm(.75))
#
mtext(paste("POS-vc",run_id, "Standardized Residual Boxplots - Incurred"),
      side = 3, line = -2.75, outer = TRUE)
#
# combined paid and incurred
#
temp=eq_s(c(std_resid_p,std_resid_i),c(ay,ay),s_size)
boxplot(split(temp$eq_r,temp$eq_y),notch=T,ylim=rng,
        col="gray",main="",xlab="Accident Year")
abline(h=0,lwd=5)
abline(h=qnorm(.25))
abline(h=qnorm(.75))
#
temp=eq_s(c(std_resid_p,std_resid_i),c(dy,dy),s_size)
boxplot(split(temp$eq_r,temp$eq_y),notch=T,ylim=rng,
        col="gray",main="",xlab="Development Year")
abline(h=0,lwd=5)
abline(h=qnorm(.25))
abline(h=qnorm(.75))
#
temp=eq_s(c(std_resid_p,std_resid_i),c(cy,cy),s_size)
boxplot(split(temp$eq_r,temp$eq_y),notch=T,ylim=rng,
        col="gray",main="",xlab="Calendar Year")
abline(h=0,lwd=5)
abline(h=qnorm(.25))
abline(h=qnorm(.75))
#
mtext(paste("POS-vc",run_id,
            "Standardized Residual Boxplots - Paid and Incurred"),
             side = 3, line = -2.75, outer = TRUE)
#
# plot paid and incurred gamma values
#
gamma_p_mean=colMeans((gamma_p[,keep_cal[-length(keep_cal)]]))
gamma_os_mean=colMeans((gamma_os[,keep_cal[-length(keep_cal)]]))
rng=range(gamma_p_mean,gamma_os_mean)
par(mfrow=c(1,1))
plot(keep_cal[-length(keep_cal)],gamma_p_mean,ylim=rng,type="l",
     lwd=3,col="blue",
     ylab=expression(paste("Mean(",gamma,")")),
     main="Mean Speedup Rates",
     xlab="Calendar Year")
par(new=T)
plot(keep_cal[-length(keep_cal)],gamma_os_mean,ylim=rng,type="l",
     lwd=3,col="black",
     main="",xlab="",ylab="")
abline(h=0)
legend(max(keep_cal)-3.5,max(rng),legend=c("Paid","Outstanding   "),
       lty=1:1,lwd=3,col=c("blue","black"))
#
# write output data files
#
parm_summary=round(fitPOS_summary[,1:2],4)
Runid=paste("POS-vcp",run_id)
Ultimate.Loss=POS.Estimate_p
Ultimate.SE=POS.SE_p
Reserve.Low=POS.ResLow_p
Reserve.Mean=POS.ResMean_p
Reserve.High=POS.ResHigh_p
ELPD=round(elpd_cpd$estimates[1,1],2)
res_stats=rbind(earlier_runs[,2:8],data.frame(Runid,Ultimate.Loss,Ultimate.SE,
                Reserve.Low,Reserve.Mean,Reserve.High,ELPD))
Ultimate.Loss=POS.Estimate_i
Ultimate.SE=POS.SE_i
Reserve.Low=POS.ResLow_i
Reserve.Mean=POS.ResMean_i
Reserve.High=POS.ResHigh_i
ELPD=round(elpd_inc$estimates[1,1],2)
res_stats=rbind(res_stats,data.frame(Runid,Ultimate.Loss,Ultimate.SE,
                Reserve.Low,Reserve.Mean,Reserve.High,ELPD))
if (is.null(res.stats)==F){
  write.csv(res_stats,file=res.stats)
}
if (is.null(AY.Exhibit)==F) write.csv(AY_Exhibit,file=AY.Exhibit)
if (is.null(parm.summary)==F) write.csv(parm_summary,file=parm.summary)
t1=Sys.time()
print(t1-t0)