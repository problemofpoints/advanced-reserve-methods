#
# Script to run the Mack Model on data from the CAS Loss Reserve Database
# by Glenn Meyers
# Run time on my computer - Nearly instantaneous
#
#
rm(list = ls())  		# clear workspace
t0=Sys.time()
library(data.table)
#
# function to get Schedule P triangle data given ins group and line of business
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
                      cum_pdloss,cum_incloss,bulk_loss,inc_pdloss,single,posted_reserve97)
  data.out=data.out[order(d,w)]
  return(data.out)
}
#
# set up function to run Mack model and create output
#
model_function=function(CASdata,grpcode,losstype){
  DataG <- ins.line.data(CASdata,grpcode)
  train_data=subset(DataG,DataG$cal<11)
  Premium=train_data$net_premium[1:10]
  if (losstype=="Incurred") loss=train_data$cum_incloss
  if (losstype=="Paid") loss=train_data$cum_pdloss
  w=train_data$w
  d=train_data$d
  rectangle=data.table(loss,w,d)
  #
  # run Mack model using ChainLadder
  #
  library(ChainLadder)
  triangle=as.triangle(rectangle,origin="w",dev="d",value="loss")
  mcl=MackChainLadder(triangle,est.sigma="Mack")
  #
  # Calculate summary statistics for Mack model
  #
  # first by accident year
  W=1:10
  Mack.Estimate=round(summary(mcl)$ByOrigin[,3])
  Mack.SE=round(summary(mcl)$ByOrigin[,5])
  Mack.CV=round(Mack.SE/Mack.Estimate,4)
  if (losstype=="Incurred") aloss=DataG$cum_incloss
  if (losstype=="Paid") aloss=DataG$cum_pdloss
  Actual=subset(aloss,DataG$d==10)
  Pct.Mack=rep(NA,10)
  risk=data.table(W,Mack.Estimate,Mack.SE,Mack.CV,Actual,Pct.Mack)
  #
  W="Total"
  Mack.Estimate=round(sum(summary(mcl)$ByOrigin[1:10,3]))
  Mack.SE=round(summary(mcl)$Totals[5,1])
  Mack.CV=round(Mack.SE/Mack.Estimate,4)
  Mack.sig=sqrt(log(1+Mack.CV^2))
  Mack.mu=log(Mack.Estimate)-Mack.sig^2/2
  Pct.Mack=round(plnorm(sum(Actual),Mack.mu,Mack.sig)*100,2)
  Actual=sum(Actual)
  Group=grpcode
  risk.for.total=data.table(W,Mack.Estimate,Mack.SE,Mack.CV,Actual,Pct.Mack)
  risk=rbind(risk,risk.for.total)
  SumStats=data.table(Group,Mack.Estimate,Mack.SE,Mack.CV,Actual,Pct.Mack)
  output=list(risk=risk,
              SumStats=SumStats,
              triangle=triangle)
  return(output)
}
#
# Single triangle

CASdata = read.csv("~/Dropbox/CAS Loss Reserve Database/comauto_pos.csv")
# Location of files in the CAS Loss Reserve Database
#    http://www.casact.org/research/reserve_data/comauto_pos.csv
#    http://www.casact.org/research/reserve_data/ppauto_pos.csv
#    http://www.casact.org/research/reserve_data/wkcomp_pos.csv
#    http://www.casact.org/research/reserve_data/othliab_pos.csv
linobus="CA"
grpcode=353
losstype="Paid"
co_model=model_function(CASdata,grpcode,losstype)
print(" ")
print(co_model$risk)
print("")
print(co_model$SumStats)
plot(co_model$triangle)
