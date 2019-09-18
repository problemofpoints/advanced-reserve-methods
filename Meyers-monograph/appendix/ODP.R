#
# Script to run the ODP Model on data from the CAS Loss Reserve Database
# by Glenn Meyers
# Run time on my computer - < 10 seconds
#
rm(list = ls())  		# clear workspace
t0=Sys.time()
library(data.table)
#
# function to get Schedule P triangle data given ins group and line of business
#
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
model_function=function(CASdata,grpcode){
  set.seed(12345)
  DataG <- ins.line.data(CASdata,grpcode)
  train_data=subset(DataG,DataG$cal<11)
  Premium=train_data$net_premium[1:10]
  loss=train_data$cum_pdloss
  w=train_data$w
  d=train_data$d
  rectangle=data.table(loss,w,d)
  #
  # run Mack model using ChainLadder
  #
  library(ChainLadder)
  numsim=10000
  triangle=as.triangle(rectangle,origin="w",dev="d",value="loss")
  odp=BootChainLadder(triangle, R = numsim, process.distr="od.pois")
  #
  # Calculate summary statistics for ODP model
  #
  ODP.Estimate=round(summary(odp)$ByOrigin[,2])
  ODP.S.E=round(summary(odp)$ByOrigin[,4])
  ODP.CV=round(ODP.S.E/ODP.Estimate,4)
  ODP.total.ult=round(summary(odp)$Totals[2,1])
  ODP.total.se=round(summary(odp)$Totals[4,1])
  ODP.total.cv=round(ODP.total.se/ODP.total.ult,4)
  act10d=subset(DataG$cum_pdloss,DataG$d==10)
  acttot=sum(act10d)
  le.act=(odp$IBNR.Totals+summary(odp)$Totals[1,1])<=acttot
  pct.ODP=round(sum(le.act)/numsim*100,2)
  W=c(1:10,"Total")
  ODP.Estimate=c(ODP.Estimate,ODP.total.ult)
  ODP.SE=c(ODP.S.E,ODP.total.se)
  ODP.CV=c(ODP.CV,ODP.total.cv)
  Actual=c(act10d,acttot)
  Percentile=c(rep("",10),pct.ODP)
  risk=data.frame(W,ODP.Estimate,ODP.SE,ODP.CV,Actual,Percentile)
  Group=grpcode
  ODP.Estimate=ODP.Estimate[11]
  ODP.SE=ODP.SE[11]
  Actual=Actual[11]
  ODP.Percentile=Percentile[11]
  SumStats=data.frame(Group,ODP.Estimate,ODP.SE,Actual,ODP.Percentile)
  output=list(risk=risk,
              SumStats=SumStats,
              triangle=triangle)
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
linobus="CA"
grpcode=353
co_model=model_function(CASdata,grpcode)
print(" ")
print(co_model$risk)
print("")
print(co_model$SumStats)
t1=Sys.time()
print(t1-t0)

