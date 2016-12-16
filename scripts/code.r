
## get raw estimates

load("composite.normal.bychr.rda")
load("cases.rda")

library(seqCBS)
dyn.load("seqCBS.so")
source("ScanCBS3.R")

ncontrol = 0
for (i in 1:22){
  ncontrol = ncontrol + length(composite.normal.bychr[[i]])
}

Nsamp = 4
statHat.all = c()
window = 1e4
for (i in 1:Nsamp){
  mycase = cases[[i]]
  ncase = dim(mycase)[1]
  for (j in 1:22){
    mycontrol = composite.normal.bychr[[j]]
    case.id = which(mycase[,1]==paste("chr",j,sep=""))
    temp = ScanCBS(mycase[case.id,2],mycontrol,p=ncase/(ncase+ncontrol))
    m = dim(temp$statHat)[1]
    statHat.all = rbind(statHat.all, cbind(rep(i,m),rep(j,m),temp$statHat))
  }
  print(i)
}
colnames(statHat.all) = c("sample", "chr", "cptL", "cptR", "cptLReadInd", "cptRReadInd", "stat", "parentL", "parentR", "BIC")

# save(statHat.all, file="statHat.all.Rdata")


# apply thresholds
load("cutoffs.rda")
len.bp = statHat.all[,4]-statHat.all[,3]
len.ind = statHat.all[,6]-statHat.all[,5]
stat = statHat.all[,7]

id1 = which(stat>stat.cut2*2)

id2 = which(len.ind>ind.cut0)
id2.1 = which(len.bp>len.cut2 & stat>stat.cut0)
id2.2 = which(len.bp>len.cut1 & stat>stat.cut1)
id2.3 = which(len.bp>len.cut0 & stat>stat.cut2)
id2.123 = unique(c(id2.1,id2.2,id2.3))

id2s = id2[which(!is.na(match(id2,id2.123)))]
statHat.remain1 = statHat.all[unique(sort(c(id1,id2s))),]

# get relCN
result3 = c()
for (j in 1:22){
  mycontrol = composite.normal.bychr[[j]]
  ids = which(statHat.remain1[,2]==j)
  statHat = statHat.remain1[ids,]
  conMin = min(mycontrol)
  conMax = max(mycontrol)
  for (i in 1:Nsamp){
    id2 = which(statHat[,1]==i)
    mycase = cases[[i]]
    ncase = dim(mycase)[1]
    case.id = which(mycase[,1]==paste("chr",j,sep=""))
    caseMin = min(mycase[case.id,2])
    caseMax = max(mycase[case.id,2])
    tauhat = sort(unique(c(min(conMin,caseMin), max(conMax,caseMax),statHat[id2,3],statHat[id2,4])))
    for (k in 1:(length(tauhat)-1)){
      a = tauhat[k]
      b = tauhat[k+1]-1
      if (k==(length(tauhat)-1)) b = tauhat[k+1]
      p = length(which(mycase[case.id,2]>=a & mycase[case.id,2]<=b))/length(which(mycontrol>=a & mycontrol<=b)) * ncontrol/ncase
      result3 = rbind(result3, c(i,j,tauhat[k],tauhat[k+1],p))
    }
  }
}

result3[,5] = 2*result3[,5]
colnames(result3) = c("sample", "chr", "start.bp", "end.bp", "CN")
write.table(result3, file="raw_CN.txt", quote=F, row.names=F)



