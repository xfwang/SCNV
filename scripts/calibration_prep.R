
#### calibration for the data

load("composite.normal.bychr.rda")
load("calibration.normal.rda")

library(seqCBS)
dyn.load("seqCBS.so")
source("ScanCBS3.R")

ncontrol = dim(composite.normal)[1]
ncali = rep(0,length(calibration.normal))
for (i in 1:length(calibration.normal)){
  ncali[i] = dim(calibration.normal[[i]])[1]
}


## initial CNV calling

result = c()
for (i in c(1:22)){
  control = composite.normal.bychr[[i]]
  for (j in 1:length(calibration.normal)){
    cali = calibration.normal[[j]]
    case = cali[which(cali[,1]==paste("chr",i,sep="")),2]
    temp = ScanCBS(case,control,p=ncali[j]/(ncali[j]+ncontrol))
    result = rbind(result, cbind(rep(i,length(temp$tauHat)), rep(j,length(temp$tauHat)), temp$tauHat, c(temp$relCN,1)))
  }
}

colnames(result) = c("chr", "sample", "tauhat", "relCN")
write.table(result, file="cali_CNV_result.txt", quote=F, row.names=F)


## determine centromere region

window = 1e4
windowCounts = c()
for (i in c(1:22)){
  control = composite.normal.bychr[[i]]
  control.w = ceiling(control/window)
  temp = table(control.w)
  count = rep(0,ceiling(max(control)/window))
  count[as.numeric(names(temp))] = temp
  windowCounts = rbind(windowCounts, cbind(rep(i,length(count)),count))
}

colnames(windowCounts) = c("chr", "count")
write.table(windowCounts, file="composite.normal.windowCounts.txt", quote=F, row.names=F)


centro = matrix(0,22,2)
for (i in 1:22){
  id = which(windowCounts[,1]==i)
  count = windowCounts[id,2]
  id1 = which(count==0)
  dif = diff(id1)
  id2 = c(1,which(dif>1),length(dif))
  dif2 = diff(id2)
  a = which.max(dif2)
  b1 = id2[a]
  b2 = id2[a+1]
  d1 = id1[b1+1]-1
  d2 = id1[b2]
  centro[i,] = c(d1,d2)
}

save(centro, file="centro.rda")


## remove inherited CNV and centro regions

result = read.table("cali_CNV_result.txt", header=T)
event.len = diff(result$tauhat)
N = dim(result)[1]
ids = which(event.len<0)
start.bp = result$tauhat[-c(ids,N)]
end.bp = result$tauhat[-c(1,ids+1)]
result2 = cbind(result[-c(ids,N),1:2],start.bp, end.bp, event.len[-ids],result[-c(ids,N),4])
colnames(result2) = c("chr", "sample", "start.bp", "end.bp", "event.len", "relCN")
write.table(result2, file="cali_CNV_result2.txt", quote=F, row.names=F)

# hist(result2$relCN*2,xlim=c(0,4),breaks=100000)

N = dim(result2)[1]
relCN = result2[,6]
event.len = result2[,5]
ids = which(relCN>0.9 & relCN<1.1 & event.len>1e6)
remain0 = result2[ids,]

# remove overlap
overlap = function(v1, v2){
  if (v1[2]<v2[1] || v1[1]>v2[2]){
    nSeg = 1
    r = v1
  }else if (v1[1]<v2[1] && v1[2]>v2[2]){
    nSeg = 2
    r = c(v1[1],v2[1]-1, v2[2]+1, v1[2])
  }else if (v1[1]>v2[1]){
    nSeg = 1
    r = c(v2[2],v1[2])
  }else{
    nSeg = 1
    r = c(v1[1],v2[1])
  }
  return(list(nSeg=nSeg, r=r))
}

remain = c()
N = dim(remain0)[1]
for (i in 1:N){
  mychr = match(remain0[i,1],1:22)
  temp = overlap(as.numeric(remain0[i,3:4]), c(centro[mychr,1]*window, (centro[mychr,2]+1)*window))
  if (temp$nSeg==1){
    remain = rbind(remain, c(remain0[i,2], mychr, remain0[i,6], temp$r))
  }else{
    remain = rbind(remain, c(remain0[i,2], mychr, remain0[i,6], temp$r[1:2]), c(remain0[i,2], mychr, remain0[i,6], temp$r[3:4]))
  }
}


len = as.numeric(remain[,5])-as.numeric(remain[,4])
minLen = 2e6
ids = which(len>minLen)
remain1 = cbind(remain[ids,],len[ids])

## re-search in the remaining region
statHat = c()
for (j in 1:22){
  print(j)
  id1 = which(remain1[,2]==j)
  if (length(id1)>0){
    control = composite.normal.bychr[[j]]
    normal.chr = remain1[id1,]
    for (i in 1:length(calibration.normal)){
      id2 = which(normal.chr[,1]==i)
      if (length(id2)>0){
        cali = calibration.normal[[i]]
        case = cali[which(cali[,1]==paste("chr",i,sep="")),2]
        for (k in 1:length(id2)){
          v = as.numeric(normal.chr[id2[k],4:5])
          m = floor((v[2]-v[1])/minLen)
          delta = floor((v[2]-v[1]-m*minLen)/2)
          for (ii in 1:m){
            a = v[1]+delta+(ii-1)*minLen
            b = a + minLen
            control.id2 = which(control>=a & control<=b)
            case.id2 = which(case>=a & case<=b)
            mycase2 = case[case.id2]
            mycontrol2 = control[control.id2]
            r = ScanCBS(mycase2, mycontrol2, p=ncali[i]/(ncali[i]+ncontrol))
            
            stat = r$statHat
            if (length(stat)>0){
              statHat = rbind(statHat, cbind(rep(j,dim(stat)[1]), rep(i,dim(stat)[1]), stat))
            }
          }
        }
      }
    }
  }
}

colnames(statHat) = c("chr", "sample", "cptL", "cptR", "cptLReadInd", "cptRReadInd", "stat", "parentL", "parentR", "BIC")
a = statHat[,4]-statHat[,3]
b = statHat[,6]-statHat[,5]
stat = statHat[,7]

stat.cut0 = quantile(stat,0.5)
stat.cut1 = quantile(stat,0.995)
stat.cut2 = max(stat)

len.cut0 = quantile(a,0.5)
len.cut1 = quantile(a,0.995)
len.cut2 = max(a)

ind.cut0 = 1000

save(stat.cut0,stat.cut1,stat.cut2,len.cut0,len.cut1,len.cut2,ind.cut0, file="cutoffs.rda")
