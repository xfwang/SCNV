
## This function is to determine the thresholds we need to use for the batch of the data.
## Inputs: (1) the reads from normal cells used as control; (2) the reads from some tumor cells for calibration.
## In the files, the first column is the chromosome ID of the reads, i.e., chr1, chr2, ..., chr22, chrX; the second column is the locations of the reads.

NormalFile = "all.WG.seq"
TumorFiles = c("1185.WG.seq", "1190.WG.seq", "1197.WG.seq", "1200.WG.seq", "1210.WG.seq")


###############################################################################################################

## get raw estimates
library(seqCBS)
source("ScanCBS3.R")
control.all = read.table(NormalFile)
nControl = dim(control.all)[1]
Nsamp = length(TumorFiles)

result1 = c()
window = 1e4
control.windowCounts = list()
for (i in 1:Nsamp){
  case.all = read.table(TumorFiles[i])
  nCase = dim(case.all)[1]
  for (j in 1:22){
    case.id = which(case.all[,1]==paste("chr",j,sep=""))
    control.id = which(control.all[,1]==paste("chr",j,sep=""))
    controls = control.all[control.id,2]
    temp = ScanCBS(case.all[case.id,2],controls,p=nCase/(nCase+nControl))
    result1 = rbind(result1, cbind(rep(i,length(temp$tauHat)), rep(j,length(temp$tauHat)), temp$tauHat, c(temp$relCN,1)))
    if (i==1){
      maxVal = max(controls)
      control.w = ceiling(controls/window)
      control.count = rep(0,ceiling(maxVal/window))
      temp = table(control.w)
      control.count[as.numeric(names(temp))] = temp
      control.windowCounts[[j]] = control.count
    }
  }
}

colnames(result1) = c("sample", "chr", "tauhat", "relCN")

event.len = diff(result1[,3])
N = dim(result1)[1]
ids = which(event.len<0)
start.bp = result1[-c(ids,N),3]
end.bp = result1[-c(1,ids+1),3]
result2 = cbind(result1[-c(ids,N),1:2],start.bp, end.bp, event.len[-ids],result1[-c(ids,N),4])
colnames(result2) = c("sample", "chr", "start.bp", "end.bp", "event.len", "relCN")


## determine centromere region

centro = matrix(0,22,2)
for (j in 1:22){
  count = control.windowCounts[[j]]
  ids = which(count==0)
  dif = diff(ids)
  id2 = c(1,which(dif>1))
  if (length(id2)==1){
    id2 = c(1,length(dif))
  }
  dif2 = diff(id2)
  a = which.max(dif2)
  b1 = id2[a]
  b2 = id2[a+1]
  d1 = ids[b1]
  d2 = ids[b2]
  centro[j,] = c(d1,d2)
}


## select normal regions
N = dim(result2)[1]
relCN = result2[,6]
event.len = result2[,5]
ids = which(relCN>0.9 & relCN<1.1 & event.len>1e6)
normal = result2[ids,]

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


normal2 = c()
N = dim(normal)[1]
for (i in 1:N){
  temp = overlap(as.numeric(normal[i,3:4]), c(centro[normal[i,2],1]*window, (centro[normal[i,2],2]+1)*window))
  if (temp$nSeg==1){
    normal2 = rbind(normal2, c(normal[i,c(1:2,6)], temp$r))
  }else{
    normal2 = rbind(normal2, c(normal[i,c(1:2,6)], temp$r[1:2]), c(normal[i,c(1:2,6)], temp$r[3:4]))
  }
}

len = as.numeric(normal2[,5])-as.numeric(normal2[,4])
ids = which(len>1e7)
normal3 = cbind(normal2[ids,],len[ids])

## re-search in the normal region
statHat = c()
for (j in 1:22){
  id1 = which(normal3[,2]==j)
  if (length(id1)>0){
    control.id = which(control.all[,1]==paste("chr",j,sep=""))
    controls = control.all[control.id,2]
    normal.chr = normal3[id1,]
    for (i in 1:Nsamp){
      id2 = which(normal.chr[,1]==i)
      if (length(id2)>0){
        case.all = read.table(TumorFiles[i])
        nCase = length(case.all)
        case.id = which(case.all[,1]==paste("chr",j,sep=""))
        cases = case.all[case.id,2]
        for (k in 1:length(id2)){
          v = as.numeric(normal.chr[id2[k],4:5])
          m = floor((v[2]-v[1])/1e7)
          delta = floor((v[2]-v[1]-m*1e7)/2)
          for (ii in 1:m){
            a = v[1]+delta+(ii-1)*1e7
            b = a + 1e7
            control.id = which(controls>=a & controls<=b)
            case.id = which(cases>=a & cases<=b)
            case = cases[case.id]
            control = controls[control.id]
            r = ScanCBS(case, control, p=nCase/(nCase+nControl))
            
            case.w = ceiling(case/window)
            control.w = ceiling(control/window)
            
            case.count = control.count = rep(0,ceiling(max(c(case,control))/window))
            temp1 = table(case.w)
            temp2 = table(control.w)
            case.count[as.numeric(names(temp1))] = temp1
            control.count[as.numeric(names(temp2))] = temp2
            
            R = case.count/control.count * nControl/nCase
            tauhat.w = ceiling(r$tauHat/window)
            relCN = rep(1,length(case.count))
            for (k in 1:(length(tauhat.w)-1)){
              startid = tauhat.w[k]
              endid = tauhat.w[k+1]
              relCN[startid:endid] = r$relCN[k]
            }
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


