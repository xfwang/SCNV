#Exploration analysis of read count matrix file
#xuefeng.wang@stonybrook.edu


#Violin plot function to compare the distribution of counts per bins between the cell samples
violin<-function (x,scale="RPKM",norm=TRUE,fname="rc_explore") {
  #x=rc_matrix[,1:10];scale="RAW"; norm=TRUE
  require(vioplot)
  
  rc<-as.matrix(x[,c(-1:-4)])
  names<-colnames(rc)
  if (scale=="RPKM"){
	  length<-(x$end-x$start)/1000000 #length 
	  km<-colSums(rc)/1000000 #total reads
      rc<- sweep(rc,2,km,`/`)
	  rc<- sweep(rc,1, length, `/`)
	  #QC
	  rc[rc>quantile(rc,0.99)]=NA
  } else if (scale=="RAW") {
      rc[rc>quantile(rc,0.99)]=NA
  }
  
  if (norm==TRUE){
	  rc.list<-lapply(as.list(as.data.frame(rc)),function(x){x<-x[!is.na(x)]/median(x,na.rm=TRUE)})
  } else {
	   rc.list<-lapply(as.list(as.data.frame(rc)),function(x){x<-x[!is.na(x)]})
  }
    names(rc.list)[1]<-"x"
  
  postscript(paste(fname, ".ps",sep=""))
	do.call(vioplot,c(rc.list,list(names=names,horizontal=TRUE,col="lightblue")))
  dev.off()
}

rc_data<-read.table("/hms/scratch1/wang/work/rc_matrix.txt",head=TRUE,colClasses = "numeric")
rc_matrix<-rc_data[,1:14]
violin(rc_matrix,scale="RAW",norm=FALSE,fname="rc_explore_RAW")