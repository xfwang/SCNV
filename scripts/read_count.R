
##
INFILE_RC="cell.rc"   #readcount file
INFILE_BO="hg19.bin.boundaries.5k.bowtie.k50.txt"  #boundary file
OUTFILE=""

RC<-read.table(INFILE_RC, header=F,as.is=T); names(RC)<-c("chr","pos")
BO<-read.table(INFILE_BO, header=F,as.is=T)[,c(1,2,4)]; names(BO)<-c("chr","start","end")

chrcode<-function(x) {
	if (class(x)!="numeric") {
	if (length(x)==1) {
		x<-strsplit(x,"r")[[1]][2]	
		if (x=="X") {x=23}; if (x=="Y") {x=24}; if (x=="M") {x=25}
		} else if (length(x)>1) {
		x=do.call(rbind,strsplit(x,"r"))[,2]
		x[x=="X"]=23; x[x=="Y"]=24; x[x=="M"]=25
	}}	
	return(as.numeric(x))
}

RC$chr<-chrcode(RC$chr)
BO$chr<-chrcode(BO$chr)

RC1<-split(RC,RC$chr)
coverage<-apply(BO,1,function(arow){
	sum((arow[[2]]<RC1[[arow[[1]]]]$pos)*(RC1[[arow[[1]]]]$pos<arow[[3]]))
})

#write.table(coverage,file=OUTFILE)
