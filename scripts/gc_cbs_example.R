
####Example##
source(gc_cbs_functions.R)

rc_data<-read.table("/ycga-ba/home/xw322/knouse/SCNV/rc_matrix.txt",head=TRUE,colClasses = "numeric")
gc.data<-read.table("/ycga-ba/home/xw322/knouse/SCNV/hg19.varbin.gc.content.50k.bowtie.k50.txt",header=T)  #gc file
setwd("/ycga-ba/home/xw322/knouse/SCNV/output/") #results directory
require(DNAcopy)

for (i in 1:(ncol(rc_data)-4)) {
	#i: sample index
	# i=1
	coverage<-rc_data[,c(1,2,4,(i+4))]
	NormRatio<-gc.norm(badbins=NA, gc.data=gc.data, coverage=coverage)
	cell.id<-names(coverage)[4]
	Seg<-CBS(NormRatio, sample.id=cell.id, alpha=0.05, nperm=1000, undo.SD=1.0, min.width=5)
	SegLong<-Seg[[1]]; SegShort<-Seg[[2]];

	pdf(paste(cell.id, ".pdf", sep=""),width=30, height=8)
	layout(matrix(1:2, 1, 2,byrow=T),widths=c(2,1))
		plot.baslan(SegLong, sample.id=cell.id)
		copy<-cne(df=SegLong, dfs=SegShort)
		save(copy,file=paste(cell.id,".Rdata",sep=""))
	dev.off()
}
