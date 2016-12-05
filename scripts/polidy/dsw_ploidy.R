##################
##################
dsw_plot<- function (rc.data,gc.data){
  #need to pre-load functions from gc_cbs_batch.R
  require(DNAcopy)
  for (i in 1:(ncol(rc_data)-4)) {
  #i: sample index
  # i=1
  coverage<-rc_data[,c(1,2,4,(i+4))]
  NormRatio<-gc.norm(badbins=NA, gc.data=gc.data, coverage=coverage)
  cell.id<-names(coverage)[4]

  #DSW plot
  pdf(paste(cell.id, ".dsw.pdf", sep=""))
    for (alpha in c(0.05, 0.0001)) {
        Seg<-CBS(NormRatio, sample.id=cell.id, alpha=alpha, nperm=1000, undo.SD=1.0, min.width=5)
        SegLong<-Seg[[1]]; SegShort<-Seg[[2]];

        #DSW (Damped sine wave)
        grid<- c(seq(1.5, 6.4, by=0.05))
        SCNP<- SegLong$seg.mean.LOWESS %o% grid
        FCNP <- round(SCNP)
        Cost <- (SCNP - FCNP) ^ 2
        SumCost <- colSums(Cost, na.rm = FALSE, dims = 1)
        min.grid <- round(grid[which.min(SumCost)])

        y=-log(SumCost,base=exp(1))
        y=y-mean(y)
        x=grid
          if (alpha==0.05) {
          plot(x, y,type="b",col="turquoise4",xlab="ploidy = argmax(-log(SumCostSquare))",ylab="-log(SumCostSquare)",main=cell.id)
          abline(v=min.grid,col="blue")
          #nfit<-nls(y~a*exp(-b*x)*(cos(f*x+phi)+sin(f*x+phi)),start=list(a=1, b=1, f=6, phi=10))
          #lines(x,predict(nfit,seq(1.6,6.4,by=0.05)),col="blue")
          } else {lines(x,y,col="grey")}
    }
  dev.off()

  ##
  pdf(paste(cell.id, ".mode.pdf", sep=""))
    #plot.baslan(SegLong, sample.id=cell.id)
    copy<-cne(df=SegLong, dfs=SegShort)
    #save(copy,file=paste(cell.id,".Rdata",sep=""))
  dev.off()
}
}
