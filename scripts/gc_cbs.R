##Initial segmenation and copy number estimation (CBS--Circular binary Segmentation)
##The following fucntions are wrapper functions modified based on Baslan 2012 Nat Protol-- SRR054616.cbs.r & SRR054616.copynumber.r


##GC normalizaiton
gc.norm<-function (badbins=NA, gc.data, coverage) {
	#bad.binsF
	#gc.data  GC content File: Refere to Baslan 2012 for the format
	#coverage Coverage File:
	#badbins=NA; gc.data="hg19.varbin.gc.content.5k.bowtie.k50.txt"; coverage="cell.bin.rc"

	lowess.gc <- function(jtkx, jtky) {
			jtklow <- lowess(jtkx, log(jtky), f=0.05)
			jtkz <- approx(jtklow$x, jtklow$y, jtkx)
			return(exp(log(jtky) - jtkz$y))
	}

	gc<-read.table(gc.data,header=T)
	if (is.na(badbins)){bad<-0 } else {bad <-read.table(bad.bins, header=F)}
	
	thisRatio <- read.table(coverage, header=T)[,c("chr","start","abs.pos","coverage")]
	names(thisRatio) <- c("chrom", "chrompos", "abspos", "bincount")
	a <- thisRatio$bincount + 1  #+1?
	thisRatio$ratio <- a / mean(a)
    
	if (all.equal(gc$bin.start.abspos, thisRatio$abspos)!=T){warning("files not matched") }
	thisRatio$gc.content <- gc$gc.content
	thisRatio$lowratio <- lowess.gc(thisRatio$gc.content, thisRatio$ratio)
	#plot(thisRatio$gc.content,thisRatio$lowratio)
	#plot(thisRatio$gc.content,thisRatio$ratio)
  return(thisRatio)
}

##CBS wrapper
CBS<-function(thisRatio, sample.id="cell", alpha=0.05, nperm=1000, undo.SD=1.0, min.width=5) {
  #sample.id="cell";alpha=0.05; nperm=1000; undo.SD=1.0; min.width=5
  require("DNAcopy")
  	CNA.object <- CNA(log(thisRatio$lowratio, base=2), thisRatio$chrom, thisRatio$chrompos, 
		data.type="logratio", sampleid=sample.id)
	smoothed.CNA.object <- smooth.CNA(CNA.object) 
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=min.width) 
	thisShort <- segment.smoothed.CNA.object[[2]]
	m <- matrix(data=0, nrow=nrow(thisRatio), ncol=1)	
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
		prevEnd = thisEnd
	}
	thisRatio$seg.mean.LOWESS <- m[, 1]
  OUT<-list(thisRatio,thisShort)
  return(OUT)
}

##Profile plots
plot.baslan<-function (thisRatio, sample.id="cell") {
	chr <- thisRatio$chrom
	chr.shift <- c(chr[-1], chr[length(chr)])
	vlines <- c(1, thisRatio$abspos[which(chr != chr.shift) + 1], thisRatio$abspos[nrow(thisRatio)])
	hlines <- c(0.5, 1.0, 1.5, 2.0)
	chr.text <- c(1:22, "X", "Y")
	chr.text <- chr.text[c(1:14, 16, 18, 20, 22:24)]
	vlines.shift <- c(vlines[-1], 4*10^9)
	chr.at <- vlines + (vlines.shift - vlines) / 2
	chr.at <- chr.at[c(1:14, 16, 18, 20, 22:24)]
	x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
	x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")
	y.at <- c(0.005, 0.020, 0.100, 0.500, 2.000)
	y.labels <- c("0.005", "0.020", "0.100", "0.500", "2.000")

	#png(paste(sample.id, ".wg.png", sep=""), height=400, width=600)
	par(mar=c(5.1,4.1,4.1,4.1))
	plot(x=thisRatio$abspos, y=thisRatio$lowratio, log="y", main=paste(sample.id), xaxt="n", xlab="Genome Position Gb", yaxt="n", ylab="Ratio", col="#CCCCCC", cex=0.5)
	axis(1, at=x.at, labels=x.labels)
	axis(2, at=y.at, labels=y.labels)
	lines(x=thisRatio$abspos, y=thisRatio$lowratio, col="#CCCCCC", cex=0.5)
	points(x=thisRatio$abspos, y=thisRatio$seg.mean.LOWESS, col="turquoise4", cex=0.5)
	lines(x=thisRatio$abspos, y=thisRatio$seg.mean.LOWESS, col="turquoise4", cex=0.5)
	abline(h=hlines)
	abline(v=vlines)
	mtext(chr.text, at = chr.at)
	#dev.off()  
}


##copy number estimation
##Caution: not working for diploid cells
cne<- function (df,dfs) {
	#df=SegLong; dfs=SegShort

	peaks <- function(series, span=3, ties.method = "first") {
		###  This peaks function from Petr Pikal https://stat.ethz.ch/pipermail/r-help/2007-February/125241.html
	if((span <- as.integer(span)) %% 2 != 1) stop("'span' must be odd")
	z <- embed(series, span)
	s <- span%/%2
	v <- max.col(z, ties.method=ties.method) == 1 + s
	pad <- rep(FALSE, s)
	result <- c(pad, v, pad)
	result
	}

	starts=ends=NULL
	prevEnd <- 0
	len <- nrow(dfs)
	for (j in 1:len) {
		thisStart = prevEnd + 1
		thisEnd = thisStart + dfs$num.mark[j] - 1
		starts <- c(starts, thisStart)
		ends <- c(ends, thisEnd)
		prevEnd = thisEnd
	} 
	
	##well, basically pairwise differences between seg.mean.LOWESS...
	amat <- matrix(data=0, nrow=1500000, ncol=1)
	counter <- 1
	for (j in 1:(len-1)) {
		for (k in (j+1):len) {
			N <- round((starts[j] - ends[j] + 1) * (starts[k] - ends[k] + 1)/1000)
			D <- abs(2^dfs$seg.mean[j] - 2^dfs$seg.mean[k])
			#cat(N, "\t")
			if (N > 0) {
				amat[(counter:(counter+N-1)), 1] <- rep.int(D, N)
				counter <- counter+N
			}
		}
	}

	a3 <- amat[(1:counter),1]
	#a3.95 <- sort(a3)[round(.95*counter)]
	a3.95 <- quantile(a3,0.95)
	a3d <- density(a3[which(a3 < a3.95)], n=1000)
	cn0 <- a3d$x[which(peaks(as.vector(a3d$y), span=59))][1]
	cn1 <- a3d$x[which(peaks(as.vector(a3d$y), span=59))][2]

	df$cn.ratio <- df$lowratio / cn1
	df$cn.seg <- df$seg.mean.LOWESS / cn1
	df$copy.number <- round(df$cn.seg)

	par(mar=c(5.1,4.1,4.1,4.1))
	plot(a3d, main="seg.mean difference density")

  return(df)
}


####
####
NormRatio<-gc.norm(badbins=NA, gc.data="hg19.varbin.gc.content.5k.bowtie.k50.txt", coverage="cell.bin.rc")
Seg<-CBS(NormRatio, sample.id="cell", alpha=0.05, nperm=1000, undo.SD=1.0, min.width=5)
SegLong<-Seg[[1]]; SegShort<-Seg[[2]]; 
plot.baslan(SegLong)
