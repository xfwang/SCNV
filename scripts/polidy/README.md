
Example
--------

  ```
  rc_data<-read.table("/data/rc_matrix.txt",head=TRUE,colClasses = "numeric")
  gc.data<-read.table("/data/hg19.varbin.gc.content.50k.bowtie.k50.txt",header=T)
  dsw_plot(rc_data,gc.data)
  ```
