
Example
--------

  ```
  rc_data<-read.table("/data/rc_matrix.txt",head=TRUE,colClasses = "numeric")
  gc.data<-read.table("/data/bin50k.txt",header=T)
  dsw_plot(rc_data,gc.data)
  ```
