A. File Preperation
-----------
**(1)** If it is a *.fastq* file, alignment need to be done using tools such as bowtie (not included in this toolkit), then go to step (2).

**(2)** If it is a well cleaned *.bam* file, go directly to step (2d)

- (2a) (Optional) Make sure to sort and reorder .bam files (according to the referen genome) based on the script **bamGATKsort.sh** . The human reference genome file is required and can be prepared using the script **hg19_reference.sh**  
  Input:cell.bam Output: cell.final.bam
  ```
  Usage: ./bamGATKsort.sh cellname.bam  
  ```

- (2b) (Optional) Run **DepthOfCoverage** (Input cell.final.bam)
  ```
java -jar GenomeAnalysisTK.jar \-omitBaseOutput \ -T DepthOfCoverage \ -R hg19.ucsc.fa \ -I cell.final.bam \ -o cellname.coverage
  ```

- (2c) (Optional) Draw coverage histogram and sample statistics

- (2d) Convert .bam file to .bed file using Bedtools (make sure [Bedtools] (http://bedtools.readthedocs.io/en/latest/content/installation.html) is installed  )
  ```
bamToBed -i cellname.bam > cellname.bed
  ```
  
B. Initial CNV discovery (Bin-based segmentation)
-----------
**One can use online tool [Ginkgo] (http://qb.cshl.edu/ginkgo/) for intial profiling (need to upload all .bed files, minimum 3 cell files) or follow following steps for running a large number of cells on local computers:** 

**(1)** Define boundaries based on mappable positions and calculate the GC content in each bin. This step can be done using the scripts ``hg19.bin.bondaries.50k.py`` and ``hg19.varbin.gc.content.50k.bowtie.k50.py``(Baslan Nat. Protoc. 2012).
Refer to **hg19_reference.sh** for the reference genome.

**(2)**  Count the number of reads in each defined bin in each file (batch mode) using the codes implented in **read_count.sh** and **read_count.Rscript** . A plot can be very helpful here by comparing the distribution of counts per bins between the cell samples as  implemented in **rc_plot.R**. The plot function also offers options to convert raw read count to "RPKM" and with median normalization. 

**(3)** GC correction and Initial CBS Segmentation (use background read depth as control): apply wrapper R functions in **gc_cbs.R**


C. Binless Segmentation
-----------


D. Example ([SRA] (https://www.ncbi.nlm.nih.gov/sra) cell file [SRR1548983](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1548983))
-----------






