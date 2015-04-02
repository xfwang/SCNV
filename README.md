Sinsat
=======


xuefeng.wang@stonybrook.edu 

hxchen@ucdavis.edu


Sinsat provides various utilities for manipulating and analyzing data generated from Single-cell DNA sequencing. Current pipeline majorly facilitates the binless segemenation (CNA break point detection) of single-cell DNA sequencing data. 


1. Preprocessing
-----------
System Requirements: R,GATK,samtools,picard, bedtools(optional)

**Steps:**    
0) .fastq --> .bam (alignment using tools such as bowtie, not included in this toolkit)    
1) Sort and reorder .bam files (according to the referen genome) based on the script **bamGATKsort.sh** . The human reference genome files can be prepared using the script **hg19_reference.sh**  
Input:cell.bam Output: cell.final.bam
```
Usage: ./bamGATKsort.sh cell.bam  
```
2) (Optional) Run **DepthOfCoverage** (Input cell.final.bam)
```
java -jar GenomeAnalysisTK.jar \-omitBaseOutput \ -T DepthOfCoverage \ -R hg19.ucsc.fa \ -I cell.final.bam \ -o cell.coverage
```
3) (Optional) Draw coverage histogram and sample statistics



2. Initial CNV discovery
-----------
**Steps:** 

1) Define boundaries based on mappable positions and calculate the GC content in each bin. This step can be done using the scripts ``hg19.bin.bondaries.50k.py`` and ``(hg19.varbin.gc.content.50k.bowtie.k50.py``( Baslan Nat. Protoc. 2012).
Refer to **hg19_reference.sh** for the reference genome.

2)  Count the number of reads in each defined bin in each file (batch mode) using the codes implented in **read_count.sh** and **read_count.Rscript** . A plot can be very helpful here by comparing the distribution of counts per bins between the cell samples as  implemented in **rc_plot.R**. The plot function also offers options to convert raw read count to "RPKM" and with median normalization. 

3) GC correction and Initial CBS Segmentation (use background read depth as control): **gc_cbs.R**


3. Segmentation based on nonhomogeneous poisson process 
-----------


