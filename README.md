SinSeq: Single-cell sequencing analysis toolkit
=======

The Single-cell Sequencing Analysis Toolkit (or SinSeq) provides various utilities for manipulating and analyzing data generated from Single-cell sequencing.

System Requirements: R,GATK,samtools,picard




1. Preprocessing
-----------


**Steps:**    
0) .fastq --> .bam (alignment using tools such as bowtie, not included in this toolkit)    
1) Sort and reorder .bam files (according to the referen genome) based on the script **bamGATKsort.sh** . The human reference genome files can be prepared using the script **hg19_reference.sh**  
Input:cell.bam Output: cell.final.bam
```
Usage: ./bamGATKsort.sh cell.bam  
```
2) Run **DepthOfCoverage** (Input cell.final.bam)
```
java -jar GenomeAnalysisTK.jar \-omitBaseOutput \ -T DepthOfCoverage \ -R hg19.ucsc.fa \ -I cell.final.bam \ -o cell.coverage
```
3) Draw coverage histogram and sample statistics




2. Initial CNV discovery (segmentation)
-----------
**Steps:** 
1) Convert the .bam file to read-coverage (.rc) file using the following one-liner:   

```
samtools view cell.final.bam | awk '{print $3,$4}' > cell.rc
```



3. CNV-calling with controls
-----------




4. Intra-tumor heterogeneity
-----------
