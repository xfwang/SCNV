SinSeq: Single-cell sequencing analysis toolkit
=======

The Single-cell Sequencing Analysis Toolkit (or SinSeq) provides various utilities for manipulating and analyzing data generated from Single-cell sequencing.


1. Preprocessing
-----------


Steps:    
0) .fastq --> .bam (alignment using tools such as bowtie, not included in this toolkit)    
1) Sort and reorder .bam files (according to the referen genome) based on the script **bamGATKsort.sh** . The human reference genome files can be prepared using the script **hg19_reference.sh**
```
Usage: ./bamGATKsort.sh cellname.bam  
```
2) Run DepthOfCoverage
```
java -jar GenomeAnalysisTK.jar \-omitBaseOutput \ -T DepthOfCoverage \ -R hg19.ucsc.fa \ -I CJA0748.final.bam \ -o CJA0748.coverage
```





2. Initial CNV discovery (control-free calling)
-----------



3. Control-based CNV-calling 
-----------




4. Intratumor heterogeneity (clustering)
-----------
