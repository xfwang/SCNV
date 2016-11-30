A. File Preperation
-----------
**(1)** If it is a *.fastq* file, alignment need to be done using tools such as bowtie (not included in this toolkit), then go to step (2).

**(2)** If it is a cleaned *.bam* file, go to step (2d)

(2a) (Optional) Make sure to Sort and reorder .bam files (according to the referen genome) based on the script **bamGATKsort.sh** . The human reference genome files can be prepared using the script **hg19_reference.sh**  
Input:cell.bam Output: cell.final.bam
```
Usage: ./bamGATKsort.sh cell.bam  
```
(2b) (Optional) Run **DepthOfCoverage** (Input cell.final.bam)
```
java -jar GenomeAnalysisTK.jar \-omitBaseOutput \ -T DepthOfCoverage \ -R hg19.ucsc.fa \ -I cell.final.bam \ -o cell.coverage
```
(2c) (Optional) Draw coverage histogram and sample statistics

(2d)
