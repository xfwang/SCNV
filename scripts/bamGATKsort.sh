#!/bin/bash

#bam file preprocessing -->GATK-ready
##xuefeng.wang@moffitt.org

date #echo the time at start
picard_path="~/libraries/picard-tools-1.107" #picard directory
ref_fa="hg19.ucsc.fa"   #referene .fasta file name  (also need .fai and .dict files)

###bam file preprocessing
#split and merge and modify header (split.sh)
for SAMPLE in "@1"   #.bam file names
    do
		##subset (remove 'random' pieces)
        samtools index ${SAMPLE}.bam #index
        for C in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M
        do
        samtools view ${SAMPLE}.bam chr${C} -b > ${SAMPLE}.${C}.chr.bam
        done
        samtools merge ${SAMPLE}.merged.bam  *.chr.bam
        rm *.chr.bam
		samtools view -h ${SAMPLE}.merged.bam | sed '/gl00/d' > ${SAMPLE}.merged.sam
		samtools view -bS ${SAMPLE}.merged.sam > ${SAMPLE}.merged_new.bam
        #samtools view -H ${SAMPLE}.merged.bam | sed '/gl00/d' | samtools reheader - ${SAMPLE}.merged.bam > ${SAMPLE}.merged_new.bam
        rm ${SAMPLE}.merged.bam ${SAMPLE}.merged.sam
		
		##sort by coordiante
		java -jar ${picard_path}/SortSam.jar INPUT=${SAMPLE}.merged_new.bam OUTPUT=${SAMPLE}.sort.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
		##order according to the refernen genome
		java -jar ${picard_path}/ReorderSam.jar INPUT=${SAMPLE}.sort.bam OUTPUT=${SAMPLE}.reorder.bam REFERENCE=${ref_fa} VALIDATION_STRINGENCY=SILENT
		##add read group...
		java -jar ${picard_path}/AddOrReplaceReadGroups.jar INPUT=${SAMPLE}.reorder.bam OUTPUT=${SAMPLE}.final.bam SORT_ORDER=coordinate RGLB=8 RGPL=Illumina RGPU=1 RGSM=1804
	    ##add index file
		samtools index ${SAMPLE}.final.bam    #.-----> final output file
		##clean up
		rm ${SAMPLE}.merged_new.bam ${SAMPLE}.sort.bam ${SAMPLE}.reorder.bam
done

date#echo the time at end



