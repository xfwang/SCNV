#!/bin/bash

#.bam file --> coverage file
##
##
work_dir="/hms/scratch1/wang/work"  #working directory
in_dir="/hms/scratch1/wang/batch"  #input file directory
bound_file="hg19.bin.boundaries.50k.bowtie.k50.sorted.txt" #boundary file name

date

cd ${work_dir}
##
	if [ -f rc_matrix.txt ]; then
	   rm rc_matrix.txt
	fi
	touch rc_matrix.txt #empty file
	
for i in `cat bam_list.txt`  #.bam file list
	do
	SAMPLE=$i  #cell file name
	samtools view ${in_dir}/${SAMPLE}.final.bam | awk '{print $3,$4}' > ${SAMPLE}.rc
	Rscript read_count.Rscript ${SAMPLE}.rc ${bound_file}

	#merge	
	if [ -s rc_matrix.txt ]; then
		gawk '{print $5}' ${SAMPLE}.rc.out > ${SAMPLE}.rc5.out
		paste -d ' ' rc_matrix.txt ${SAMPLE}.rc5.out > rc_matrix.txt
		rm ${SAMPLE}.rc5.out
	else 
        paste -d ' ' rc_matrix.txt ${SAMPLE}.rc.out > rc_matrix.txt 
	fi 
	
	#clean up
	rm ${SAMPLE}.rc
	rm ${SAMPLE}.rc.out 
	echo ${SAMPLE} "done"
	
done

