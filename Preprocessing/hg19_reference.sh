
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*'
gunzip *.gz
for C in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M
do
cat /hms/scratch1/wang/chr${C}.fa >> hg19.ucsc.fa
done
#creat .dict
java -jar picard-tools-1.107/CreateSequenceDictionary.jar R= hg19.ucsc.fa O= hg19.ucsc.dict
#creat .fai
samtools faidx hg19.ucsc.fa