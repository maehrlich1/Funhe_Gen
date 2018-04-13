#!/bin/bash
echo "The script runs Trimmomatic in paired end mode for all reads in the working directory of the form *R1*.gz (forward) matching *_R2*.gz (reverse read)."

for i in *R1*; do
	[[ $i =~ (.*)R1.fastq.gz ]]
	java -jar /Applications/Bioinf/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 $i ${BASH_REMATCH[1]}R2.fastq.gz /Users/Moritz/Documents/Academic/GEOMAR/Thesis/Seq/miseq/trim_reads_best/${BASH_REMATCH[1]}fw_p.fq.gz /Users/Moritz/Documents/Academic/GEOMAR/Thesis/Seq/miseq/trim_reads_best/${BASH_REMATCH[1]}fw_u.fq.gz /Users/Moritz/Documents/Academic/GEOMAR/Thesis/Seq/miseq/trim_reads_best/${BASH_REMATCH[1]}rv_p.fq.gz /Users/Moritz/Documents/Academic/GEOMAR/Thesis/Seq/miseq/trim_reads_best/${BASH_REMATCH[1]}rv_u.fq.gz ILLUMINACLIP:/Applications/Bioinf/Trimmomatic-0.36/adapters/TruSeq3-PE-2_NexteraPE.fa:1:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36;
done

# use of Trimmomatic:
# java -jar <path to trimmomatic.jar> PE [-threads <threads] [-phred33 | -phred64] [-trimlog <logFile>] >] [-basein <inputBase> | <input 1> <input 2>] [-baseout <outputBase> | <unpaired output 1> <paired output 2> <unpaired output 2> ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>	