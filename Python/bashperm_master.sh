#! /bin/bash
### This script performs permutations on random data
### It requires script rand_ind_master.py to be in the same directory
### To silence console, run "bashperm_all.sh &> /dev/null"

for file in *.vcf ## VCF with individuals and genotypes
do
	for run in {1..3} ## Change {1..100} to any number of permutations i.e. {1..500}
	do
		OUTPUT=$(python 1rand_ind_master.py) ### Run your python script to randomize samples
		
		### Now run the test on randomized groups
		/Users/admin/Programs/vcftools/bin/vcftools --vcf $file --weir-fst-pop Atest.txt --weir-fst-pop Btest.txt --weir-fst-pop Ctest.txt --weir-fst-pop Dtest.txt --out mktemp
		awk 'BEGIN{FS="\t"}{getline f1 <"FstVal.txt" ;print f1,$3}' OFS="\t" mktemp.weir.fst > tmptxt.txt
		mv tmptxt.txt FstVal.txt
	done
	> fst_$file.txt
	mv FstVal.txt fst_$file.txt
done
