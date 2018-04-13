for i in *.bam; do
	[[ $i =~ (.*)_merge_mt_sorted.bam ]]
	mkdir ${BASH_REMATCH[1]}_mt_ass
	mkdir ${BASH_REMATCH[1]}_mt_ass/Scaffold
	samtools sort -n $i | velveth ${BASH_REMATCH[1]}_mt_ass 31 -bam -shortPaired -
	velvetg ${BASH_REMATCH[1]}_mt_ass -cov_cutoff auto -exp_cov auto -ins_length 500 -min_contig_lgth 31 -alignments yes
	cd ${BASH_REMATCH[1]}_mt_ass/Scaffold
	abacas -r ../../../../reference_genome/Mnemiopsis_leidyi.GCA_000226015.1.31.dna.mt.fa -q ../contigs.fa -p nucmer
	cd ../..;
done