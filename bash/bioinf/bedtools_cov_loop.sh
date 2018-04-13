for i in *sorted.bam; do
	[[ $i =~ (.*)_sorted.bam ]]
	bedtools genomecov -ibam $i -g ../../reference_genome/Mnemiopsis_leidyi.GCA_000226015.1.31.dna.genome.fa.gz | tee ../coverage/${BASH_REMATCH[1]}_cov_contigs.txt | grep 'genome' > ../coverage/${BASH_REMATCH[1]}_cov_genome.txt;
done