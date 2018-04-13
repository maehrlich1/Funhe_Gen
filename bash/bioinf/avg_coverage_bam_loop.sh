> ../coverage/avg_cov.txt
for i in *sorted.bam; do
	[[ $i =~ (.*)_sorted.bam ]]
	paste <(echo ${BASH_REMATCH[1]}) <(samtools depth $i | awk '{sum+=$3}END{print sum/158476310}') >> ../coverage/avg_cov.txt;
done