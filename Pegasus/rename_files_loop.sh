for i in *; do
	[[ $i =~ (.*)16_(.*)_S(.*)_(.*)_001.fastq.gz ]]
	mv $i ${BASH_REMATCH[1]}_${BASH_REMATCH[2]}_${BASH_REMATCH[4]}.fastq.gz
done