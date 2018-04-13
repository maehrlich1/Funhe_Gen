for i in *R1*; do

    [[ $i =~ (.*)_fc(.*)_R1.fastq.gz ]]

    ~/software/bwa-0.7.15/bwa mem -M -t 16 \
        -R "@RG\tID:${BASH_REMATCH[1]}.fc${BASH_REMATCH[2]}\tSM:${BASH_REMATCH[1]}\tPL:illumina\tLB:${BASH_REMATCH[1]}.lib1\tPU:fc${BASH_REMATCH[2]}" \
        ../reference/Mnemiopsis_leidyi.GCA_000226015.1.31.dna.genome.red10kb.fa \
        $i ${BASH_REMATCH[1]}_fc${BASH_REMATCH[2]}_R2.fastq.gz \
        > ../mapped/${BASH_REMATCH[1]}_fc${BASH_REMATCH[2]}.sam

done