#RSEM v1.3.1
index="Human/RSEM/hg38"
bamdir="../STAR_mapped"
outdir="../RSEM"
file=$1

rsem-calculate-expression --paired-end \
        --num-threads 32 \
        --alignments \
        --seed 77 \
        --no-bam-output \
        ${bamdir}/${file}Aligned.toTranscriptome.out.bam \
        ${index} \
        ${outdir}/${file}
