#STAR v2.7.10b
gDir="STAR"
gtf="Homo_sapiens.GRCh38.108.chr.gtf"
fastq="shRNA"
outDir="STAR_Mapped"
file=$1

/work/users/g/g/ggiri/Tools/STAR/bin/Linux_x86_64_static/STAR \
        --runThreadN 16 \
        --sjdbGTFfile ${gtf} \
        --genomeDir ${gDir} \
        --readFilesIn ${fastq}/${file}_R1_001.fastq.gz ${fastq}/${file}_R2_001.fastq.gz \
        --readFilesCommand zcat \
        --genomeLoad NoSharedMemory \
        --outFileNamePrefix ${outDir}/${file} \
        --outSAMtype BAM Unsorted \
        --quantMode TranscriptomeSAM

samtools sort -@ 16 ${outDir}/${file}Aligned.out.bam > ${outDir}/${file}Aligned.sorted.bam
samtools index -@ 16 ${outDir}/${file}Aligned.sorted.bam
