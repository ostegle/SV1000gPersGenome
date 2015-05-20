#BSUB -n 8
#BSUB -q research-rh6 
#BSUB -R "rusage[mem=65000,tmp=1000] select[gpfs]"
#BSUB -M 65000
#out dir
#BSUB -o /homes/stegle/research/projects/1000GenomesRNASeq/data/personalized_genome/scripts/cluster_out

annotation=personal_transcriptome/gencode.v19.annotation_trunc.gtf
genome_dir=STAR
genome=b37_g1k_phase2.fa

fasta_1=/homes/stegle/research/users/stegle/statgenom-100grnaseq/data/data/personalized_genome/rna-seq/raw/ERR356372_1.fastq
fasta_2=/homes/stegle/research/users/stegle/statgenom-100grnaseq/data/data/personalized_genome/rna-seq/raw/ERR356372_2.fastq

alignment_base=alignment/

nthreads=8
mem_thread=65000
multi_map_max=10

#run STAR
STAR --genomeDir $genome_dir --outFilterMultimapNmax $multi_map_max --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --outFilterScoreMinOverLread 0.33 --outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --runThreadN $nthreads --readFilesIn $fasta_1 $fasta_2 --outFileNamePrefix $alignment_base

#index
samtools index $alignment_base/Aligned.sortedByCoord.out.bam
