#BSUB -n 6
#BSUB -q research-rh6 
#BSUB -R "rusage[mem=65000,tmp=1000]"
#BSUB -M 65000
#out dir
#BSUB -o /homes/stegle/research/projects/1000GenomesRNASeq/data/personalized_genome/scripts/cluster_out

annotation=personal_transcriptome/gencode.v19.annotation.paternal.bed.fixed.gtf
genome_dir=STAR_paternal
genome=NA12878_paternal.fa
tmp_dir=_STAR_paternal_tmp

STAR --runMode genomeGenerate --outTmpDir $tmp_dir --genomeDir $genome_dir --genomeFastaFiles $genome --runThreadN 6 --sjdbGTFfile $annotation --sjdbOverhang 100
