#BSUB -n 6
#BSUB -q research-rh6 
#BSUB -R "rusage[mem=65000,tmp=1000] select[gpfs]"
#BSUB -M 65000
#out dir
#BSUB -o /homes/stegle/research/projects/1000GenomesRNASeq/data/personalized_genome/scripts/cluster_out

annotation=personal_transcriptome/gencode.v19.annotation_trunc.gtf
genome_dir=STAR
genome=b37_g1k_phase2.fa

STAR --runMode genomeGenerate --genomeDir $genome_dir --genomeFastaFiles $genome --runThreadN 6 --sjdbGTFfile $annotation --sjdbOverhang 100
