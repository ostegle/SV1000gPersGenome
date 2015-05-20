annotation=/homes/stegle/research/users/stegle/statgenom-100grnaseq/data/data/personalized_genome/reference/NA12878_diploid_1kgp3/maternal/gencode.v19.annotation.maternal.bed.fixed.gtf
alignment=/homes/stegle/research/users/stegle/statgenom-100grnaseq/data/data/personalized_genome/rna-seq/alignments/GRCH37_NA12878_maternal/Aligned.sortedByCoord.out.bam
out_dir=/homes/stegle/research/users/stegle/statgenom-100grnaseq/data/data/personalized_genome/rna-seq/counts/GRCH37_NA12878_maternal

python ./../pysrc/counting.py $alignment $annotation $out_dir
