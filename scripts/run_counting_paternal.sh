annotation=/homes/stegle/research/users/stegle/statgenom-100grnaseq/data/data/personalized_genome/reference/NA12878_diploid_1kgp3/paternal/gencode.v19.annotation.paternal.bed.fixed.gtf
alignment=/homes/stegle/research/users/stegle/statgenom-100grnaseq/data/data/personalized_genome/rna-seq/alignments/GRCH37_NA12878_paternal/Aligned.sortedByCoord.out.bam
out_dir=/homes/stegle/research/users/stegle/statgenom-100grnaseq/data/data/personalized_genome/rna-seq/counts/GRCH37_NA12878_paternal

python ./../pysrc/counting.py $alignment $annotation $out_dir
