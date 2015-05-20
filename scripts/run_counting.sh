annotation=/homes/stegle/research/users/stegle/statgenom-100grnaseq/data/data/personalized_genome/reference/GRCH37/gencode.v19.annotation.gtf
alignment=/homes/stegle/research/users/stegle/statgenom-100grnaseq/data/data/personalized_genome/rna-seq/alignments/GRCH37/Aligned.sortedByCoord.out.bam
out_dir=/homes/stegle/research/users/stegle/statgenom-100grnaseq/data/data/personalized_genome/rna-seq/counts/GRCH37

python ./../pysrc/counting.py $alignment $annotation $out_dir
