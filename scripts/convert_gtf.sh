awk -F'\t' -v OFS="\t" '{ print $1, $4, $5, $2, $3, $6, $7, $8, $9}' gencode.v19.annotation.maternal.bed.gtf > gencode.v19.annotation.maternal.bed.fixed.gtf

awk -F'\t' -v OFS="\t" '{ print $1, $4, $5, $2, $3, $6, $7, $8, $9}' gencode.v19.annotation.paternal.bed.gtf > gencode.v19.annotation.paternal.bed.fixed.gtf
