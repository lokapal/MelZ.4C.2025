#!/bin/sh
featureCounts -a /usr/local/genomes/hg38.112.gtf -o counts.txt -T 30 -t gene -g gene_id -O -p --countReadPairs --fraction --readExtension5 2500 --readExtension3 2500 -p \
../Matrigel/intersect.nodfam/rep1.inter.bam ../Matrigel/intersect.nodfam/rep2.inter.bam \
../Plastic/intersect.nodfam/rep1.inter.bam ../Plastic/intersect.nodfam/rep2.inter.bam
