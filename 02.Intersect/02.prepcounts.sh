#!/bin/sh
# Dependency tools:
# 1. featureCounts  https://subread.sourceforge.net/
# set cell treatment/growth type before use. E.g. plastic, matrigel, untreated, hemin etc.
TYPE=plastic
featureCounts -a /usr/local/genomes/hg38.112.gtf -o counts.txt -T 30 -t gene -g gene_id -O -p --countReadPairs --fraction --readExtension5 2500 --readExtension3 2500 -p $TYPE.rep1.inter.bam $TYPE.rep2.inter.bam
