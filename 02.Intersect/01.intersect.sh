#!/bin/sh
# Dependency tools:
# 1. bedtools
# 2. samtools
# 3. R
# set cell treatment/growth type before use. E.g. plastic, matrigel, untreated, hemin etc.
TYPE=plastic
# convert noDFAM BAMs to BED files and compress to gzip
bedtools bamtobed -i ../$TYPE/$TYPE.rep1.nodfam.bam > $TYPE.rep1.nodfam.bed
pigz -f $TYPE.rep1.nodfam.bed
bedtools bamtobed -i ../$TYPE/$TYPE.rep2.nodfam.bam > $TYPE.rep2.nodfam.bed
pigz -f $TYPE.rep2.nodfam.bed
# find overlaps between $TYPE.rep1.nodfam.bed.gz and $TYPE.rep2.nodfam.bed.gz
Rscript overlap_reps_names.R $TYPE. nodfam
# find uniq entries in overlap file rep1
pigz -dc $TYPE.nodfam.rep1.overlap.txt.gz | sort -k1,1 | uniq > $TYPE.nodfam.rep1.overlap.txt
# remove non-intersected entries from rep1 noDFAM alignment file 
samtools view -@30 -N $TYPE.nodfam.rep1.overlap.txt $TYPE.rep1.nodfam.bam -o $TYPE.rep1.inter.bam
samtools index -@30 $TYPE.rep1.inter.bam
# find uniq entries in overlap file rep2
pigz -dc $TYPE.nodfam.rep2.overlap.txt.gz | sort -k1,1 | uniq > $TYPE.nodfam.rep2.overlap.txt
# remove non-intersected entries from rep2 noDFAM alignment file 
samtools view -@30 -N $TYPE.nodfam.rep2.overlap.txt $TYPE.rep2.nodfam.bam -o $TYPE.rep2.inter.bam
samtools index -@30 $TYPE.rep2.inter.bam
