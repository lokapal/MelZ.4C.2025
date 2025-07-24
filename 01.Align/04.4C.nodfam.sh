#!/bin/sh
# remove chrMT and chrY from DFAM database to pass bedtools limitations
pigz -dc hg38_dfam_full.bed.gz | sed '/chrMT/d;/chrY/d' | pigz > hg38_dfam.bed.gz
# set cell treatment/growth type before use. E.g. plastic, matrigel, untreated, hemin etc.
TYPE=plastic
# set replicate number before use. rep1 for the first replicate, rep2 for the second etc.
REP=rep1
bedtools intersect -sorted -wa -v -f 1.0 -a $TYPE.$REP.bam -b hg38_dfam.bed.gz > $TYPE.$REP.nodfam.bam
samtools index -@30 $TYPE.$REP.nodfam.bam
