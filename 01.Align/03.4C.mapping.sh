#!/bin/bash
# script to align to genome previously filtered from adapters, primers and low quality single end reads
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Input:  1. GSEXXXXX: GSMXXXX, GSMXXXX 	MelZ untreated 4C-rDNA reads after adapter removing
# Output: 1. table_hg38.plastic.rep1.txt the file with the genome-wide hg38 DSB mappings: chromosome coordinates, alignment length, coverage, reads, sequence
#         2. plastic.rep1.bam       sorted alignment file that doesn't contain unaligned reads
#            plastic.rep1.bam.bai   index file for plastic.rep1.bam
#
# Dependency tools:
# 1. bwa                  http://bio-bwa.sourceforge.net/
# 2. samtools             http://www.htslib.org
# 3. tabix from samtools 
# 4. Perl with BIO::DB::Fasta BioPerl library
#
# set cell treatment/growth type before use. E.g. plastic, matrigel, untreated, hemin etc.
TYPE=plastic
# set replicate number before use. rep1 for the first replicate, rep2 for the second etc.
REP=rep1
# samtools flags:
# -F 2052 filters out both unaligned reads and supplementary alignments
# -F 4    filters out unaligned reads only
# -F 2048 filters out supplementary alignments only
# 4C anchor and similar regions should be removed from the alignment, i.e.:
# chrMT, chrY (MelZ line is a female cell line), chrGL000220 (that contains rDNA genes), and 
# some regions on chr21 that contain rDNA genes too
bwa mem -t 20 /usr/local/genomes/hg38.mfa $TYPE.$REP.filtered.R1.fastq.gz $TYPE.$REP.filtered.R2.fastq.gz |\
samtools view --threads 8 -bS -F 2052 | samtools view --threads 8 - -h -o inchr21.bam -U - -L rDNA21.bed | sed '/chrMT/d;/chrGL000220/d;/chrY/d' |\
samtools sort -@ 10 -O BAM -o hits.bam
samtools mpileup -f /usr/local/genomes/hg38.mfa hits.bam -o pileup.txt
samtools view -@ 10 -O SAM hits.bam -o hits.sam
mv hits.bam $TYPE.$REP.bam
samtools index -@ 10 $TYPE.$REP.bam
perl ./lib/makeTablePileup.pl pileup.txt > table.txt
perl ./lib/addNumOfReads.pl hits.sam table.txt > tableA.txt
perl ./lib/addSubseq.pl /usr/local/genomes/hg38.mfa tableA.txt > table_hg38.$TYPE.$REP.txt
rm -f *.sam pileup.* inchr21.bam
