#!/bin/bash
# TYPE should be set to plastic and/or matrigel
TYPE=platic
# REP should be set to rep1 and/or rep2
REP=rep1
repair.sh -Xmx50g in1=$TYPE.$REP.R1.fastq.gz in2=$TYPE.$REP.R2.fastq.gz \
out1=4C.$TYPE.$REP.filtered.R1.fastq.gz out2=4C.$TYPE.$REP.filtered.R2.fastq.gz \
outs=4C.MelZ.filtered.singletons.fastq.gz repair overwrite=t
