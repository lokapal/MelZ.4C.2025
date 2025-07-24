#!/bin/sh
# dependence deepTools2 https://deeptools.readthedocs.io
# set cell treatment/growth type before use. E.g. plastic, matrigel, untreated, hemin etc.
TYPE=plastic
bamCoverage --bam $TYPE.rep1.inter.bam -o rep1.inter.bg -of bedgraph -bs 10 -p 32 --effectiveGenomeSize 2913022398 \
            --normalizeUsing RPKM --exactScaling --ignoreForNormalization chrY chrMT chrGL000220 chr21 --skipNAs
bamCoverage --bam $TYPE.rep2.inter.bam -o rep2.inter.bg -of bedgraph -bs 10 -p 32 --effectiveGenomeSize 2913022398 \
            --normalizeUsing RPKM --exactScaling --ignoreForNormalization chrY chrMT chrGL000220 chr21 --skipNAs
./lib/joinbedgraph.sh rep1.inter.bg rep2.inter.bg MelZ.4C.$TYPE.bedGraph
