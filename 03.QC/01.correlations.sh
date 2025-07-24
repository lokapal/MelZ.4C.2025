#!/bin/bash
#hg19 2864785220
#hg38 2913022398
# dependence deepTools2 https://deeptools.readthedocs.io
# set cell treatment/growth type before use. E.g. plastic, matrigel, untreated, hemin etc.
TYPE=plastic
bamCoverage --effectiveGenomeSize 2913022398 --normalizeUsing RPKM --ignoreForNormalization chrY chrMT chrGL000220 chr21 \
--skipNAs --exactScaling -p 30 -b ../01.Align/$TYPE.rep1.bam -o rep1.bw
bamCoverage --effectiveGenomeSize 2913022398 --normalizeUsing RPKM --ignoreForNormalization chrY chrMT chrGL000220 chr21 \
--skipNAs --exactScaling -p 30 -b ../01.Align/$TYPE.rep2.bam -o rep2.bw
multiBigwigSummary bins -p 12 -b rep1.bw rep2.bw -o res12hs.npz
plotCorrelation -in res12hs.npz -c spearman --removeOutliers --skipZeros -p scatterplot -o 4C.MelZ.$TYPE.spearman.pdf --labels 4C.$TYPE.rep1 4C.$TYPE.rep2 --log1p
plotCorrelation -in res12hs.npz -c pearson  --removeOutliers --skipZeros -p scatterplot -o 4C.MelZ.$TYPE.pearson.pdf  --labels 4C.$TYPE.rep1 4C.$TYPE.rep2 --log1p
plotPCA -in res12hs.npz -o MelZ.4C.PCA.pdf -T "PCA of MelZ 4C read normalized counts"
