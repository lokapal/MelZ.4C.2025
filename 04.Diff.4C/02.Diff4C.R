#!/usr/bin/Rscript
# script to perform differential 4C analysis 
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Input:  counts.txt pre-computed file by featureCounts
# Output: 1. MelZ.4C.results.tsv          MelZ plastic-matrigel 4C Differential analysis tabbed file
#         2. MelZ.4C.scatterplots.pdf     Scatterplots to control replicates consistency
#         3. MelZ.4C.volcanoplot.pdf      VolcanoPlot to display the most prominent 4C down/upcontacted genes
#         4. MelZ.4C.PCA.pdf              PCA to control replicates consistency
#
# Dependency tools & libraries:
# 1. R
# Bioconductor packages:
# 2. DESeq2 R             https://bioconductor.org/packages/release/bioc/html/DESeq2.html
# 3. EnhancedVolcano      https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html
# 4. PCAExplorer          https://bioconductor.org/packages/release/bioc/html/pcaExplorer.html
# Common R libraries:
# 5. dplyr, ggplot2, tibble, RColorBrewer, gplots, ggrepel

# Import data from featureCounts
countdata <- read.table("counts.txt", header=TRUE, row.names=1, stringsAsFactors = F)

# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]

# Remove .bam or .sam from filenames
colnames(countdata) <- gsub("\\.inter\\.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub("intersect.nodfam", "", colnames(countdata))
colnames(countdata) <- gsub("...", "", colnames(countdata), fixed=T)
colnames(countdata) <- gsub("..", ".", colnames(countdata), fixed=T)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(powerjoin))

# Convert to matrix
countdata <- round(as.matrix(countdata))

# Assign condition (first two are expansion, second two are control)
(condition <- factor(c(rep("exp", 2), rep("ctl", 2))))

# Analysis with DESeq2 ----------------------------------------------------

suppressPackageStartupMessages(library(DESeq2))

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- dds[ rowSums(counts(dds)) > 0, ]


vsd <- varianceStabilizingTransformation(dds, blind=FALSE, fitType="local")

rld <- rlogTransformation(dds, fitType="local")

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("tibble"))

ddc <- estimateSizeFactors(dds)
head(assay(ddc))

df <- bind_rows(
  as_tibble(log2(counts(ddc, normalized=TRUE)[, 1:2]+1)) %>% mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "VarianceStabilizingTransformation"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rLogTransformation (DESeq2 default)"))

colnames(df)[1:2] <- c("x", "y")

cairo_pdf("qc-scatterplots.pdf",width=15,height=10,antialias="default",fallback_resolution = 300,onefile=T)
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation) + labs(x = "rep1.Matrigel", y="rep2.Matrigel")

rm (df)
df <- bind_rows(
  as_tibble(log2(counts(ddc, normalized=TRUE)[, 3:4]+1)) %>% mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, 3:4]) %>% mutate(transformation = "VarianceStabilizingTransformation"),
  as_tibble(assay(rld)[, 3:4]) %>% mutate(transformation = "rLogTransformation (DESeq2 default)"))

colnames(df)[1:2] <- c("x", "y")
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation) + labs(x="rep1.Plastic", y="rep2.Plastic")
invisible(dev.off())

# Run the DESeq pipeline
dds <- DESeq(dds, fitType="local")

# Plot dispersions
#png("qc-dispersions.png", 1000, 1000, pointsize=20)
cairo_pdf("qc-dispersions.pdf",width=15,height=10,antialias="default",fallback_resolution = 300)
plotDispEsts(dds, main="Dispersion plot")
invisible(dev.off())

suppressPackageStartupMessages(library(RColorBrewer))
mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))]

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
suppressPackageStartupMessages(library(gplots))
#png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
cairo_pdf("qc-heatmap.pdf",width=15,height=10,antialias="default",fallback_resolution = 300)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
invisible(dev.off())

suppressPackageStartupMessages(require(pcaExplorer))
cairo_pdf("qc-pca.pdf",width=15,height=10,antialias="default",fallback_resolution = 300,onefile=T)
pcaplot(rld,intgroup="condition",title = "RlogStabilizingTrans 4C (default")
pcaplot(vsd,intgroup="condition",title = "VarianceStabilizingTrans 4C")
invisible(dev.off())

# Get differential expression results
res <- results(dds)

resD <- as.data.frame(res)

finalnames <- as.data.frame(read.table("genes.allnames", sep = "\t", header=T, row.names=NULL, quote="", comment.char = "", stringsAsFactors = F))

resD$rnames <- rownames(res)
resD <- resD %>% dplyr::left_join(finalnames, by = c("rnames" = "gene_id")) %>%
                 dplyr::rename(GeneID = rnames) %>%
                 select (GeneID, gene_name, everything())

resD <- data.frame(resD, row.names = 'GeneID')

## Order by adjusted p-value
resD <- resD[order(resD$padj), ]

## Merge with normalized count data
resdata <- merge(resD, as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
#resdata <- merge(resD, as.data.frame(countdata), by="row.names", sort=FALSE)
names(resdata)[1] <- "GeneID"
resdata <- resdata %>% select(GeneID, gene_name, everything())

## Write results
write.table(subset(resdata, !is.na(pvalue)), file="MelZ.4C.results.tsv",sep='\t', row.names=F, quote=F)

upreg   <- subset(resdata, log2FoldChange >=  1.5 & pvalue <= 0.055 )
downreg <- subset(resdata, log2FoldChange <= -1.5 & pvalue <= 0.055 )
write.table(upreg,   file="MelZ.4C.upreg.1.5.tsv", sep='\t', quote=FALSE, row.names=FALSE)
write.table(downreg, file="MelZ.4C.downreg.1.5.tsv", sep='\t', quote = FALSE, row.names=FALSE)

upreg   <- subset(resdata, log2FoldChange >=  2.0 & pvalue <= 0.055 )
downreg <- subset(resdata, log2FoldChange <= -2.0 & pvalue <= 0.055 )
write.table(upreg,   file="MelZ.4C.upreg.2.0.tsv", sep='\t', quote=FALSE, row.names=FALSE)
write.table(downreg, file="MelZ.4C.downreg.2.0.tsv", sep='\t', quote = FALSE, row.names=FALSE)

#png("MelZ.4C.maplot.png", 1500, 1000, pointsize=20)
cairo_pdf("MelZ.4C.maplot.pdf",width=15,height=10,pointsize=16,antialias="default",fallback_resolution = 300,onefile=T)
DESeq2::plotMA(res, ylim=c(-1,1))
invisible(dev.off())

suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("EnhancedVolcano"))

cairo_pdf("MelZ.4C.volcanoplot.pdf",width=15,height=10,pointsize=16,antialias="default",fallback_resolution = 300,onefile=T)
EnhancedVolcano(resD, lab=resD$gene_name, x='log2FoldChange', y='pvalue',
    xlab = bquote(~Log[2]~ 'fold change'),
    pCutoff = 5e-2,
    FCcutoff = 2.0,
    pointSize = 1.0,
    labSize = 4.0,
    colAlpha = 4/5,
    legendLabels=c('NS', 'Log2 FC', 'P value', 'P value & Log2 FC'),
    legendPosition = 'bottom',
    legendLabSize = 10,
    legendIconSize = 3.0,
#    drawConnectors = TRUE,
#    widthConnectors = 0.2,
#    colConnectors = 'grey30',
) + coord_cartesian(xlim=c(-14,14),ylim=c(0,310)) + scale_x_continuous(breaks=seq(-30,30,5), minor_breaks=seq(-29,29,1)) #+ scale_y_continuous(breaks=seq(0,40,10), minor_breaks=seq(1,40,5))


EnhancedVolcano(resD,lab=resD$gene_name,x='log2FoldChange',y='padj',
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~-Log[10]~adjusted~italic(P)),
    pCutoff = 5e-2,
    FCcutoff = 2.0,
    pointSize = 1.0,
    labSize = 4.0,
    colAlpha = 4/5,
    legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
    legendPosition = 'bottom',
    legendLabSize = 10,
    legendIconSize = 3.0,
#    drawConnectors = TRUE,
#    widthConnectors = 0.2,
#    colConnectors = 'grey30'
) + coord_cartesian(xlim=c(-14,14),ylim=c(0,300)) + scale_x_continuous(breaks=seq(-30,30,5), minor_breaks=seq(-29,29,1)) #+ scale_y_continuous(breaks=seq(0,40,10), minor_breaks=seq(1,40,5))

invisible(dev.off())
