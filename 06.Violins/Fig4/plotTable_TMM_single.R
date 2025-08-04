#!/usr/local/bin/Rscript
# script to draw violin plots from table made by 
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Inputs:  1. counts.byTMMs.full.names with differential 4C, made by edgeR
#          2. 4C.Venn.docx with gene lists made by Venn server
# Outputs: 1. Violin_TMM_single.pdf with violins
#          2. MWtest.txt table with Mann-Whitney pair-wise statistics table
#
# Dependency tools & libraries:
# 1. R
# Common R libraries:
# 2. dplyr
# 3. docxtractr
# 4. ggplot2
# 5. scales
# 6. ggsci

give.n <- function(x){
  return(c(y = log10(50000), label = length(x))) 
}

library(docxtractr)
library(dplyr)

doc <- read_docx("4C.Venn.docx")
doctbl <- docx_extract_tbl(doc, header=TRUE)
head(doctbl)

Plastic   <- doctbl %>% filter(Names == "Plastic") %>% pull(elements)
Decreased <- doctbl %>% filter(Names == "Decreased") %>% pull(elements)
Increased <- doctbl %>% filter(Names == "Increased") %>% pull(elements)
PlastDecr <- doctbl %>% filter(Names == "Plastic.decreased") %>% pull(elements)
PlastIncr <- doctbl %>% filter(Names == "Plastic.increased") %>% pull(elements)

Plastic   <- strsplit(Plastic, " ")[[1]]
Decreased <- strsplit(Decreased, " ")[[1]]
Increased <- strsplit(Increased, " ")[[1]]
PlastIncr <- strsplit(PlastIncr, " ")[[1]]
PlastDecr <- strsplit(PlastDecr, " ")[[1]]


F4Call <- as.data.frame(read.table("counts.byTMMs.full.names", sep = "\t", header = T, stringsAsFactors = FALSE, row.names = NULL))
F4Call$MelZ.M.TPM <- replace(F4Call$MelZ.M.TPM, F4Call$MelZ.M.TPM<0.001, 0.001) # replace expression values below 0.001 to 0.001 to use log10 scale
F4Call$MelZ.P.TPM <- replace(F4Call$MelZ.P.TPM, F4Call$MelZ.P.TPM<0.001, 0.001) # replace expression values below 0.001 to 0.001 to use log10 scale
F4Call$MelZ.Plastic.Mean.TMM <- replace(F4Call$MelZ.Plastic.Mean.TMM, F4Call$MelZ.Plastic.Mean.TMM<0.1, 0.1) # replace TMM values below 0.1 to 0.1 to use log10 scale
F4Call$MelZ.Matrigel.Mean.TMM <- replace(F4Call$MelZ.Matrigel.Mean.TMM, F4Call$MelZ.Matrigel.Mean.TMM<0.1, 0.1) # replace TMM values below 0.1 to 0.1 to use log10 scale

Plastic   <- as.data.frame(Plastic)
colnames(Plastic)[1] <- "Genes"
PlastFinal <- dplyr::left_join(Plastic, F4Call, by = c("Genes" = "GeneName"))
PlastPlot <- PlastFinal %>% dplyr::select (MelZ.Plastic.Mean.TMM, MelZ.Matrigel.Mean.TMM) %>% dplyr::rename (P.TMM = MelZ.Plastic.Mean.TMM, M.TMM = MelZ.Matrigel.Mean.TMM)
TP <- data.frame("data"=PlastPlot$P.TMM,key="Plastic.P")
TM <- data.frame("data"=PlastPlot$M.TMM,key="Plastic.M")
Tests <- rbind (TP,TM)

Decreased  <- as.data.frame(Decreased)
colnames(Decreased)[1] <- "Genes"
DecrFinal <- dplyr::left_join(Decreased, F4Call, by = c("Genes" = "GeneName"))
DecrPlot  <- DecrFinal %>% dplyr::select (MelZ.Plastic.Mean.TMM, MelZ.Matrigel.Mean.TMM) %>% dplyr::rename (P.TMM = MelZ.Plastic.Mean.TMM, M.TMM = MelZ.Matrigel.Mean.TMM)
TP <- data.frame("data"=DecrPlot$P.TMM,key="Decreased.P")
TM <- data.frame("data"=DecrPlot$M.TMM,key="Decreased.M")
Tests <- rbind (Tests,TP,TM)

Increased  <- as.data.frame(Increased)
colnames(Increased)[1] <- "Genes"
IncrFinal <- dplyr::left_join(Increased, F4Call, by = c("Genes" = "GeneName"))
IncrPlot  <- IncrFinal %>% dplyr::select (MelZ.Plastic.Mean.TMM, MelZ.Matrigel.Mean.TMM) %>% dplyr::rename (P.TMM = MelZ.Plastic.Mean.TMM, M.TMM = MelZ.Matrigel.Mean.TMM)
TP <- data.frame("data"=IncrPlot$P.TMM,key="Increased.P")
TM <- data.frame("data"=IncrPlot$M.TMM,key="Increased.M")
Tests <- rbind (Tests,TP,TM)

PlastIncr <- as.data.frame(PlastIncr)
colnames(PlastIncr)[1] <- "Genes"
PlastIncrFinal <- dplyr::left_join(PlastIncr, F4Call, by = c("Genes" = "GeneName"))
PlastIncrPlot  <- PlastIncrFinal %>% dplyr::select (MelZ.Plastic.Mean.TMM, MelZ.Matrigel.Mean.TMM) %>% dplyr::rename (P.TMM = MelZ.Plastic.Mean.TMM, M.TMM = MelZ.Matrigel.Mean.TMM)
TP <- data.frame("data"=PlastIncrPlot$P.TMM,key="PlastIncr.P")
TM <- data.frame("data"=PlastIncrPlot$M.TMM,key="PlastIncr.M")
Tests <- rbind (Tests,TP,TM)

PlastDecr <- as.data.frame(PlastDecr)
colnames(PlastDecr)[1] <- "Genes"
PlastDecrFinal <- dplyr::left_join(PlastDecr, F4Call, by = c("Genes" = "GeneName"))
PlastDecrPlot  <- PlastDecrFinal %>% dplyr::select (MelZ.Plastic.Mean.TMM, MelZ.Matrigel.Mean.TMM) %>% dplyr::rename (P.TMM = MelZ.Plastic.Mean.TMM, M.TMM = MelZ.Matrigel.Mean.TMM)
TP <- data.frame("data"=PlastDecrPlot$P.TMM,key="PlastDecr.P")
TM <- data.frame("data"=PlastDecrPlot$M.TMM,key="PlastDecr.M")
Tests <- rbind (Tests,TP,TM)

#head(Tests)
#FinalPlot  <- tidyr::pivot_longer(Tests, cols = 1:2)
#head(FinalPlot)
#stop()

#PlastDecrPlot  <- tidyr::pivot_longer(PlastDecrPlot, cols = 1:2)
#PlastDecrPlot  <- PlastDecrPlot %>% mutate (Condition='PlasticDecreased')
#PlastDecrPlot$name <- factor(PlastDecrPlot$name, levels = c("P.TMM", "M.TMM"))

FinalPlot <- Tests %>% arrange(data) %>%
  mutate(key = factor(key, levels=c("Plastic.P", "Plastic.M", "Decreased.P", "Decreased.M", "Increased.P", "Increased.M", "PlastDecr.P", "PlastDecr.M", "PlastIncr.P", "PlastIncr.M")))


suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(hrbrthemes))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggsci))

cairo_pdf("violin_TMM_single.pdf",width=10,height=5,antialias="default",fallback_resolution = 600,onefile=T)
ggplot(FinalPlot,aes(x=key, y=data, fill=key)) +
  theme_bw() + 
  geom_violin(position="identity", alpha=0.6, lwd=0.4, width=2.0,trim=T) + 
  geom_boxplot(width=0.08, color="Black", fill="White", alpha=0.6, outlier.size = 0.1, show.legend=F, lwd=0.2) +
  scale_y_log10(name="4C contacts, TMM", breaks=c(0.001,0.01, 0.1, 1, 10, 100, 1000), labels=scales::number,
                minor_breaks = c(0.005,0.05,0.5,5,50,500,5000)) +
  stat_summary(fun.data = give.n, geom = "text", fun = median, position = position_dodge(width = 0.75), color="black") +
  scale_fill_npg()
#  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
#  scale_fill_manual(values=c("#56B4E9", "#E69F00")) 
invisible(dev.off())


# Alternative version with splitted violins
#
#suppressPackageStartupMessages(library(gghalves))
#ggplot(data=FinalPlot,aes(x=Condition, y=value, split=name, fill=name)) +
#  theme_bw() +
#  geom_half_violin(position="identity", alpha=0.6, lwd=0.4, width=2.8) + 
#  geom_boxplot(position = position_dodge(0.25), alpha=0.5, width=0.1, color="Black", fill="White", outlier.size = 0.1, show.legend=F, lwd=0.3) +
#  scale_y_log10(name="4C contacts, TMM", breaks=c(0.001,0.01, 0.1, 1, 10, 100, 1000, 10000, 100000), labels=scales::number,
#                minor_breaks = c(0.5,5,50,500,5000,50000)) +
#  stat_summary(fun.data = give.n, geom = "text", fun = median, position = position_dodge(width = 0.75), color="black") 
#invisible(dev.off())

#https://r-graph-gallery.com/267-reorder-a-variable-in-ggplot2.html

# p.adjust.method
# method to adjust p values for multiple comparisons. Used when pairwise comparisons are performed. 
# "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". 
# If you don't want to adjust the p value (not recommended), use p.adjust.method = "none"

df <- pairwise.wilcox.test(Tests$data, Tests$key, p.adjust.method = "holm")
write.table (df$p.value, file="MWtest.txt", sep="\t", row.names=T, col.names=T, quote = F)
