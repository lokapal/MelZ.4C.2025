#!/usr/local/bin/Rscript
# script to chart circular genomes with top 150 4C-rDNA associated genes
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Input:  4C.genes.txt, presenting results of differential 4C analysis
#         genome markup in GTF (that was utilized for differential analysis): hg38.112.gtf.gz
# Output: 1. MelZ.circos.padj.pdf          MelZ plastic and matrigel 4C rDNA-associated 
#                                          top 150 (by padj) genes with the most significant 
#                                          difference in the number of contacts between 
#                                          plastic and matrigel
#
# Dependency tools & libraries:
# 1. R
# Bioconductor & CRAN packages:
# 2. circlize             https://cran.r-project.org/web/packages/circlize/index.html
# 3. plyranges            https://bioconductor.org/packages/release/bioc/html/plyranges.html
# Common R libraries:
# 4. dplyr, Polychrome
GTFile <- "hg38.112.gtf.gz"

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(Polychrome))
suppressPackageStartupMessages(library(plyranges))

cairo_pdf("MelZ.circos.padj.pdf",width=15,height=10,pointsize=16,antialias="default",fallback_resolution = 300,onefile=T)
circos.clear()
circos.par("start.degree" = 90)
par(cex = 1.0)
chromolist <- paste0("chr", c(1:22,'X')) 
circos.initializeWithIdeogram(species = "hg38", chromosome.index = chromolist, major.by=50000000, ideogram.height=0.05)

allgenes <- as.data.frame(read.table("4C.genes.txt", header = T, row.names=NULL, sep = "\t", stringsAsFactors = FALSE))
# allgenes <=> GeneID

gr <- read_gff(GTFile)
gr <- filter (gr, type == "gene" ) %>% select(gene_id, gene_name)
genetable<-as.data.frame(gr)
genetable<-genetable %>% select (seqnames, start, end, strand, gene_id) %>%
                         dplyr::rename (Chr=seqnames, Start=start, End=end)
#genetable<-genetable %>% rename (Chr=seqnames)
#                         rename (Chr=seqnames, Start=start, End=end)
#, gene_name)
# genetable <=> gene_id

allgenes  <- allgenes %>% dplyr::left_join(genetable, by = c("GeneID" = "gene_id"))

allgenes  <- subset(allgenes, padj <= 0.05)
allgenes  <- allgenes %>%  mutate(Matrigel.Mean = rowMeans(select(., Matrigel.rep1, Matrigel.rep2))) %>%
                           mutate(Plastic.Mean  = rowMeans(select(., Plastic.rep1, Plastic.rep2))) %>%
                           mutate(AbsLog = abs(log2FoldChange)) %>%
                           dplyr::rename (GeneName=gene_name)

allgenes  <- allgenes[order(allgenes$padj, -allgenes$AbsLog),]
selgenes  <- head(allgenes, 150)
head(selgenes)

PlasticGenes  <- selgenes %>% dplyr::rename (score=Plastic.Mean) %>%
                              dplyr::filter (score>0)                       
MatrigelGenes <- selgenes %>% dplyr::rename (score=Matrigel.Mean) %>%
                              dplyr::filter (score>0)            

#PlasticGenes <- subset(PlasticGenes, score > 0)
# or df %>% filter(column1=='A' | column2 > 8)
#MatrigelGenes <- subset(MatrigelGenes, score > 0)

PlastBed <- PlasticGenes %>% select(Chr,Start,End,GeneName)
PlastBed <- PlastBed[order(PlastBed$Chr, PlastBed$Start, PlastBed$End), ]
#nrow(PlastBed)

MatriBed <- MatrigelGenes %>% select(Chr,Start,End,GeneName)
MatriBed <- MatriBed[order(MatriBed$Chr, MatriBed$Start, MatriBed$End), ]

circos.genomicLabels(PlastBed,side="inside",labels.column=4,cex=0.5,line_lwd=0.25)

if(!exists(".Random.seed")) set.seed(NULL)
if (file.exists("seeddata")) {
    load("seeddata")
    .Random.seed <- x        } else {
     x <- .Random.seed
     save(x, file = "seeddata") 
     }
# create palette of N distinct colors
P23 = createPalette(length(chromolist), c("#ff0000", "#00ff00", "#0000ff"))

#chrom_colors <- setNames(adjustcolor(rainbow(length(chromolist)),alpha.f=0.5), chromolist)
chrom_colors <- setNames(adjustcolor(P23,alpha.f=0.5), chromolist)

#             col = chrom_colors[as.character(from_chr)]
#             link.lwd = Plastlinks$score / max(Plastlinks$score) * 5  # Adjust link width by score

     maxPlastic  <- max(PlasticGenes$score)
     maxMatrigel <- max(MatrigelGenes$score)
     multiplier  <- 4e7
     multiplier2 <- 2000

Plastlinks <- PlasticGenes %>% select(Chr,Start,End,score) %>%
     mutate( middl = (End - Start )/2.0 + Start,
#             widt  = 0.5 * score/maxPlastic * multiplier,
             widt  = score * multiplier2,
             ANCchr = 'chr14',
             Ancbeg = 0,
             Ancend = 43000,
             Nstart = middl - widt,
             Nend   = middl + widt )
     
Plastlinks$Nstart[Plastlinks$Nstart < 0] <- 0

bedG <- Plastlinks %>% select (Chr, Nstart, Nend) # %>% rename (chr=Chr, start=Nstart, end=Nend)
bedA <- Plastlinks %>% select (ANCchr, Ancbeg, Ancend) # %>% rename (chr=ANCchr, start=Ancbeg, end=Ancend)
#nrow(bedG)

circos.genomicLink(bedG, bedA, col = chrom_colors[as.character(bedG$Chr)], border = NA)

# second page
circos.clear()
circos.par("start.degree" = 90)
par(cex = 1.0)
circos.initializeWithIdeogram(species = "hg38", chromosome.index = chromolist, major.by=50000000, ideogram.height=0.05)
circos.genomicLabels(MatriBed,side="inside",labels.column=4,cex=0.5,line_lwd=0.25)

Matrilinks <- MatrigelGenes %>% select(Chr,Start,End,score) %>%
     mutate( middl = (End - Start )/2.0 + Start,
#             widt  = 0.5 * score/maxPlastic * multiplier,
             widt  = score * multiplier2,
             ANCchr = 'chr14',
             Ancbeg = 0,
             Ancend = 43000,
             Nstart = middl - widt,
             Nend   = middl + widt )
     
Matrilinks$Nstart[Matrilinks$Nstart < 0] <- 0

bedG <- Matrilinks %>% select (Chr, Nstart, Nend) # %>% rename (chr=Chr, start=Nstart, end=Nend)
bedA <- Matrilinks %>% select (ANCchr, Ancbeg, Ancend) # %>% rename (chr=ANCchr, start=Ancbeg, end=Ancend)
#nrow(bedG)

circos.genomicLink(bedG, bedA, col = chrom_colors[as.character(bedG$Chr)], border = NA)

invisible(dev.off())

