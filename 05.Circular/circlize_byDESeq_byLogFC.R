#!/usr/local/bin/Rscript
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(Polychrome))
suppressPackageStartupMessages(library(plyranges))

cairo_pdf("circ1.pdf",width=15,height=10,pointsize=16,antialias="default",fallback_resolution = 300,onefile=T)
circos.clear()
circos.par("start.degree" = 90)
par(cex = 1.0)
chromolist <- paste0("chr", c(1:22,'X')) 
circos.initializeWithIdeogram(species = "hg38", chromosome.index = chromolist, major.by=50000000, ideogram.height=0.05)

allgenes <- as.data.frame(read.table("4C.genes.txt", header = T, row.names=NULL, sep = "\t", stringsAsFactors = FALSE))
# allgenes <=> GeneID

gr <- read_gff("hg38.112.gtf.gz")
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

allgenes  <- allgenes[order(-allgenes$AbsLog, allgenes$padj),]
selgenes  <- head(allgenes, 250)
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

chrom_colors <- setNames(adjustcolor(P23,alpha.f=0.5), chromolist)

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
             widt  = 0.5 * score/maxPlastic * multiplier,
#             widt  = score * multiplier2,
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
