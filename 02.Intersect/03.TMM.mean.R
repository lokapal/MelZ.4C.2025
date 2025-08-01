#!/usr/local/bin/Rscript
# script to calculate TMM (Trimmed Mean of M-values) normalized values from counts.txt table
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Input:  counts.txt generated by featureCounts
# Output: 1. counts.byTMMs.txt     
#
# Dependency tools & libraries:
# 1. R
# 2. edgeR	https://bioconductor.org/packages/release/bioc/html/edgeR.html
# Common R libraries:
# 2. dplyr


suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(dplyr))

# we cannot guess the number of replicates, set it manually below
replicates <- 2

#counts a numeric matrix of raw feature counts.
#featureLength a numeric vector with feature lengths that can be obtained using biomaRt package.
#meanFragmentLength a numeric vector with mean fragment lengths, which can be calculated using the CollectInsertSizeMetrics(Picard) tool.

countsR <- as.data.frame(read.table("counts.txt", header=TRUE, sep='\t', skip = 1, row.names=NULL, stringsAsFactors = FALSE))

fcs            <- countsR #%>% arrange(Chr, Start, End)
#head(fcs)
row.names(fcs) <- fcs$Geneid
# exclude superfluous columns
fcs            <- fcs[, -c(1:6)]

totalexp       <- ncol(fcs)
totalgroups    <- totalexp/replicates

fcvector       <- vector()
for (i in 1:totalgroups) {
    for (j in 1:replicates) {
      fcvector <- append(fcvector, i)
                            }
                         }
#### in manual mode, make something like for two group of two replicates:
# fcvector <- c (1, 1, 2, 2)

fcgroup        <- factor(fcvector)

dge  <- DGEList ( counts=fcs, group = fcgroup )

dge  <- calcNormFactors(dge, method="TMM")
tmms <- as.data.frame(cpm(dge, log= FALSE))
tmms <- tmms %>% mutate (Symbol = rownames(tmms)) %>% 
                 select(Symbol, everything())

fcnames        <- vector()
namescounts    <- vector()
for (i in 2:ncol(tmms))  {
          colname <- colnames(tmms[i])
          namescounts <- append(namescounts, colname)
          colname <- gsub("(\\.[A-Za-z0-9]+)$", ".TMM", colname)
          colnames(tmms)[i] <- colname
          fcnames <- append (fcnames, colname)
                         }

final <- left_join(countsR, tmms, by = c("Geneid" = "Symbol"))

#mutate (!!sym(RPKname) := sprintf("%.4f", !!sym(RPKname)))

final <- final  %>%  mutate(MeanCounts = rowMeans(select(., all_of(namescounts)))) %>%
                     mutate(MeanTMMs  = rowMeans(select(., all_of(fcnames)))) %>%
# The select to skip unknown number of first columns except of
                       select (-c(all_of(namescounts), MeanCounts, all_of(fcnames), MeanTMMs), all_of(namescounts), MeanCounts, all_of(fcnames), MeanTMMs) %>%
#                       arrange(Chr, Start, End) %>%
                       arrange(desc(MeanTMMs)) %>%
                       filter(MeanCounts != 0)

outfile <- paste0("counts.byTMMs.txt")
write.table(final, file=outfile, row.names=F, col.names=T, sep="\t", quote=F)
