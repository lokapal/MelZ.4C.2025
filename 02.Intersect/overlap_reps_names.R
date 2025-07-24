#!/usr/local/bin/Rscript
# script to perform differential 4C analysis 
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Input:  two files like plastic.rep1.nodfam.bed.gz and plastic.rep2.nodfam.bed.gz
# Output: 1. plastic.nodfam.rep1.overlap.txt.gz     replicate1 that has intersections with replicate2
#         2. plastic.nodfam.rep2.overlap.txt.gz     replicate2 that has intersections with replicate1
#
# Dependency tools & libraries:
# 1. R
# Common R libraries:
# 2. dplyr, foreach, data.table

suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))

suffix = ".bed.gz"

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Chromosome must be supplied: chr1", call.=FALSE)
}

multiResultClass <- function(result1=NULL,result2=NULL)
{
  me <- list (
        result1 = result1,
        result2 = result2
             )
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}

prefix  <- args[1]
chrname <- args[2]

name1 <- paste0(prefix,"rep1.",chrname,suffix)
name2 <- paste0(prefix,"rep2.",chrname,suffix)

bam1 <- fread (file=name1, sep="\t", header=F, col.names=c("seqnames","start","end","qname","mapq","strand"))
bam2 <- fread (file=name2, sep="\t", header=F, col.names=c("seqnames","start","end","qname","mapq","strand"))

bam1 <- bam1 %>% mutate (qname = gsub("/\\d$","",qname))
bam2 <- bam2 %>% mutate (qname = gsub("/\\d$","",qname))

blist1 <- split(bam1, by="seqnames", flatten = TRUE)
blist2 <- split(bam2, by="seqnames", flatten = TRUE)

cat (paste(chrname, "Data loaded\n"))

names1 <- names(blist1)
names2 <- names(blist2)

matches <- match(names1, names2)
#Now matches can be used to subset
chrnames <- names1[matches > 0]

start.time <- Sys.time()

res <- foreach (i = seq_along(chrnames))  %do% { 
           currname <- chrnames[i]
           setkey(blist1[[currname]], start, end)
           setkey(blist2[[currname]], start, end)
           b1b2 <- foverlaps(blist1[[currname]], blist2[[currname]], by.x=c("start","end"), type="any", mult="first", nomatch=NULL) %>% select (i.qname)
           b2b1 <- foverlaps(blist2[[currname]], blist1[[currname]], by.x=c("start","end"), type="any", mult="first", nomatch=NULL) %>% select (i.qname)
           result <- multiResultClass()
           result$result1 <- b1b2
           result$result2 <- b2b1
           rm (b1b2, b2b1) #, lst2, lst1)
           cat(paste0(currname, " processed\n"))
           return(result)
                                               }

end.time <- Sys.time()
time.taken <- end.time - start.time
cat(paste("Processing took ", time.taken, "\n"))


start.time <- Sys.time()
res_rep1 <- list()
res_rep2 <- list()
for (i in 1:length(res)) {
     res_rep1[[i]] <- res[[i]]$result1
     res_rep2[[i]] <- res[[i]]$result2
                         }

res_rep1  <- res_rep1 %>% bind_rows()
res_rep2  <- res_rep2 %>% bind_rows()

end.time <- Sys.time()
time.taken <- end.time - start.time
cat(paste("Conversion took ", time.taken, "\n"))

outname <- paste0(prefix,chrname,".rep1.overlap.txt.gz")
fwrite(res_rep1, file=outname, row.names=F, col.names=F, sep="\t")
outname <- paste0(prefix,chrname,".rep2.overlap.txt.gz")
fwrite(res_rep2, file=outname, row.names=F, col.names=F, sep="\t")
