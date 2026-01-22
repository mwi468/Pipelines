rm(list = ls())

library(R.utils)
library(data.table)
library(maftools)
library(stringr)
library(rapport)




# load gene file 
genes <- read.table('/Users/iqbalw/Downloads/gencode.v35.basic.annotation.gtf ', header = FALSE, sep = '\t')


# seperate the gene names from the string in column in 9

genes$GeneNames <- paste(str_match(genes$V9, "gene_name \\s*(.*?)\\s*;")[,2])


# extract only the gene names that are part of the immune genne list 

#immune gene list

ilist <- read.csv ('/Users/iqbalw/thorsson_immunomodulators_list2.csv')


# getting ref. info  for onloy the immune genes

ig <- ilist[1:78,2]
ig <- as.list(ig)


genes <- genes[genes$GeneNames %in% ig,]

# getiing only the full coordinates  
# dont want individual codon info 

genes <- genes[genes$V3 == 'gene',]

# extending the TSS by 1500 depending on stran orientation

# for + strand want to subtract 1500 from TSS ( start for +)

genes_pos <- genes[genes$V7 == '+',]
genes_pos$V4 <- paste(genes_pos$V4 - 15)


# for - strand  add 1500 to TSS ( which would be the end in the file)

genes_neg <- genes[genes$V7 == '-',]
genes_neg$V5 <- paste(genes_neg$V5 + 15)



genes <- rbind(genes_pos, genes_neg)




# Formatting it so we Chr, Start, End and Probe_ID in the reference file 

genes <- cbind(genes$V1,genes$V4, genes$V5, genes$GeneNames)
write.table(genes, file="/Users/iqbalw/Downloads/Extended_ImmuneGenes_HG38ref3.bed", quote=F, sep="\t", row.names=F, col.names=F)

# fromatting ataseq data

ataq <- readRDS(file = '/Users/iqbalw/ATACseq.RDS' )
#write.table(ataqs, file="/Users/iqbalw/Downloads/ataq.bed", quote=F, sep="\t", row.names=F, col.names=F)




#bed tools intersect in terminal
#bedtools intersect -a /Users/iqbalw/Downloads/ataq.bed -b /Users/iqbalw/Downloads/Extended_ImmuneGenes_HG38ref2.bed -wo -bed >> Reference_Matched_ataq2.bed


#ataqs <- read.table(file="/Users/iqbalw/Reference_Matched_ataq3.bed", sep="\t")
ataqs <- as.data.frame(read.table(file="/Users/iqbalw/Reference_Matched_ataq3.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

#remove row with NA in this case the last one 

#ataqs <- ataqs[-326, ] 

names(ataqs)[1:(ncol(ataqs)-5)] <- names(ataq)
names(ataqs)[(ncol(ataqs)-4):ncol(ataqs)] <- c('Chr','Start','End','Gene_Name','Overlap')

#write.table(x, file="/Users/iqbalw/Downloads/ataqIntersect.bed", quote=F, sep="\t", row.names=F, col.names=F)
rm(ataq)


# load mutation samples and prep fie


D <- read.maf('/Users/iqbalw/Downloads/mc3.v0.2.8.PUBLIC.maf.likelyDriverLoose.tsv')


#subsetMaf(maf = laml, genes = c('DNMT3A', 'NPM1'), mafObj = FALSE, query = "Variant_Classification == 'Splice_Site'", fields = 'i_transcript_name')


l <- read.csv('/Users/iqbalw/Downloads/20200803_ChromatinRemodellingGeneList_WI.csv')

list <- c("ARID1A","ARID2","BAP1","SMARCB1","SMARCA1","SMARCA4","PBRM1","KMT2A","KMT2D","KDM5C","KMT2C",
          "KMT2B","KDM6A","DNMT3A","CHD4","CHD3","CHD8","CREBBP","NCOR1","TBL1XR1","SIN3A","NIPBL")


hm <- read.delim("/Users/iqbalw/Documents/HyperMutatorSamples.txt")

x1 <- subsetMaf(maf = D, genes = list, mafObj = FALSE, query = "Variant_Type  == 'SNP'")

x1$Sample <- gsub('\\.','-',substr(x1$Tumor_Sample_Barcode,1,12))


#remove hyper mutant samples from Hm list
list <- hm$TCGA_Barcode

x <- x1[!(x1$Sample %in% list),]



mlist <- x$Sample

#mutant ids

t <- read.csv(file = '/Users/iqbalw/ATACseqmetadata.txt', sep = '\t')

mi <- t[(t$submitter_id %in% mlist),]

t <- x[(x$Sample %in% mi$submitter_id),]

#save 

saveRDS(file = '/Users/iqbalw/ATACseqSampleMutationsin ChromModulators')
