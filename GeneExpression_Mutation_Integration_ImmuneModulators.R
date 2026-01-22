
## clearing workspace

rm(list = ls())

## loading some initial packages for reading and subsetting files 

library(rapport)
library(maftools)


### file locations might be different depending on where and when this is being run 

## TCGA mutation MAF file 

D <- read.maf('/Users/iqbalw/Downloads/mc3.v0.2.8.PUBLIC.maf.likelyDriverLoose.tsv')


## TCGA sample Mutations to rmove
# in this case hyper mutator samples

hm <- read.delim("/Users/iqbalw/Documents/HyperMutatorSamples.txt")

## Gene list to use for mutations

#l <- read.csv('/Users/iqbalw/Downloads/20200803_ChromatinRemodellingGeneList_WI.csv')


#subset of genes from list above we are interested in 

list <- c("ARID1A","ARID2","BAP1","SMARCB1","SMARCA1","SMARCA4","PBRM1","KMT2A","KMT2D","KDM5C","KMT2C",
          "KMT2B","KDM6A","DNMT3A","CHD4","CHD3","CHD8","CREBBP","NCOR1","TBL1XR1","SIN3A","NIPBL", "BRD1", "BRD2", "BRD3", "BRD4", "BRD5", "BRD6", "BRD7", "BRD8")


## subset the maf file using maftools package

x1 <- subsetMaf(maf = D, genes = list, mafObj = FALSE, query = "Variant_Type  == 'SNP'")

## Seperate the long format, since sample names in list have different format
#otherwise code below wont work, and samples will not be rmoved

x1$Sample <- gsub('\\.','-',substr(x1$Tumor_Sample_Barcode,1,12))


#remove hyper mutant samples from Hm list
list <- hm$TCGA_Barcode

#keep only samples not in list
x <- x1[!(x1$Sample %in% list),]



#load the gene expression file 
#already correctly formatted previously

c <- readRDS('ExpressiondOrganized.rds')


#Removing columns with extra Id info, seperate parts of original barcode
c <- c[,5:ncol(c)]

#Adding Gene Annotation for Probe
pa <- readRDS(file = '/Users/iqbalw/Downloads/ImmuneProbeAnnotations.rds')

pa$Cancer_Type <- 'NA'
pa$formatted_ID <- 'NA'

#combine annotations, not counting the gene names row in c, otherwise will be repeated 
c <- rbind(pa,c[2:nrow(c),])

#Add 'NA' in cancer type column for annotation rows
c$Cancer_Type[1:5] <- 'NA'

#Removing samples from hyper mutated set

#removing hyper mutated samples from the methylation data

c$Sample <-  trimws(paste('TCGA',trimws(sapply(strsplit(rownames(c),"_"),"[", 1)),trimws(sapply(strsplit(rownames(c),"_"),"[", 2)), sep="-"))

c <-c[!(c$Sample %in% list),]

#Remove the last column c$sample, last column

c<- c[,1:(ncol(c) -1)]

#Changing Tumor Barcode in mutation file to match that of methylation file

slist <- as.character(x$Tumor_Sample_Barcode)


# format IDs for matching

x1a <- sapply(strsplit(slist,"-"), "[", 2) #Study code
x2a <- sapply(strsplit(slist, "-"), "[", 3) #Pateient ID
x3a <-sapply(strsplit(slist, "-"), "[", 5) # Sample Type

# want only the nuerical part and not the letter as in 11D
# other files contains Ds (DNA), this file contains Rs (RNA)

x4a <- sapply(strsplit(slist, "-"), "[", 4)
x4a <- as.numeric(gsub("([0-9]+).*$", "\\1", x4a))



x$Tumor_Sample_Barcode <- paste(x1a,"_", x2a, "_", x4a)

c2 <- as.data.frame(t(c[,1:(ncol(c) -2)]), make.names = FALSE)

cadd <- c[, (ncol(c) -1)] 


t <- c2[with(c2, order(geneNames)),]

t <- as.data.frame(t(t), make.names = FALSE)

c <- cbind(t,cadd)
names(c)[ncol(c)] <- 'Cancer_Type'

#Adding 'NA' to cancerType column for Annotation rows


clist <- unique(c$Cancer_Type[6:nrow(c)])
glist <- unique(x$Hugo_Symbol)


#create a dataframe to store
# rows is the number of probes we are working with
# columns is gene names associated with each Probe + # of modulator genes for which drivers were found
# columns is where rank sum test for immune and 

# to store data for differentialy expressed Gene data
edata <- data.frame(matrix(NA, nrow = 0, ncol = 10))
names(edata)[1] <- 'ImmuneGeneID'
names(edata)[2] <- 'MutatedChromatinModulatorGene'
names(edata)[3] <- 'CancerType'
names(edata)[4] <- 'P_value'
names(edata)[5] <- 'MutatedSamples'
names(edata)[6] <- 'NonMutatedSamples'
names(edata)[7] <- 'MedianExpression_Mutated'
names(edata)[8] <- 'MedianExpression_NonMutated'
names(edata)[9] <- 'MedianDifference_MutatedVs.Non'
names(edata)[10] <- 'Effect_on_Expression'



library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(DT)
library(dplyr)
library(reshape2)
library(gtools)
library(cowplot)
library(gridExtra)

CombinedData <- readRDS('/Users/iqbalw/DT/CombinedData2.rds')

c$rownames <- rownames(c)
tc <- rbind(c[1,],c[6:nrow(c),])
tc <- data.frame(t(tc))
tc <- data.frame(t(distinct(tc)))
rownames(tc) <- tc$rownames

#drop rownames column
tc <- tc[ , !names(tc) %in% "rownames"]

## Padding  by 1

# convert the number part to numerical
# ignore the first row which hase gene names
# ignore last column which has cancer info
ntc <- data.frame(apply(tc[2:nrow(tc), 1:(ncol(tc)-1)], 2, function(x) as.numeric(as.character(x))))
rownames(ntc) <- rownames(tc)[2:nrow(tc)]

# now add 1 to the converted part 
ntc <- ntc  + 1

ntc <- cbind(ntc, tc[2:nrow(tc), ncol(tc)] )
names(ntc)[ncol(ntc)] <- "Cancer_Type"
ntc <- rbind(tc[1,], ntc)

tc <- ntc
#probe list 
pl <- data.frame(t(c[1,1:(ncol(c)-1)]))

Upcount <- data.frame(matrix(NA, nrow = 23, ncol = 0))
Downcount <- data.frame(matrix(NA, nrow = 23, ncol = 0))

pa <- readRDS(file = '/Users/iqbalw/Downloads/ImmuneProbeAnnotations.rds')
plist <- cbind(data.frame(t(pa[1,])), data.frame(t(pa[5,])))

# Gene names with gene family
#missing gene family changed to 'Other" in gene family col ( col. 2)
# used later in code to reorder heatmap columns


newdata <- plist[1:(nrow(plist) -1),]
newdata <- distinct(newdata)
newdata[newdata$Gene.Family == "",2] <- paste('Other')

#Alphabeticallly so all family genes are together
newdata<- newdata[order(newdata$Gene.Family),]

# used for row annotation
annf <- data.frame(newdata$Gene.Family)
rownames(annf) <- newdata$geneNames

for (i in 1:length(clist)) {
  Mcount <- data.frame(matrix(NA, nrow = 1, ncol = length(glist)))
  
  add <- data.frame(matrix(NA, nrow = (ncol(tc) -1), ncol = (length(glist) + 1)))
  #subset by cancer type
  v <- tc[tc$Cancer_Type == clist[i],]
  colnames(v)[1:(ncol(v) - 1)] <- tc[1,1:(ncol(tc) -1)]
  names(add)[1] <- 'GeneName'
  print(i)
    # Compare within each Probe by  getting Methylation data for each Probe
        for (k in (1:(ncol(v)-1))){
          add[k,1] <- tc[1,k]
          rownames(add)[k] <- colnames(tc)[k]
          Pdata <- as.data.frame(v[,k])
          rownames(Pdata) <- rownames(v)
          names(Pdata) <- names(v)[k]
    # comparison with mutated and not mutated samples (one chromatin gene at a time)
                for (z in 1:length(glist)){
                      #subset mutation data for this gene
                      j <- x[x$Hugo_Symbol == glist[z],]
                      #list of mutated samples
                      mlist <- as.character(j$Tumor_Sample_Barcode)
                      #annotate Probe data column if samples mutated or not for this gene 
                      Pdata$Mutation <- ifelse(rownames(Pdata) %in% mlist,'Mutated','NotMutated')
                      # extract the methylation data for mutated and non mutated samples
                      muts <- Pdata[Pdata$Mutation == 'Mutated',]
                      muts <- as.numeric(muts[,1])
                      nmuts <- Pdata[Pdata$Mutation == 'NotMutated',]
                      nmuts <- as.numeric(nmuts[,1])
                      #rownames(add)[k] <- names(Pdata)[1]
                      #Add the correlation value, it can be that in certain cancer set,no mutated samples for this gene 
                      
                      mm <- median(muts)
                      mn <- median(nmuts)
                      
                      
                      diff <- mm - mn
                   
                      
                              if(length(muts) >= 5 && abs(diff) >= 3){
                                res <- wilcox.test(muts, nmuts)
                                colnames(add)[z+1] <- glist[z]
                                #add[k,z+1] <- res$p.value
                                # median mutationvalue
                                # want to avoid cases of nan, inf or -inf
                                # in the foldchange output 
                          
                                        if (res$p.value <= .05)
                                            {
                                              Mcount[1,z] <- paste(length(muts), "|", length(nmuts))
                                              names(Mcount)[z] <- glist[z]
                                              add[k,z+1] <- foldchange(mm, mn)
                                              
                                              
                                              temp <- data.frame(matrix(NA, nrow = 1, ncol = 10))
                                              names(temp) <- names(edata)
                                              
                                              temp$ImmuneGeneID <- c[1,k]
                                              temp$MutatedChromatinModulatorGene <- glist[z]
                                              temp$CancerType <- clist[i]
                                              temp$P_value <- res$p.value 
                                              temp$MutatedSamples <- length(muts)
                                              temp$NonMutatedSamples <- nrow(v) - length(muts)
                                              temp$MedianExpression_Mutated <- median(muts)
                                              temp$MedianExpression_NonMutated <- median(nmuts)
                                              temp$MedianDifference_MutatedVs.Non <- diff
                                              
                                              if (diff > 0) {
                                                temp$Effect_on_Expression <- 'Increased'
                                    
                                              } else{
                                                temp$Effect_on_Expression <- 'Decreased'
                                      
                                              }
                                              edata <- rbind(edata,temp)
                                              
                                            }
                                          else{}

                              }
                      # Dont count in mutations < 5
                              else {
                                colnames(add)[z+1] <- glist[z]
                              }
                      
                   
                }
          
        }


 
  #mdtemp <- mdata[mdata$CancerType == i,]
  #keep <- data.frame(cbind(mdata$ImmuneGeneID, mdata$MutatedChromatinModulatorGene, kdf$NumberofDifferentiallyMethylatedProbes))
  
  #names(keep) <- c('ImmuneGenes', 'Chr.Modulator', 'DifferentialProbes')
  #ndf <- dcast(keep,ImmuneGenes~Chr.Modulator)
  
  
  
  #Remove any columns with all NA
  #add <- distinct(add)
  edata <- distinct(edata)
  
  #temp add to 
  
  tadd <- add[,2:ncol(add)]
  tadd[is.na(tadd)] <- 0
  
  Uc <- data.frame((colSums(tadd > 0)))
  names(Uc)[ncol(Uc)] <- clist[i]
  
  Upcount <- cbind(Upcount,Uc)
  
  Dc <- data.frame((colSums(tadd < 0)))
  names(Dc)[ncol(Dc)] <- clist[i]
  
  Downcount <- cbind(Downcount,Dc)
  
  
  
  add <- data.frame(add[,colSums(is.na(add))<nrow(add)])
  names(add)[1] <- 'GeneName'
  
  
  add <- data.frame(add[order(add$GeneName),])
  names(add)[1] <- 'GeneName'
  
  #load summarry sheet to add numbers on heatmap

  if (ncol(add) > 2){
    y <- datatable(add, rownames = FALSE) %>%
     formatStyle(columns = names(add))
                #  background = styleInterval(c(0, 0.05), c("red", "red", "white")))
    
    DT::saveWidget(y, paste('/Users/iqbalw/DT/GeneExpression/FoldChange_',clist[i],'GeneExpression.html'))
    
    # add[,2:ncol(add)] <- ifelse(add[,2:ncol(add)]<=.05 & add[,2:ncol(add)] > 0 ,1, ifelse(add[,2:ncol(add)]<=.05 & add[,2:ncol(add)] < 0, -1,0 ))
    
    #also remove any column that have all 0s (no sig in any of the rows)
    # add <- add[, colSums(add != 0) > 0]
    
    #seperate out the numerical parts from annotation parts
    
    anns <- data.frame(add[,1])
    rownames(anns) <- rownames(add)
    names(anns) <- 'Mutated_ImmuneGenes'
    
    test <- data.matrix(add[,2:ncol(add)])
    colnames(test) <- colnames(add)[2:ncol(add)]
    rownames(test) <- add$GeneName
    
    colourCount = length(unique(anns$ImmuneGenes))
    getPalette = colorRampPalette(brewer.pal(12, "Paired"))
    
    
    myColors <- getPalette(colourCount)
    names(myColors) <- levels(anns$ImmuneGene)
    colScale <- scale_colour_manual(name = "Mutated_ImmuneGenes",values = myColors)
    
    #addind the number of differentialy methylated probes to heatmapp
    kdf <- CombinedData[CombinedData$CancerType == clist[i],]
    keep <- data.frame(cbind(kdf$ImmuneGeneID, kdf$MutatedChromatinModulatorGene, kdf$NumberofProbes_IncreasedMethylation, kdf$NumberofProbes_DecreasedMethylation))
    names(keep) <- c('ImmuneGenes', 'Chr.Modulator', 'IncreasedMeth.Probes', 'DecreasedMeth.Probes')
    ndf1 <- dcast(keep,ImmuneGenes~Chr.Modulator, value.var = 'IncreasedMeth.Probes' )
    ndf2 <- dcast(keep,ImmuneGenes~Chr.Modulator, value.var = 'DecreasedMeth.Probes' )
    test2 <- matrix(NA, nrow = nrow(test), ncol = ncol(test))
    colnames(test2) <- colnames(test)
    rownames(test2) <- add[,1]
    for (ii in 1:nrow(ndf1)){
      for (kk in 2:ncol(ndf1)){
           if(!is.na(ndf1[ii,kk])){
              col <- colnames(ndf1)[kk]
              row <- ndf1[ii,1]
              c2 <- which(colnames(test2)==col)
              r <- which(rownames(test2)==row)
              #total probes for the probe gene
              tp <- pl[pl$geneNames == row,]
              test2[r,c2] <- paste("\u2191",ndf1[ii,kk], "\u2193", ndf2[ii,kk], "|", length(tp), sep = " ")
              #print(paste(ii,kk,'Value:', ndf1[ii,kk]))
           }else{
          }
      }
    }
    
    test2[is.na(test2)] <- paste('0')
    
    #test2 = lapply(test2, function(x) {x[x == ""] <- 0})
    
    # get median and std dev. before converting Na to 0 

    
    #mt <- min(test[!is.na(test)])
    
    mt <- 0 - 2*sd(test[!is.na(test)])
    
    test[is.na(test)] <- 0
    Mcount <- data.frame(Mcount[,colSums(is.na(Mcount))<nrow(Mcount)])
  

    
    rownames(Mcount) <- ''
    
   
    # Reorder rows by gene family 
    
    test <- test[order(match(rownames(test),newdata$geneNames)),]
    test2 <- test2[order(match(rownames(test2),newdata$geneNames)),]
    
   #breaksList = seq(-val, val, by = 1)
    breaksList = seq(-4, 4, by = 1)
    u <- pheatmap(test, cluster_cols = FALSE, cluster_rows = FALSE, annotation_row =  annf, display_numbers = test2, cellheight = 25, cellwidth = 85, fontsize_number = 15, main = paste(clist[i],"(n = ",nrow(v),") Gene Expression Differences"), breaks = breaksList, show_colnames= TRUE, show_rownames= TRUE, fontsize = 20, legend = TRUE,, legend_labels = c('DownRegulated', 'UpRegulated'), height = 25, width = 30, annotation_colors = colScale, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), silent = TRUE) #+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) #filename = paste('/Users/iqbalw/DT/GeneExpression/FoldChange/HeatMap/',clist[i],'GeneExpressionHeatmap44.jpeg'))
    
    #jpeg(paste('/Users/iqbalw/DT/TEST/',clist[i],'t.jpeg'), height = 2500, width = 2000)
    #grid.arrange( u[[4]], tableGrob(Mcount), nrow = 2)
    
    
    z <- data.frame(table(plist$geneNames[1:nrow(plist)- 1]))
    rownames(z) <- NULL
    names(z) <- c('Gene', 'Total_Probes')
    
    j <- tableGrob(z, rows = NULL, theme=ttheme_default(base_size = 20))
    
    k <- tableGrob(Mcount, theme=ttheme_default(base_size = 25))

    g <- plot_grid(u[[4]], j, k, NULL, ncol = 2, axis = "b", rel_widths = c(8, 1), rel_heights = c(20, 1 ), align = "hv") 
    ggsave(file= paste('/Users/iqbalw/DT/TEST/',clist[i],'Test2.jpeg'), g, width = 30, height = 30, limitsize = TRUE, scale = TRUE)
    #dev.off()
    
    
    

    
#   }else if (ncol(add) ==  1 ) {
#     
#   }else if (ncol(add) == 2 && min(add[,2]) < .05 ){
#     y <- datatable(add, rownames = FALSE) %>%
#       formatStyle(columns = names(add))
#                   #background = styleInterval(c(0, 0.05), c("red", "red", "white")))
#     
#     DT::saveWidget(y, paste('/Users/iqbalw/DT/GeneExpression/FoldChange_',clist[i],'GeneExpression.html'))
#     
#     add[,2:ncol(add)] <- ifelse(add[,2:ncol(add)]<=.05 & add[,2:ncol(add)] > 0 ,1, ifelse(add[,2:ncol(add)]<=.05 & add[,2:ncol(add)] < 0, -1,0 ))
#     
#     anns <- data.frame(add[,1])
#     rownames(anns) <- rownames(add)
#     names(anns) <- 'Mutated ImmuneGenes'
#     
#     test <- data.matrix(add[,2:ncol(add)])
#     colnames(test) <- colnames(add)[2:ncol(add)]
#     rownames(test) <- rownames(add)
#     
# 
#     
#     
#     pheatmap(test, cluster_cols = FALSE, cluster_rows = FALSE, annotation_row =  anns, main = paste(clist[i],"(n = ",nrow(v),") Gene Expression Differences"), show_colnames= TRUE, show_rownames= FALSE, fontsize = 11, legend = TRUE,legend_breaks= c(-1,1), legend_labels = c('DownRegulated', 'UpRegulated'), height = 17, width = 15, filename = paste('/Users/iqbalw/DT/GeneExpression/FoldChange/HeatMap/',clist[i],'GeneExpressionHeatmap.jpeg'))
#     
   }else {}
 }

Upcount <- cbind(Upcount, z$Total_Probes)
up <- datatable(Upcount, rownames = TRUE) %>%
        formatStyle(columns = names(Upcount))

DT::saveWidget(up, paste('/Users/iqbalw/DT/TEST/Upcount.html'))


Downcount <- cbind(Downcount, z$Total_Probes)
down <- datatable(Downcount, rownames = TRUE) %>%
  formatStyle(columns = names(Downcount))

DT::saveWidget(down, paste('/Users/iqbalw/DT/TEST/Downcount.html'))



c <- readRDS('/Users/iqbalw/Downloads/FormattedMethylationData_ImmuneGenes')

saveRDS(edata, file = '/Users/iqbalw/DT/GeneExpression/FoldChangeDifferences.rds' )
