

#Read Files

b <- readRDS('/mnt/isilon/zhou_lab/projects/20191212_GEO_datasets/GSE110554/5.22.2020_OpenSesameBetas.rds')
c <- readRDS('/mnt/isilon/zhou_lab/projects/20191212_GEO_datasets/GSE110554/samples_GEOLoadSeriesMatrix.rds')

#Getting the proper names
#Check to make sure sample names are the same
#otherwise rely on sorting

colnames(b) <- c$title

#Select only subset of data contianing th samples

df <- b[ , grepl( "Th" , colnames( b ) ) ]

#Remove NA

#transpose so that columns are probe names
#Rows are sample names

df <- t(df)

#Create Copy of df to work with bimodal index because dont want to lose it

df2 <- df

#remove NA from the data since bimodalIndex does not work with NAs
df2 <- na.omit(df2)

# Run BimodalIndex on the data
library(BimodalIndex)
bi <- bimodalIndex(df2, verbose=FALSE)

#Change the values of the beta file into 1, 0 or -1
#for methylated, in between or unmethylated




#Colums are Probes
#samples are rows

#change the data in original df (that has NA values also) based on results for bimodal Index

#Code Will not work with NA values
#Convert all NA to .5  (they will then be converted to 0 in last statement) or add in loop below

df1 <- t(df)
"5.22.2020_GettingThSampleMethylationConsensus" 109L, 2411C                                                                                                                               1,1           Top
{
  for(k in 1:ncol(df1))
  {
    if (df1[i,k] <= (bi[i,1] + 2*bi[i,3])) {
      df1[i,k] = -1
    } else if (df1[i,k] >= (bi[i,2] - 2*bi[i,3])) {
      df1[i,k] = 1
    } else {
      df1[i,k] = 0
    }
  }
}



##Get Consensus Vector

#Take and store column sum


s <- colSums(df1[,-1])

#add column names
colnames(s) <-colnames(df1)

#If 6/7 of the samples were methylated sum would be 6 -> Convert to 1 if sum 6 or more
#if sum -6 or less convert to 0
# else convert to -1

#In this case we will get 1 for methylated, 0 for unmethylated and -1 for all else (non relevant probes)

s <- colSums(df1[,-1])
for(k in 1:length(s))
{
  if (s[k] <= -6 ){
    s[k] = 0
  } else if (s[k] >= 6 ) {
    s[k] = 1
  } else {
    s[k] = 0
  }
}

#Save the binary file along with new conssus vector

saveRDS(s, file = '/mnt/isilon/zhou_lab/projects/20191212_GEO_datasets/GSE110554/5.22.2020_ThSampleConsensusVector.rds')
saveRDS(df1, file = '/mnt/isilon/zhou_lab/projects/20191212_GEO_datasets/GSE110554/5.22.2020_AllThSampleBinary.rds')
