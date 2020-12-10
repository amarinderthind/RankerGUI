###### Example how to intersect the data using R sctipt
#In this example we are merging multiple (13) cancer data files together for the comparisons
###### Path for the Directory

setwd("./")

#### keep all the files in one foler

dat_files <- list.files("./", pattern='*.txt')

### you can trim the names of files 

trim_names <- lapply(basename(dat_files), gsub, pattern = "-Full_table.*", replacement = "")

#### load files and assign into dataFrames

library(data.table)
l <- lapply( dat_files, read.table, sep = '\t' ) ## loading files into l

for (i in 1:length(l)) {
  
  d  <- as.data.frame(l[[i]])  #d is temp variable
   ##extracting specified coulmns from the files ## these col names should be same in all the files 
  names(d) <- c("Gene_symbol","Gene_ID","Log2_FC","adj_P") ##col with gene/probe_id is must; either Log2_FC or pvalue/padj can be used for this analysis
  assign(trim_names[[i]], d)   #Assigning names back 
}


Common_genes <- Reduce(intersect, list(ACC$Gene_symbol, BRCA$Gene_symbol, BLCA$Gene_symbol,HNSC$Gene_symbol,KIRC$Gene_symbol,LAML$Gene_symbol,LGG$Gene_symbol, LIHC$Gene_symbol, LUAD$Gene_symbol, LUSC$Gene_symbol,PAAD$Gene_symbol, SKCM$Gene_symbol,THCA$Gene_symbol))

## extract logFC for common genes for each files and merge

for (i in 1:length(l)) {
  
  d  <- as.data.frame(l[[i]])  #d is temp variable
  names(d) <- c("Gene_symbol","Gene_ID","Median_tumor","Median_normal","Log2_FC","adj_P")
  d <- d[!duplicated(d[ , c("Gene_symbol")]),]
  d <- d[d$Gene_symbol %in% Common_genes,c("Gene_symbol","Log2_FC")]
  names(d) <- c("Gene_symbol",paste(trim_names[[i]],"Log2FC", sep = '_'))
  
  #assign(trim_names[[i]], d)
  assign(paste('d',i,sep = ''), d)
}


#final_data <- Reduce(function(x, y) merge(x, y,by='Gene_symbol', all=TRUE), list(ACC$Gene_symbol, BRCA$Gene_symbol, BLCA$Gene_symbol,HNSC$Gene_symbol,KIRC$Gene_symbol,LAML$Gene_symbol,LGG$Gene_symbol, LIHC$Gene_symbol, LUAD$Gene_symbol, LUSC$Gene_symbol,PAAD$Gene_symbol, SKCM$Gene_symbol,THCA$Gene_symbol))

final_data <- Reduce(function(x, y) merge(x, y,by='Gene_symbol', all=TRUE), list(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13))

write.table(final_data,"ALL_data.tsv", sep = '\t', row.names = F)
