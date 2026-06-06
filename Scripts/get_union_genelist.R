library(dplyr)
library(magrittr)
library(data.table)

rm(list=ls())

file_list <- list.files(path="./input/gene_list/",pattern="5aQTL_tested_gene.txt")

gene.list <- list()
tissue_list <- c()
for(i in 1:length(file_list)){
  tissueName <- strsplit(file_list[i],split = ".",fixed = T)[[1]][1]
  cat(i,": ",tissueName,"\n")
  tissue_list[i] <- tissueName
  dat <- fread(paste0("./input/gene_list/",file_list[i]),header = F,sep = "\t")
  gene.list[[i]] <- dat$V1
}

names(gene.list) <- tissue_list

common_genes <- Reduce(intersect, gene.list)

write.table(data.frame(phenotype_id=common_genes),file = "./input/Shared_5aQTL_genes.txt",quote=F,sep = "\t",row.names=F,col.names = F)
