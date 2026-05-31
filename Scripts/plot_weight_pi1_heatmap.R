setwd("/home/workspace/")
rm(list=ls())
df_tissues <- read.table("./input/GTEx_49tissues.txt",header=F,sep="\t")
tissueList <- unique(df_tissues$V1)
rm(df_tissues)

df_full <- expand.grid(
  discovery = tissueList,
  replicates = tissueList,
  stringsAsFactors = FALSE
)

dat <- fread("./output/all.sharing.pi.lambda0.9.txt",header=T)
dat$name <- NULL
dat$mid <- paste(dat$leadqtl,dat$allqtl,sep=":")


df_full$mid <- paste(df_full$discovery,df_full$replicates,sep = ":")
df_full$replicates <- NULL
df_full$discovery <- NULL
df_full$pi <- NULL

df_full <- merge(df_full,dat,by="mid",all.x = T)


disc_tissue <- function(x){
  return(strsplit(x,split=":",fixed=T)[[1]][1])
}


rep_tissue <- function(x){
  return(strsplit(x,split=":",fixed=T)[[1]][2])
}

df_full$discovery <- sapply(df_full$mid,disc_tissue)
df_full$replicate <- sapply(df_full$mid,rep_tissue)

df_plot <- df_full %>% select(discovery,replicate,pro_pi)
df_plot$pro_pi[df_plot$discovery==df_plot$replicate] <- 1

library(pheatmap)


weight_matrix <- dcast(df_plot, replicate ~ discovery, value.var = "pro_pi")
rownames(weight_matrix) <- weight_matrix$replicate
weight_matrix$replicate <- NULL
weight_matrix <- as.matrix(weight_matrix)

tmp_dat <- weight_matrix
tmp_mat[is.na(tmp_mat)] <- 0
hc <- hclust(
  dist(tmp_dat),
  method = "average"
)

pheatmap(
  weight_matrix,
  
  cluster_rows = hc,
  cluster_cols = hc,
  
  color = colorRampPalette(
    c("#4575B4",
      "#F7F7F7",
      "#FEE08B",
      "#D73027")
  )(100),
  
  na_col = "grey90"
)
