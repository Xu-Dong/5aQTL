library(data.table)
library(dplyr)

setwd("./worksplace/")

qtl_file_list <- list.files("./input",pattern=".cis_5aQTL.FDR05.lead.txt")

for(i in 1:length(qtl_file_list)){
        qtl <- qtl_file_list[i]
        tissue <- strsplit(qtl,split=".",fixed=T)[[1]][1]
        qtl_dat <- fread(paste0("input/",qtl),header=T,sep="\t")
        names(qtl_dat)<-c("id","pvalue")
        leadqtl <- merge(leadqtl,qtl_dat,by="id",all.x=TRUE)
        names(leadqtl)[i+1] <- tissue
}


saveRDS(leadqtl,"./output/cis_atQTL.FDR05.lead.merge.rds")
