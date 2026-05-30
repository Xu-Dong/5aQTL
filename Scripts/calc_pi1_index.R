library(optparse)
library(qvalue)
library(data.table)
library(dplyr)


option_list <- list(
        make_option(c("-t","--tissue"),type="character",default="NA",action="store",help="specify an input file")
)
opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))


# define functions
bootstrap_pi1 <- function(pvalues, pi0.method="bootstrap"){

  lambda <- seq(0.05, 0.9, 0.05)

  pvalues <- na.omit(pvalues)

  pi1_hat <- (1 - pi0est(pvalues, lambda = lambda, pi0.method = pi0.method)$pi0)

  return(pi1_hat)
}

leadqtl<-readRDS("output/cis_5aQTL.FDR05.lead.merge.rds")


allqtl <- fread(paste0(allpath, opt$tissue,"_all_pairs.txt"),header=T,sep="\t")
allqtl$id <- paste0(allqtl$gene_id,"_",allqtl$variant_id)

for(i in 2:length(leadqtl)){
        pvalues <- leadqtl[,c(1,..i)]
        pvalues <- na.omit(pvalues)
        message(i,": ",names(leadqtl)[i]," No.leadQTL ", nrow(pvalues))
        if(names(leadqtl)[i] != opt$tissue){
                idx_select <- allqtl[id %in% pvalues$id]
                tmp <- paste0("./output/",names(leadqtl)[i],"_.in._",opt$tissue,".rds")
                saveRDS(idx_select,tmp)
                pvalues1 <-idx_select$pvalue
                message(opt$tissue," No.allQTL ", length(pvalues1))
                fit <- try(pi1_bootstrap <- bootstrap_pi1(pvalues1, pi0.method = "bootstrap"))
                if("try-error" %in% class(fit))
                {
                        next
                }else
                {
                        fwrite( data.table(leadqtl = names(leadqtl)[i], allqtl=opt$tissue, pi_value=pi1_bootstrap), paste("./output/",opt$tissue,".sharing.pi.txt"), row.names = FALSE, quote = FALSE, append = TRUE, sep="\t")
                        message("leadqtl: ", names(leadqtl)[i], " & allqtl: ", opt$tissue, "=", pi1_bootstrap," has been writen")
                }
        }
}
