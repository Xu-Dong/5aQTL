library(ggplot2)
library(mashr)

rm(list=ls())

source("./bin/normfuncs.R")
thresh <- 0.05 #0.05
hetthresh <- 0.5 #0.5
nsigthresh <- 0 #0

# snv_eQTL ----------------------------------------------------------------

# thresh <- 0.05
out <- readRDS("./fastqtl_to_mash_output/eQTL_file_list.mash.rds")
maxbeta <- out$strong.b
maxz    <- out$strong.z
cname <- colnames(maxbeta) %>% gsub(".txt", "", .)


lfsr <- read.table("./eQTL.lfsr.txt", header = T)

out2 <- readRDS("./mashr_flashr_workflow_output/eQTL.nominal.mash.final.results.rds")

pm.mash <- out2[["result"]][["PosteriorMean"]] # Z-score
standard.error <- maxbeta/maxz
pm.mash.beta   <- pm.mash * standard.error

#Histograms of sharing by magnitude
sigmat <- (lfsr <= thresh)
nsig   <- rowSums(sigmat)
c <- hist((het.func(het.norm(effectsize=pm.mash.beta[nsig>nsigthresh,]),threshold=hetthresh)),main="",xlab="",ylab="Proportion of nteQTL",
          breaks=0.5:49.5,col=rgb(1,0.5,0.5,alpha=0.8),freq=FALSE, xaxt="n",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,border = NA)
axis(1,at = seq(1, 49, by=1),labels = c(1:49))
mtext("All Tissues")
maglist <- het.func(het.norm(effectsize=pm.mash.beta[nsig>nsigthresh,]), threshold=hetthresh)
new=data.frame("eQTL",c$mids,c$density)
colnames(new)=c("Group",	"Number",	"Density")
new <- new %>% mutate(type="Magnitude")

p1 <- ggplot(new,aes(x=Number,y=Density)) + geom_bar(stat = "identity",width=.9,fill="grey") + theme_classic();p1


#Histogram of sharing by sign
sign.func <- function (normeffectsize){
  apply(normeffectsize,1,function(x)(sum(x>0)))}
sigmat <- (lfsr<=thresh)
nsig   <- rowSums(sigmat)
siglist <- sign.func(het.norm(effectsize=pm.mash.beta[nsig>nsigthresh,]))
a <- hist(sign.func(het.norm(effectsize=pm.mash.beta[nsig>nsigthresh,])),main="",xlab="",
          breaks=0.5:49.5,col="black",freq=FALSE,xaxt="n",ylim=c(0,1)) 
newsiglist <- sign.func(het.norm(effectsize=pm.mash.beta[nsig>nsigthresh,]))
axis(1, at=seq(1, 49, by=1), labels=c(1:49))
mtext("All Tissues")
newsig=data.frame("eQTL",a$mids,a$density)
colnames(newsig)=c("Group",	"Number",	"Density")
newsig <- newsig %>% mutate(type="Sign")

p2 <- ggplot(newsig,aes(x=Number,y=Density)) + geom_bar(stat = "identity",width=.9,fill="grey") + theme_classic();p2
df_all <- rbind(new, newsig)
saveRDS(df_all,file = "Sharing.hist.5aQTL.RDS")

library(cowplot)
cowplot::plot_grid(p1,p2,ncol=2,align = "h")
