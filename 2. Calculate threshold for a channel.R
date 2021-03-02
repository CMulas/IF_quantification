
library(plyr)
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)


#### 1. Determine threshold for one channel (e.g. Ch02)
#### run each time for one channel
#### requires: MyExpt_Nuclei-annotated.txt = output of 1. Import files.R
#### Column Names...ImageNumber, ObjectNumber, Intensity_MeanIntensity_Ch01, 
### Intensity_MeanIntensity_Ch02, FileName_DAPI, Sample_Name

raw2<-read.table("MyExpt_Nuclei-annotated.txt", sep="\t", header=T)
raw2$Sample_Name<-as.factor(raw2$Sample_Name)

# make smaller table, containing only positive and negatives (to calculate threshold)
# and the column of interests.
# 
positive<-"2i"
negative<-"Nanog"
channel<-"Intensity_MeanIntensity_Ch02"


subset_raw<-function(df, pos, neg, ch) {
  
  keep<-c("Sample_Name", channel)
  keep_columns<- paste(keep, collapse = "|")
  
  #clean table to include only positive and negative controls
  raw2_small<-rbind(df[(df$Sample_Name==pos),], df[(df$Sample_Name==neg),])
  raw2_small$Sample_Name<-droplevels(raw2_small$Sample_Name, exclude = if(anyNA(levels(raw2_small$Sample_Name))) NULL else NA)
  
  raw2_small2<-raw2_small %>% select(matches(keep_columns))
  colnames(raw2_small2)<-c("MeanIntensity", "Sample_Name")
  return(raw2_small2)
}

raw2_threshold<-subset_raw(raw2, positive, negative, channel)



####Receiver Operating Characteristic code
## from:https://www.r-bloggers.com/illustrated-guide-to-roc-and-auc/
# df = raw2_threshold; colnames "MeanIntensity", "Sample"
# cost_of_fp = ... 2= if false positive is less costly; 1=false postive is more costly
# cost_of_fn = ... 2= if false negative is less costly; 1=false negative is more costly
# PosSample = ... string of positive staining control (e.g. "2i")
# NegSample = ... string of negative staining control (e.g. "Nanog")
# n = 100 ... number of thresholds to calculate
# important! data range= 0-1, but it can be changed


calculate_roc<-function(df, cost_of_fp, cost_of_fn, PosSample, NegSample, n=200) {
  tpr <- function(df, threshold) {
    sum(df$MeanIntensity >= threshold & df$Sample_Name == PosSample) / sum(df$Sample_Name == PosSample)
  }
  
  fpr <- function(df, threshold) {
    sum(df$MeanIntensity >= threshold & df$Sample_Name == NegSample) / sum(df$Sample_Name == NegSample)
  }
  
  cost <- function(df, threshold, FP, FN) {
    sum(df$MeanIntensity >= threshold & df$Sample_Name == NegSample) * FP + 
      sum(df$MeanIntensity >= threshold & df$Sample_Name == PosSample) * FN
  }
  
  # important! data range default is 0,1, but any number can be used if the parameters below are changed
  # to change, substitute seq(x, y, length.out=n)
  roc <- data.frame(threshold = seq(0,1,length.out=n), tpr=NA, fpr=NA)
  roc$tpr <- sapply(roc$threshold, function(th) tpr(df, th))
  roc$fpr <- sapply(roc$threshold, function(th) fpr(df, th))
  roc$cost <- sapply(roc$threshold, function(th) cost(df, th, cost_of_fp, cost_of_fn))
  
  
  ##calculate distance to optimal 
  #the ideal threshold is the one closest to false positive rate of 0, and a true positive rate of 1
  
  library(raster)
  optimal <- c(0,1)
  FvsTPR<- cbind(roc$fpr, roc$tpr)
  roc$dist_to_min<- pointDistance(FvsTPR, optimal, lonlat=FALSE)
  
  return(roc)
}

# test parameters
# df = raw2_threshold
# cost_of_fp <- 2
# cost_of_fn <- 1
# PosSample <- positive ... = 2i
# Negative <- negative ... = Nanog
# n <- 100


roc_table<-calculate_roc(raw2_threshold, 1, 2, positive, negative, n=100)
table_title<-sprintf("ROC table for %s.txt", channel)
write.table(roc_table, table_title, sep="\t", quote=F)

new_threshold<- roc_table[roc_table$dist_to_min == min(roc_table$dist_to_min),1]


###plot ROC with new threshold
library(ggplot2)
plot_roc <- function(roc, threshold, cost_of_fp, cost_of_fn) {
  library(grid)
  library(gridExtra)
  
  norm_vec <- function(v) (v - min(v))/diff(range(v))
  
  idx_threshold = which.min(abs(roc$threshold-threshold))
  
  col_ramp <- colorRampPalette(c("green","orange","red","black"))(100)
  col_by_cost <- col_ramp[ceiling(norm_vec(roc$cost)*99)+1]
  p_roc <- ggplot(roc, aes(fpr,tpr)) + 
    geom_line(color=rgb(0,0,1,alpha=0.3)) +
    geom_point(color=col_by_cost, size=4, alpha=0.5) +
    coord_fixed() +
    geom_line(aes(threshold,threshold), color=rgb(0,0,1,alpha=0.5)) +
    labs(title = sprintf("ROC")) + xlab("FPR") + ylab("TPR") +
    geom_hline(yintercept=roc[idx_threshold,"tpr"], alpha=0.5, linetype="dashed") +
    geom_vline(xintercept=roc[idx_threshold,"fpr"], alpha=0.5, linetype="dashed")
  
  p_cost <- ggplot(roc, aes(threshold, cost)) +
    geom_line(color=rgb(0,0,1,alpha=0.3)) +
    geom_point(color=col_by_cost, size=4, alpha=0.5) +
    labs(title = sprintf("cost function")) +
    geom_vline(xintercept=threshold, alpha=0.5, linetype="dashed")
  
  sub_title <- sprintf("Threshold at %.2f - of channel %s", threshold, channel)
  
  grid.arrange(p_roc, p_cost, ncol=2, top=textGrob(sub_title, gp=gpar(cex=1), just="top"))
}

pdf_title<-sprintf("Threshold estimation_ROC and cost function for %s.pdf", channel)

pdf(pdf_title)
plot_roc(roc_table, threshold = new_threshold, cost_of_fp = cost_of_fp, cost_of_fn = cost_of_fn)
dev.off()



# plot distribution of positive and negative control after analysis
titles<-sprintf("Threshold at %.2f - for channel %s", new_threshold, channel)
pdf_title<-sprintf("Positive vs Negative controls for %s.pdf", channel)

pdf(pdf_title)
z<-ggplot(raw2_threshold, aes(x=MeanIntensity, fill=Sample_Name))+
  geom_density(alpha=0.3)+
  geom_vline(aes(xintercept=new_threshold, color="threshold"), linetype = "longdash")+
  ggtitle(titles)
z
dev.off()


# re-plot individual points
pdf_title<-sprintf("Positive vs Negative dotplot for %s.pdf", channel)

pdf(pdf_title)
w<-ggplot(raw2_threshold, aes(x=Sample_Name, y=MeanIntensity))+
  geom_violin(aes(color=Sample_Name), alpha=0.6)+
  geom_jitter(aes(color=Sample_Name), alpha=0.6)+
  geom_hline(yintercept=new_threshold, color="red", alpha=0.6)+
  labs(title=titles)
w
dev.off()


##plot density distributions for each sample for 'channel', red line=threshold

pdf_title<-sprintf("Distribution of %s.pdf", channel)

pdf(pdf_title)
qq<-ggplot(raw2_threshold, aes(x=MeanIntensity, fill=Sample_Name))+
  geom_density(alpha=0.3)+
  facet_grid(Sample_Name ~ ., scales="free")+
  xlim(0,0.4)+
  geom_vline(aes(xintercept=new_threshold, color="Threshold"), linetype = "longdash")+
  labs(title=titles)
qq
dev.off()



