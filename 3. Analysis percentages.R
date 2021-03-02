#### Analyse fluorescent data based on thresholds, determined automatically or manually
#### input1: raw_data = "MyExpt_Nuclei-annotated.txt"
####
# ImageNumber ObjectNumber Intensity_MeanIntensity_Ch01 Intensity_MeanIntensity_Ch02   FileName_DAPI Sample_Name
# 1           1            1                   0.09778267                   0.08431373 -ve_1a_ch00.tif         -ve
# 2           1            2                   0.16577860                   0.24153698 -ve_1a_ch00.tif         -ve
# 3           1            3                   0.17214422                   0.22188489 -ve_1a_ch00.tif         -ve
# 4           1            4                   0.15087458                   0.20857766 -ve_1a_ch00.tif         -ve
# 5           1            5                   0.16284182                   0.23585686 -ve_1a_ch00.tif         -ve
# 6           1            6                   0.17007807                   0.24846027 -ve_1a_ch00.tif         -ve
### input2: threshold_Nanog, threshold_Esrrb, threshold_Ch03, threshold_Ch04, etc... calculated from 2. Calculate threshold for a channel.R

#### 1. score each value for being above or below threshold

AboveThreshold<-function(df, channel, threshold) {
  column<-match(channel, names(df))
  as.numeric(df[,column]>=threshold)
}

raw_data$Ch01_positive<-AboveThreshold(raw_data, "Intensity_MeanIntensity_Ch01", threshold_Esrrb)
raw_data$Ch02_positive<-AboveThreshold(raw_data, "Intensity_MeanIntensity_Ch02", threshold_Nanog)
# raw_data$Ch03_positive<-AboveThreshold(raw_data, "Intensity_MeanIntensity_Ch03", threshold_Ch03)
# raw_data$Ch04_positive<-AboveThreshold(raw_data, "Intensity_MeanIntensity_Ch04", threshold_Ch04)


#### 2. count number of positive cells per condition, per image

summary_data<-ddply(raw_data, .(Sample_Name, ImageNumber), summarize, 
                    positives.Ch01 = 100*sum(Ch01_positive)/length(Ch01_positive),
                    positives.Ch02 = 100*sum(Ch02_positive)/length(Ch01_positive),
#                    positives.Ch03 = 100*sum(Ch03_positive)/length(Ch01_positive),
#                    positives.Ch04 = 100*sum(Ch04_positive)/length(Ch01_positive),
                    totals = length(Ch01_positive))

write.table(summary_data, "summary_data.txt", sep="\t", quote=F)


#### 3. Reshape table and calculate average per conditions

summarise_channel<-function (sum_table, pos_channel) {

  column<-match(pos_channel, names(sum_table))
  
  library(tidyverse)

  small_ch<-as.data.frame(cbind(sum_table$Sample_Name, sum_table[,column]))
  small_ch$V1<-as.factor(sum_table$Sample_Name)

  #all samples must have equal number of images to form a table. Add extra NAs to 

  sum_ch<-small_ch %>% 
   group_by_at(vars(-V2)) %>%
   dplyr::mutate(row_id=1:n()) %>% ungroup() %>%
   spread(key=V1, value=V2, fill=NA) %>%
   dplyr::select(-row_id)

  summ<-t(sum_ch)
  
  #calculate mean and standard deviation
  mu<-rowMeans(summ, na.rm=TRUE)
  sd<-apply(summ, 1, function(x) {sd(x, na.rm=TRUE)})

  summary_ch<-cbind(summ, mu, sd)
  
  table_name<-sprintf("Summary of %s.txt", pos_channel)
  write.table(summary_ch, table_name, sep="\t", quote=F)

  return(summary_ch)
}

Ch01_final<-summarise_channel(summary_data, "positives.Ch01")
Ch02_final<-summarise_channel(summary_data, "positives.Ch02")
# Ch03_final<-summarise_channel(summary_data, "positives.Ch03")
# Ch04_final<-summarise_channel(summary_data, "positives.Ch04")


#### 4. Plot bar chart for each condition + dot plot on top 
# input: summary_data (tidy)

#current sample order
t(levels(summary_data$Sample_Name))

sample_order<-c(2,1,5,4,6,3)

levels(summary_data$Sample_Name)[sample_order]

plot_bar_dot<-function(df, order){
  
  library(ggplot2)
  library(plyr)
  library(forcats)
  
  #make smaller table and re-order based on 'order' string
  dff<- df %>% 
    filter(Sample_Name %in% levels(summary_data$Sample_Name)[order])
  
  #generate summary data
  ord<-as.vector(levels(dff$Sample_Name)[order])
  sum_data<-ddply(dff, .(Sample_Name), summarize, 
                  mean = mean(get(channel)), sd = sd(get(channel)))
 
  sum_data$Sample_Name<-factor(sum_data$Sample_Name, ord)
  
  limits <- aes(ymax=mean+sd, ymin=mean-sd)

  pdftitle<-sprintf("Percentage %s cells.pdf", channel)
  
  gg<-ggplot(sum_data, aes(x=Sample_Name, y=mean, fill=Sample_Name))+
    geom_bar(stat="identity", width=0.5)+
    geom_errorbar(limits, width=0.25)+
    geom_jitter(data=dff, aes(x=Sample_Name, y=get(channel)), 
                 size=2,width = 0.1)+
    ggtitle(sprintf("Percentage of %s across conditions", channel))+
    ylab(sprintf("Percentage of %s", channel))+
    xlab("Condition")+
    theme_bw()+
   theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=10),
         axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
         panel.border = element_blank(), panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.line = element_line(colour = "black")) +
    scale_y_continuous(expand = c(0, 0))
  ggsave(pdftitle)
  gg
  }

channel<-"positives.Ch01"
plot_bar_dot(summary_data,  sample_order)

channel<-"positives.Ch02"
plot_bar_dot(summary_data,  sample_order)

# channel<-"positives.Ch03"
# plot_bar_dot(summary_data,  sample_order)
# 
# channel<-"positives.Ch04"
# plot_bar_dot(summary_data,  sample_order)




