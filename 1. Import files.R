##### Analysis of CellProfiler results
##
## Requirements: 
## 1) meta = MyExpt_Image.txt
## Contains Image names
## Name: MyExpt_Image.txt = contains list of image conditions
## Columns: FileName_DAPI, ImageNumber
##
## 2) Fluorescent values
## Name: MyExpt_Nuclei.txt
## Columns: ImageNumber, ObjectNumber, Intensity_MeanIntensity_Ch00, 
## Intensity_MeanIntensity_Ch01, Intensity_MeanIntensity_Ch02, 
## Intensity_MeanIntensity_Ch03, 

### 1. Set working environment
setwd("...")

### 2. Import files
meta<-read.table("MyExpt_Image.txt", sep="\t", header=T)
raw<-read.table("MyExpt_Nuclei.txt", sep="\t", header=T)

### 3. Extract sample name from meta - FileName_DAPI
#3a. IF there is a constant number of characters
first_ch<-1       #first character
end_ch<-12        #last character from the end
meta[,1]<-as.character(meta[,1])
Sample_Name<-substr(meta[,1], first_ch, (nchar(meta[,1])-end_ch))    #extract name string
meta<-cbind(meta, Sample_Name)    #generate meta table

#3b. IF there is a regulat expression
end_ch<-16
Sample_Name<-sub('.*488_', '', meta[,1])
Sample_Name<-substr(Sample_Name, 1, (nchar(Sample_Name)-end_ch))    #extract name string
meta<-cbind(meta, Sample_Name)


### 4. Add sample name to raw data
raw_data<-merge(raw, meta, by="ImageNumber")
write.table(raw_data, "MyExpt_Nuclei-annotated.txt", sep="\t", quote=F)

