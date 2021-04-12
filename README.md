# IF_quantification
Pipelines to quantify and analysed immunofluorescence images

## 1. Segment and quantify fluorescent images with CellProfiler
Approach: identify primary objects on DAPI channel, then measure object intensity across all channels. 
The exact parameters used to segment images for further quantification vary every time. However, I have uploaded a test pipelines for 10x images and quantification of nuclear signal. Adaptive Otsu thresolding worked better in my hands. The following parameters are particularly important:
1. Typical diameter of objects, in pixel units (Min, Max) = 12, 30 [small variations can affect nuclei merging/splitting]
2. Lower and upper bounds on threshold = 0.003, 1 [ to set the lower bound, enter test mode and use cursor to explore pixel intensities around nuclei.]

Tips for accurate quantification across experiments:
1. Extensive pre-analysis and image filtering can be done. However, images are typically pre-screened to avoid very uneven background, as the quantification is not very accurate in those cases. Ideally, modify acquision settings to avoid uneven background, or ensure the background is regular across all images.
2. Many microscopes have different focal planes depending on the fluorescent channel imaged. Use autofocus or refocus on each channel manually. 
3. Ensure the channel is not saturated and adjust imaging parameters to obtain a good dynamic range (ideally 12bit +).
  
Output files required for quantification: Select measurments to export 
1. All/Image/FileName/Ch00 (or any other channel) = "MyExpt_Image.txt"
2. All/Nuclei/Intensity/MeanIntensity/(all channels of interest) = "MyExpt_Nuclei.txt"
  
## 2.  Use ROC curves to calculate threshold.
I use a simple version of ROC curves in order to determine the automatic theshold for immonofluorescent images. This allows me to quantify across conditions and experiments the percentage of positive cells in a given channel (=immunostaining) to determine if a particular treatment has an effect or not.

### Data formating
File: "1. Import files.R"

Description: Adds file names to results tables

Input: "MyExpt_Nuclei.txt" and "MyExpt_Image.txt"

Output: "MyExpt_Nuclei-annotated.txt"

    ImageNumber= uniquely identifies each image

    ObjectNumber= identifies each object(=nuclei) within an image

    Intensity_MeanIntensity_Ch01 and Intensity_MeanIntensity_Ch01 = mean intensity of object (=nuclei) for a given channel (Ch01 or Ch02). Channels are immunofluorescence for a particular protein.

    FileName_DAPI = unique image ID, matches image names to look back if anything looks weird

    Sample_Name = Name of sample

<img width="788" alt="image" src="https://user-images.githubusercontent.com/61800079/109635160-f648f800-7b41-11eb-8799-a76fa5446081.png">

### ROC code
File: "2. Calculate threshold for a channel.R"
The code is adapted from this blog: https://www.r-bloggers.com/illustrated-guide-to-roc-and-auc/

You need to provide the program with:

    positive<-"2i"    = an element of Sample_Name (=i.e. what samples should be positive for the marker).

    negative<-"Nanog"   = an element of Sample_Name (=i.e. what samples should be negative for the marker).

    channel<-"Intensity_MeanIntensity_Ch02" = identifies a column name (=i.e. what channel you are comparing). 

### Output
You will generate an ROC table (saved as "ROC table for [channel]") and a threshold value (new_threshold). 

![ROC output](https://user-images.githubusercontent.com/61800079/109636458-850a4480-7b43-11eb-9d04-4b160cfdb80b.png)


The rest of the code generates diagnostic plots to make sure the the threshold makes sense.
![Threshold values](https://user-images.githubusercontent.com/61800079/109636542-9d7a5f00-7b43-11eb-93a2-9f9f97d51b72.png)

## 3. Analysis of percentage of positive cells across conditions
File "3. Analysis percentages.R" allows quantification of percentages of positive cells, to compare across experiments, calculate statistics, etc. Can merge multiple replicates at this stage.
