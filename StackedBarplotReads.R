library(vroom)
library(ggplot2)
library(stringr)
library(cowplot)
library(dplyr)
setwd("C:/Users/tyson/OneDrive/Desktop/Coronavirus Proteins/UpdatedHackensack")
data=vroom("humanAlignmentRateAcrossSamples.csv",delim=",",col_names = c("Sample","Index","Human Percentage","Other Percentage"))
data$`Sample`=str_remove_all(data$Sample,"'")

data$`Human Percentage`=str_remove_all(data$`Human Percentage`,"%")
data$`Human Percentage`=as.numeric(data$`Human Percentage`)

data$`Other Percentage`=str_remove_all(data$`Other Percentage`,"%")
data$`Other Percentage`=as.numeric(data$`Other Percentage`)
newData=as.data.frame(rbind(cbind(data$Sample,"Human",data$`Human Percentage`),
                            cbind(data$Sample,"Other",data$`Other Percentage`)))
newData$V3=as.numeric(newData$V3)
readCountData=vroom("readCountsOfFiles.txt",col_names = c("Sample","Count"))
colnames(newData)=c("Sample","Source","Percent")
fullData=full_join(newData,readCountData,by="Sample")
fullData$Sample=as.character(str_extract_all(fullData$Sample,pattern="NJ-BioR-[0-9]{2,3}"))
readCount=as.data.frame(fullData[1:158,])
p1=ggplot(readCount,aes(x=Sample,y=Count))+
  geom_bar(stat="identity",fill="darkblue") +
  theme(axis.text.x = element_text(angle = 45,size=4, vjust = 0.1, hjust=1))
p2=ggplot(fullData, aes(fill=Source, y=Percent, x=Sample)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle = 45,size=4, vjust = 0.1, hjust=1), legend.position = "none")
pdf(file="ReadCountAndPercent.pdf",width=7.2,height=4)
plot_grid(p1,p2,nrow=2,ncol=1)
dev.off()
