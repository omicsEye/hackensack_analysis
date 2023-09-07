#Note: importing BioC pkgs after dplyr requires explicitly using dplyr::select()
library(ggplot2)
library(stringr)
library(dplyr)
library(vroom)
library(Maaslin2)
setwd("~/Desktop")
### Gather all counts files
for(r in 2:79){
  args = c(r)
  myColNum=as.numeric(args[1])
  myColNum=myColNum+1
  myTable=vroom("mapping_clean.tsv")
  myTable=t(myTable)
  #counts.data=counts.data[,-as.data.frame(which(colSums(counts.data)==0))[,1]]
  humanReadPercents=vroom("humanAlignmentRateAcrossSamples.csv",col_names = c("sample","number","human","other"))
  humanReadPercents$sample=as.character(str_extract_all(humanReadPercents$sample,"NJ-BioR-[0-9]+"))
  humanReadPercents$human=as.numeric(str_remove_all(humanReadPercents$human,"%"))
  length(which(humanReadPercents$human>20))
  #counts.data2=myTable[,which(humanReadPercents$sample[which(humanReadPercents$human>20)] %in% colnames(myTable))]
  #hist(colSums(counts.data2),breaks=20)
  #counts.data3=counts.data2[,-which(colnames(counts.data2) %in% names(which(colSums(counts.data2)<1000000)))]
  counts.data=myTable
  #write.table(counts.data,file="counts.data.tsv",sep='\t',quote=FALSE)
  #row.names(counts.data)="viralLoad"
  metadata=vroom("cleanMetadataNonDiscretized.tsv")
  colnames(counts.data)=counts.data[1,]
  idx=match(metadata$ID,colnames(counts.data))
  metadata=metadata[which(!is.na(idx)),]
  idx2=match(colnames(counts.data),metadata$ID)
  counts.data=counts.data[,which(!is.na(idx2))]
  colnames(counts.data)=counts.data[1,]
  identical(colnames(counts.data),metadata$ID)
  metadata$severe="N"
  metadata$severe[which(metadata$ICU_Duration>0)]="Y"
  SEVERE=c(metadata$severe)
  counts.data[is.na(counts.data)] <- 0
  colnames(counts.data)=counts.data[1,]
  for(x in colnames(metadata)[myColNum]){
    print(x)
    metadata=as.data.frame(metadata)
    thisMeta=as.data.frame(metadata[-which(is.na(metadata[,which(colnames(metadata)==x)])),
                                    (which(colnames(metadata)==x))])
    #row.names(thisMeta)=metadata$ID[-which(is.na(metadata[,which(colnames(metadata)==x)]))]
    if(length(which(is.na(metadata[,which(colnames(metadata)==x)])))>0){
      input=as.data.frame(counts.data[,-which(is.na(metadata[,which(colnames(metadata)==x)]))])
      colnames(counts.data)=counts.data[1,]
      row.names(thisMeta)=metadata$ID[-which(is.na(metadata[,which(colnames(metadata)==x)]))]
    }else{
      input=as.data.frame(counts.data)
      colnames(input)=counts.data[1,]
      thisMeta=as.data.frame(metadata[,(which(colnames(metadata)==x))])
      row.names(thisMeta)=metadata$ID
    }
    
    colnames(thisMeta)=x
    #row.names(input)=myTable$Geneid
    write.table(input,file="input.tsv",sep='\t',quote=FALSE)
    write.table(thisMeta,file="thisMeta.tsv",sep='\t',quote=FALSE)
    print(counts.data[1,])
    input=input[-1,]
    input=as.data.frame(t(input))
    input$SARS_Count=as.numeric(input$SARS_Count)
    if(myColNum==2){
      Maaslin2(input_data=input,
               input_metadata=thisMeta,
               output=paste("Maaslin2_Viral_Load/","Maaslin2",x,"viralLoad",sep="_"),
               plot_heatmap = T,
               plot_scatter = T,
               standardize = F,
               reference="Test,Nasal")  
    }
    if(myColNum==9){
      Maaslin2(input_data=input,
               input_metadata=thisMeta,
               output=paste("Maaslin2_Viral_Load/","Maaslin2",x,"viralLoad",sep="_"),
               plot_heatmap = T,
               plot_scatter = T,
               standardize = F,
               min_prevalence = 0.1,
               reference="Race,White")  
    }
    if(myColNum==13){
      Maaslin2(input_data=input,
               input_metadata=thisMeta,
               output=paste("Maaslin2_Viral_Load/","Maaslin2",x,"viralLoad",sep="_"),
               max_significance = 0.05,
               plot_heatmap = T,
               plot_scatter = T,
               standardize = F,
               min_prevalence = 0.1,
               reference="Smoker,N")  
    }
    if(myColNum==31){
      Maaslin2(input_data=input,
               input_metadata=thisMeta,
               output=paste("Maaslin2_Viral_Load/","Maaslin2",x,"viralLoad",sep="_"),
               plot_heatmap = T,
               plot_scatter = T,
               standardize = F,
               min_prevalence = 0.1,
               reference="Cause_of_Death,COVID19")  
    }
    if(myColNum==78){
      thisMeta$Month[which(thisMeta$Month==3)]="March"
      thisMeta$Month[which(thisMeta$Month==4)]="April"
      thisMeta$Month[which(thisMeta$Month==5)]="May"
      print("Yes, hello")
      print(thisMeta$Month)
      
      Maaslin2(input_data=input,
               input_metadata=thisMeta,
               output=paste("Maaslin2_Viral_Load/","Maaslin2",x,"viralLoad",sep="_"),
               plot_heatmap = T,
               plot_scatter = T,
               standardize = F,
               min_prevalence = 0.0,
               reference="Month,March")
    }
    if(myColNum==79){
      Maaslin2(input_data=input,
               input_metadata=thisMeta,
               output=paste("Maaslin2_Viral_Load/","Maaslin2",x,"viralLoad",sep="_"),
               plot_heatmap = T,
               plot_scatter = T,
               standardize = F,
               min_prevalence = 0.0,
               reference="clade,B.1")  
    }else{
      Maaslin2(input_data=input,
               input_metadata=thisMeta,
               output=paste("Maaslin2_Viral_Load/","Maaslin2",x,"viralLoad",sep="_"),
               plot_heatmap = T,
               plot_scatter = T,
               standardize = F,
               min_prevalence = 0.0)
    }
  }
}
