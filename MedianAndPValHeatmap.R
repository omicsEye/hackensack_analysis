library(vroom)
library(Maaslin2)
library(Tweedieverse)
library(dplyr)
library(data.table)
library(stringr)
library(ComplexHeatmap)
library(robustbase)
library(circlize)
library(Rfast)
#Read in the files
setwd("C:/Users/tyson/OneDrive/Desktop/Coronavirus Proteins/UpdatedHackensack/GMGC")
geneCounts=vroom("CombinedTable.txt")

metadata=vroom("cleanMetadataNonDiscretized.tsv")
#identical(gt$Sample,metadata$Patient)

saveRowNames=geneCounts$`...1`
geneCounts=geneCounts[,-1]
row.names(geneCounts)=saveRowNames
setwd("C:/Users/tyson/OneDrive/Desktop/Coronavirus Proteins/UpdatedHackensack/GMGC/enrichment_stats/")
list_of_enrichment_stats=list.files(pattern=".*enrichment*.")
counter=1
for(x in list_of_enrichment_stats){
  name=str_remove_all(x,"Deepath_BP_bigger_db_")
  name=str_remove_all(name,"_enrichment_stats.tsv")
  assign(paste(name,"_enrichment_stats",sep=""),vroom(x))
  assign((paste(name,"_enrichment_stats",sep="")),as.data.frame(cbind(name,get(paste(name,"_enrichment_stats",sep="")))))
  thisTable=get((paste(name,"_enrichment_stats",sep="")))
  thisTable=thisTable[,c(1,3,4,5,6,7,8,9)]
  colnames(thisTable)[1]="condition"
  if(counter==1){
    fullTable=thisTable
  }else{
    fullTable=as.data.frame(rbind(fullTable,thisTable))
  }
  counter=counter+1
}
setwd("C:/Users/tyson/OneDrive/Desktop/Coronavirus Proteins/UpdatedHackensack/GMGC/stats_table/")

list_of_stats_table=list.files(pattern=".*stats_table.*")
for(x in list_of_stats_table){
  name=str_remove_all(x,"Deepath_BP_bigger_db_")
  name=str_remove_all(name,"_stats_table.tsv")
  assign(paste(name,"_stats_table",sep=""),vroom(x))
}
neutrophil=which(str_detect(fullTable$pathway,"neutrophil"))
leukocyte=which(str_detect(fullTable$pathway,"leukocyte"))
myeloid=which(str_detect(fullTable$pathway,"myeloid"))
immune=which(str_detect(fullTable$pathway,"immune"))
fullTable=fullTable[sort(unique(c(neutrophil,leukocyte,myeloid,immune))),]
#clusters=as.data.frame(rbind(cluster1,cluster2,cluster3,cluster4))
#fullClusters=clusters
clusters=fullTable
#clusters=clusters[which(clusters$pval<0.00000005),]
#clusters=clusters[which(clusters$pval<0.00005),]
#clusters=clusters[which(clusters$pval<0.05),]
pValFrame=data.frame(matrix(ncol = length(unique(clusters$pathway)), nrow = length(list_of_enrichment_stats)))
medianFrame=data.frame(matrix(ncol = length(unique(clusters$pathway)), nrow = length(list_of_enrichment_stats)))
colnames(pValFrame)=unique(clusters$pathway)
row.names(pValFrame)=list_of_enrichment_stats
x_counter=1
newColNames=c()
newRowNames=c()
for(x in colnames(pValFrame)){
  newColNames=append(newColNames,x)
  thisPathway=fullTable[which(fullTable$pathway==x),]
  y_counter=1
  for(y in row.names(pValFrame)){
    newRowNames=append(newRowNames,y)
    if(x=="immune response" & y=="Deepath_BP_bigger_db_Ventilator_Duration_Deepath_Ventilator_Duration_enrichment_stats.tsv"){
     print("here")
    }
    if(x=="neutrophil activation" & y=="Deepath_BP_bigger_db_Creatinine_Deepath_Creatinine_enrichment_stats.tsv"){
      print("test")
    }
    if(y=="Deepath_BP_bigger_db_Cause_of_Death_Deepath_Cardiac (Heart Failure, MI, Myocarditis)_enrichment_stats.tsv"){
      cleanName="Deepath_BP_bigger_db_Cause_of_Death_Deepath_Cardiac_enrichment_stats.tsv"
      listOfGenes=unlist(str_split(thisPathway$pathway_members[which(row.names(pValFrame)==y)],";"))
      whichStatsTable=list_of_stats_table[which(str_detect(list_of_stats_table,str_remove_all(cleanName,"_enrichment_stats.tsv")))]
    }else if(y=="Deepath_BP_bigger_db_Test_Deepath_Throat_enrichment_stats.tsv"){
      #cleanName="Deepath_BP_bigger_db_Cause_of_Death_Deepath_Cardiac_enrichment_stats.tsv"
      #listOfGenes=unlist(str_split(thisPathway$pathway_members[which(row.names(pValFrame)==y)],";"))
      whichStatsTable="Deepath_BP_bigger_db_Test_Deepath_Throat_stats_table.tsv"
    }else if(y=="Deepath_BP_bigger_db_Test_Deepath_Throat+Nasal_enrichment_stats.tsv"){
      #cleanName="Deepath_BP_bigger_db_Cause_of_Death_Deepath_Cardiac_enrichment_stats.tsv"
      #listOfGenes=unlist(str_split(thisPathway$pathway_members[which(row.names(pValFrame)==y)],";"))
      whichStatsTable="Deepath_BP_bigger_db_Test_Deepath_Throat+Nasal_stats_table.tsv"
    }
    else{
      listOfGenes=unlist(str_split(thisPathway$pathway_members[which(row.names(pValFrame)==y)],";"))
      whichStatsTable=list_of_stats_table[which(str_detect(list_of_stats_table,str_remove_all(y,"_enrichment_stats.tsv")))]
    }
    name=str_remove_all(whichStatsTable,"Deepath_BP_bigger_db_")
    name=str_remove_all(name,"_stats_table.tsv")
    myStatsTable=get(paste(name,"_stats_table",sep=""))
    myStatsTable=myStatsTable[which(myStatsTable$feature %in% listOfGenes),]
    medianFrame[y_counter,x_counter]=median(myStatsTable$coef)
    #y=str_extract_all(y,"chr[0-9X]*_[0-9]*")
    name=str_remove_all(y,"Deepath_BP_bigger_db_")
    name=str_remove_all(name,"_enrichment_stats.tsv")
    pValFrame[y_counter,x_counter]=fullTable[which(fullTable$pathway==x & fullTable$condition==name),]$pval
    y_counter=y_counter+1
  }
  x_counter=x_counter+1
  print(x_counter)
}
important=rep("",dim(pValFrame)[1]*dim(pValFrame)[2])
important[which(pValFrame<0.05)]="*"
important[which(pValFrame<0.005)]="**"
important[which(pValFrame<0.001)]="***"
important=matrix(important,nrow=dim(pValFrame)[1],ncol=dim(pValFrame)[2])
colnames(medianFrame)=newColNames
row.names(medianFrame)=row.names(pValFrame)
#colMins(as.matrix(medianFrame),value=TRUE)
#colsToKeep=unique(c(which(colMins(as.matrix(medianFrame),value=TRUE)>0.1),which(colMaxs(as.matrix(medianFrame),value=TRUE)>0.1)))
name=str_remove_all(row.names(medianFrame),"Deepath_BP_bigger_db_")
name=str_remove_all(name,"_enrichment_stats.tsv")
row.names(medianFrame)=name
setwd("C:/Users/tyson/OneDrive/Desktop/Coronavirus Proteins/UpdatedHackensack/GMGC/")
neutrophil=which(str_detect(colnames(medianFrame),"neutrophil"))
leukocyte=which(str_detect(colnames(medianFrame),"leukocyte"))
myeloid=which(str_detect(colnames(medianFrame),"myeloid"))
immune=which(str_detect(colnames(medianFrame),"immune"))
colsToKeep=sort(unique(c(neutrophil,leukocyte,myeloid,immune)))
rowsToKeep=c()
for(x in c("Admission_To_Diagnosis_Time",
           "Admission_To_ICU_Transfer_Time",
           "Cough",
           "Hospitalization_Duration",
           "Diagnosis_To_Discharge_Time",
           "Discharged",
           "Ventilator_Duration",
           "Fever",
           "Shortness_of_Breath",
           "Age",
           "BMI",
           "Sex",
           "Weight",
           "Diabetes",
           "Height",
           "Hypertension",
           "Month",
           "acetaminophen",
           "Azithromycin",
           "ALT",
           "AST",
           "Neutrophil_Count",
           "Bilirubin",
           "CRP",
           "BNP",
           "Creatinine",
           "Ddimer",
           "Ddimer_correction",
           "E_Cq",
           "Elevated_LFTs",
           "IL_6_Admission",
           "N2_Cq",
           "RP_Cq",
           "troponin"
           )){
  if(length(which(str_detect(row.names(medianFrame),x)))<1){
    print(x)
  }
  rowsToKeep=append(rowsToKeep,(which(str_detect(row.names(medianFrame),x))))
}
important=important[rowsToKeep,]
library(ComplexHeatmap)
library(circlize)
#pdf(file="MedianAndPVal_PVal0.00000005.pdf",height=50,width=100)
pdf(file="GMGC_Immune_Processes_Heatmap.pdf",height=50,width=100)
Heatmap(medianFrame[rowsToKeep,], rect_gp = gpar(col = "grey93", lwd = 1),
        col = colorRamp2(c(-2,0,2),c("darkblue", "grey90", "darkred")),cluster_rows = TRUE, cluster_columns = TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", important[i, j]), x, y, gp = gpar(fontsize = 60,fontface="bold",col="white"),vjust=0.7)
        },
        width=unit(dim(medianFrame[rowsToKeep,])[2],"in"),
        height=unit(dim(medianFrame[rowsToKeep,])[1]/2,"in"),
        row_names_gp = gpar(fontsize = 22),
        column_names_gp = gpar(fontsize = 44),
        #,
        heatmap_legend_param = list(
          at = c(-2,0,2),
          labels = c("-2","0","2"),
          title = "Set Enrichment Score"
        ),
        column_names_rot=45)
dev.off()
# x_counter=1
# for(x in unique(fullTable$condition)){
#   thisCondition=fullTable[which(fullTable$condition==x),]
#   y_counter=1
#   pValFrame=data.frame(matrix(ncol = 10, nrow = 1))
#   medianFrame=data.frame(matrix(ncol = 10, nrow = 1))
#   for(y in thisCondition$pathway[1:10]){
#     thisPathway=fullTable[which(fullTable$pathway==y),]
#     listOfGenes=unlist(str_split(thisPathway$pathway_members[x_counter],";"))
#     whichStatsTable=list_of_stats_table[x_counter]
#     myStatsTable=get(paste(str_extract_all(whichStatsTable,"chr[0-9X]*_[0-9]*"),"_stats_table",sep=""))
#     myStatsTable=myStatsTable[which(myStatsTable$feature %in% listOfGenes),]
#     medianFrame[1,y_counter]=median(myStatsTable$coef)
#     x=str_extract_all(x,"chr[0-9X]*_[0-9]*")
#     pValFrame[1,y_counter]=fullTable[which(fullTable$pathway==y & fullTable$condition==x),]$pval
#     coefTable=as.data.frame(myStatsTable[which(myStatsTable$feature %in% listOfGenes),])
#     coefTable=coefTable[,c(4,5)]
#     coefTable=as.data.frame(rbind(coefTable,myStatsTable[which(myStatsTable$feature %in% listOfGenes),c(4,5)]))
#     p=ggplot(coefTable, aes(x=coef, fill=value))+
#       geom_vline(xintercept = 0, linetype=2,size=2) +
#       #geom_text(aes(x=0.01, label="coef = 0", y=20), colour="black", size=5, fontface="bold", angle=45) +
#       geom_density(alpha=0.4) + theme(axis.text.x=element_blank(),
#                                       axis.ticks=element_blank(),
#                                       axis.line=element_blank(),
#                                       axis.text.y=element_blank(),
#                                       axis.title.x=element_blank(),
#                                       axis.title.y=element_blank(),legend.position="none",
#                                       panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#                                       panel.grid.minor=element_blank(),plot.background=element_blank())
#     assign(paste("p",y_counter,sep=""),p)
#     y_counter=y_counter+1
#   }
#   important=rep("",dim(pValFrame)[1]*dim(pValFrame)[2])
#   important[which(pValFrame<0.01)]="*"
#   important[which(pValFrame<0.005)]="**"
#   important[which(pValFrame<0.001)]="***"
#   important=matrix(important,nrow=1,ncol=10)
#   colnames(medianFrame)=thisCondition$pathway[1:10]
#   row.names(medianFrame)=x[[1]]
#   start=round(min(medianFrame),digits=2)
#   if(start>0){
#     start=-0.01
#   }
#   end=round(max(medianFrame),digits=2)
#   if(end<0){
#     end=0.01
#   }
#   myHeatmap=Heatmap(medianFrame, rect_gp = gpar(col = "grey93", lwd = 1),
#                     col = colorRamp2(c(start,0,end),c("darkblue", "grey90", "darkred")),cluster_rows = FALSE, cluster_columns = FALSE,
#                     cell_fun = function(j, i, x, y, width, height, fill) {
#                       grid.text(sprintf("%s", important[i, j]), x, y, gp = gpar(fontsize = 60,fontface="bold",col="white"),vjust=0.7)
#                     },
#                     width=unit(dim(medianFrame)[2],"in"),
#                     height=unit(dim(medianFrame)[1]/2,"in"),
#                     row_names_gp = gpar(fontsize = 22),
#                     column_names_gp = gpar(fontsize = 44),
#                     #,
#                     heatmap_legend_param = list(
#                       at = c(start,0,end),
#                       labels = c(as.character(start),"0",as.character(end)),
#                       title = "Set Enrichment Score"
#                     ))
#   pdf(file=paste("MedianAndPValClusters","_",x[[1]],".pdf",sep=""),height=40,width=20)
#   print(myHeatmap)
#   dev.off()
#   myCow=cowplot::plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,nrow = 1)
#   pdf(file=paste("Densities","_",x[[1]],".pdf",sep=""),height=5,width=50)
#   print(myCow)
#   dev.off()
#   x_counter=x_counter+1
# }