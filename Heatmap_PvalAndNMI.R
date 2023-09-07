library(vroom)
setwd("C:/Users/tyson/OneDrive/Desktop/Coronavirus Proteins/UpdatedHackensack/omeClustFilteredMeta")
listOfRegions=list.files()
counter=1
listOfRegions=listOfRegions[-which(listOfRegions=="output")]
region=c()
for (x in listOfRegions){
  exists=list.files(path=x)
  if(length(exists)>2){
    currentRegion=as.data.frame(vroom(paste(x,"/","clusters.txt",sep="")))
    if(counter==1){
      associationTable=as.data.frame(currentRegion[1,])
    }else{
      currentRegion=currentRegion[,match(x=colnames(currentRegion), table=colnames(associationTable))]
      associationTable=as.data.frame(rbind(associationTable,currentRegion[1,]))
    }
    region=append(region,x)
    counter=counter+1
  }
}
rownames(associationTable)=region
#View(associationTable)
associationTable=associationTable[,c(-1,-2,-3,-4,-5)]
#View(associationTable)
categories=c("Clinical",
             "Outcome",
             "Clinical",
             "Outcome",
             "Clinical",
             "Clinical",
             "Clinical",
             "Clinical",
             "Clinical",
             "Clinical",
             "Clinical",
             "Clinical",
             "Outcome",
             "Clinical",
             "Patient",
             "Clinical",
             "Clinical",
             "Patient",
             "Outcome",
             "Outcome",
             "Clinical",
             "Clinical",
             "Other",
             "Patient",
             "Outcome",
             "Outcome",
             "Other",
             "Clinical",
             "Patient",
             "Outcome",
             "Clinical",
             "Patient",
             "Patient",
             "Other",
             "Patient",
             "Medicine",
             "Outcome",
             "Medicine",
             "Patient",
             "Patient",
             "Medicine",
             "Clinical",
             "Outcome",
             "Outcome",
             "Medicine",
             "Medicine",
             "Outcome",
             "Patient",
             "Medicine",
             "Medicine",
             "Patient",
             "Medicine",
             "Medicine",
             "Patient",
             "Medicine",
             "Medicine",
             "Patient",
             "Medicine",
             "Medicine",
             "Outcome",
             "Patient",
             "Medicine",
             "Patient",
             "Patient",
             "Medicine",
             "Medicine",
             "Medicine",
             "Patient",
             "Medicine",
             "Patient",
             "Outcome",
             "Medicine",
             "Medicine",
             "Patient",
             "Patient",
             "Patient",
             "Other",
             "Medicine")
library(ComplexHeatmap)
library(RColorBrewer)
#Heatmap(as.matrix(associationTable))
row.names(associationTable)=str_remove_all(row.names(associationTable),"_output")
row_dend = dendsort(hclust(dist(t(associationTable))))
#pdf(file="heatmap_cleanNames.pdf",width=30,height=8)
Heatmap(as.matrix(associationTable),
        column_split = data.frame(categories),
        col = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
        column_title_gp = gpar(fill = c("red", "blue", "green","orange","purple")),
        column_order = sort(colnames(associationTable)),
        cluster_columns = FALSE,
        cluster_rows = FALSE
)
#dev.off()
#write.csv(associationTable,file="associationTable.csv")

pval=vroom("PERMANOVA_pval.txt")
coef=vroom("PERMANOVA_coef.txt")
saveRowNames=coef$...1
saveRowNames2=pval$...1
identical(colnames(coef),colnames(pval))
matchIDX=match(colnames(associationTable),colnames(coef))
#matchIDX2=match(colnames(associationTable),colnames(pval))
coef=coef[,matchIDX]
pval=pval[,matchIDX]
#coef=coef[,-1]
row.names(coef)=saveRowNames
row.names(pval)=saveRowNames2
library(dplyr)
coef <- mutate_all(coef, function(x) as.numeric(as.character(x)))
coef[is.na(coef)] <- 0
pval <- mutate_all(pval, function(x) as.numeric(as.character(x)))
#pval[is.na(pval)] <- 0
library(dendsort)

#matPos=match(colnames(coef),heat1@column_names_param$labels)
#saveRowNames3=row.names(pval)
#pval=pval[,-1,drop=FALSE]
#row.names(pval)=saveRowNames3
#matPos=match(colnames(pval),heat1@column_names_param$labels)
matchIDX=match(tolower(row.names(associationTable)),tolower(row.names(pval)))
pval=pval[matchIDX,,drop=FALSE]
saveColNames=colnames(pval)
saveRowNames=row.names(pval)
saveColNames2=colnames(coef)
matchIDX2=match(tolower(colnames(associationTable)),tolower(colnames(pval)))
pval=pval[,matchIDX2,drop=FALSE]
colnames(pval)=saveColNames[matchIDX2]
row.names(pval)=saveRowNames2[matchIDX]
colnames(coef)=saveColNames2[matchIDX2]
row.names(pval)=saveRowNames2[matchIDX]
row.names(coef)=saveRowNames2[matchIDX]
#matchIDX=match(tolower(row.names(associationTable)),tolower(row.names(coef)))
coef=coef[matchIDX,,drop=FALSE]
saveColNames=colnames(coef)
#=match(tolower(colnames(associationTable)),tolower(colnames(coef)))
#coef=coef[,matchIDX,drop=FALSE]
colnames(coef)=saveColNames[matchIDX2]


important=rep("",dim(pval)[1]*dim(pval)[2])
important[which(pval<0.05)]="*"
important[which(pval<0.01)]="**"
important[which(pval<0.005)]="***"
important=matrix(important,nrow=24,ncol=78)

pdf(file="PERMANOVA_heatmap_PvalAndNMI.pdf",width=30,height=8)
Heatmap(as.matrix(associationTable),
        column_split = data.frame(categories),
        col = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
        column_title_gp = gpar(fill = c("red", "blue", "green","orange","purple")),
        cluster_rows = FALSE,
        column_order = sort(colnames(associationTable)),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", important[i, j]), x, y, gp = gpar(fontsize = 8,fontface="bold",col="white"))
        })
dev.off()


pdf(file="PERMANOVA_heatmap_Pval.pdf",width=30,height=8)
Heatmap(as.matrix(pval),
        column_split = data.frame(categories),
        col = rev(colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100)),
        column_title_gp = gpar(fill = c("red", "blue", "green","orange","purple")),
        cluster_rows = FALSE,
        column_order = sort(colnames(associationTable)),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", important[i, j]), x, y, gp = gpar(fontsize = 8,fontface="bold",col="white"))
        })
dev.off()