#Do any of the clinical characteristics
#or microbial diversity 
#associate with sequencing metrics like total reads 
#or percentage of reads mapping to host?
library(vroom)
library(stringr)
library(ggcorrplot)
library(ggplot2)
library(GGally)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
readCounts=vroom("/Users/tyden46/Desktop/GWU/ReadCountCalculations/readCountsOfFiles.txt",col_names = c("sample","reads"))
humanCounts=vroom("/Users/tyden46/Desktop/GWU/ReadCountCalculations/humanAlignmentRateAcrossSamples.csv",col_names=c("sample","index","human","other"))
myData=vroom("/Users/tyden46/Desktop/GWU/Hackensack/cleanMetadataNonDiscretized.tsv")
humanCounts$sample=str_remove_all(humanCounts$sample,"'")
identical(readCounts$sample,humanCounts$sample)
humanCounts$totalReads=readCounts$reads
rm(readCounts)
humanCounts$sample=sapply(str_split(humanCounts$sample,"_"),"[",1)
myData=myData[order(myData$ID),]
humanCounts=humanCounts[order(humanCounts$sample),]
humanCounts$humanDecimal=as.numeric(sub("%", "",humanCounts$human,fixed=TRUE))/100
humanCounts=humanCounts[which(humanCounts$sample %in% myData$ID),]
myData=myData[which(myData$ID %in% humanCounts$sample),]
identical(myData$ID,humanCounts$sample)

myData$readCounts=humanCounts$totalReads
myData$humanReads=humanCounts$totalReads*humanCounts$humanDecimal
# Subset the dataframe to include only the relevant columns
myData[sapply(myData, is.character)] <- lapply(myData[sapply(myData, is.character)], 
                                               as.factor)
subsetData <- myData[, -c(ncol(myData)-1, ncol(myData))]
# Identify the numeric columns



# Visualize the correlation using a scatter plot matrix
#ggpairs(subsetData[, numericColumns], title = "Correlation Scatter Plot Matrix")

# Calculate the correlation coefficients
#correlationMatrix <- cor(subsetData)


subsetData=subsetData[,-36] #No variation in this column
factorColumns <- sapply(subsetData, is.factor)
numericColumns <- sapply(subsetData, is.numeric)
namesOfFactors=c()
pvalues=c()
for(x in which(factorColumns)){
  print(x)
  thisTable=as.data.frame(cbind(subsetData[,x],myData$readCounts))
  colnames(thisTable)=c(colnames(subsetData)[x],"readCounts")
  anova_result <- aov(as.formula(paste("readCounts","~",colnames(subsetData)[x])),thisTable)
  if(length(summary(anova_result)[[1]][["Pr(>F)"]])){
      P_Value = summary(anova_result)[[1]][["Pr(>F)"]][1]
  }else{
    P_Value=1
  }
  namesOfFactors=append(namesOfFactors,colnames(subsetData)[x])
  pvalues=append(pvalues,P_Value)
}


namesOfFactorsHumanReads=c()
pvaluesHumanReads=c()
for(x in which(factorColumns)){
  print(x)
  thisTable=as.data.frame(cbind(subsetData[,x],myData$humanReads))
  colnames(thisTable)=c(colnames(subsetData)[x],"humanReads")
  anova_result <- aov(as.formula(paste("humanReads","~",colnames(subsetData)[x])),thisTable)
  if(length(summary(anova_result)[[1]][["Pr(>F)"]])){
    P_Value = summary(anova_result)[[1]][["Pr(>F)"]][1]
  }else{
    P_Value=1
  }
  namesOfFactorsHumanReads=append(namesOfFactorsHumanReads,colnames(subsetData)[x])
  pvaluesHumanReads=append(pvaluesHumanReads,P_Value)
}
savePHumanAnova=pvaluesHumanReads
savePAnova=pvalues
anovaFrame=as.data.frame(rbind(pvalues,pvaluesHumanReads))
colnames(anovaFrame)=namesOfFactors

#pdf(file="CategoricalPValues.pdf",width=7.2,height=3)
Heatmap(anovaFrame,cluster_rows = F,cluster_columns = F,
        height = nrow(anovaFrame)*unit(5, "mm"),
        column_names_gp = gpar(fontsize = 8),
        col = colorRamp2(c(0,0.05,1),c("darkblue", "white", "red")))
#dev.off()

### Numerical

namesOfFactors=c()
pvalues=c()
for(x in which(numericColumns)){
  print(x)
  thisTable=as.data.frame(cbind(subsetData[,x],myData$readCounts))
  colnames(thisTable)=c(colnames(subsetData)[x],"readCounts")
  #anova_result <- cor.test(thisTable[,1],thisTable[,2])
  P_Value=cor.test(thisTable[,1],thisTable[,2])[3]$p.value
  namesOfFactors=append(namesOfFactors,colnames(subsetData)[x])
  pvalues=append(pvalues,P_Value)
}


namesOfFactorsHumanReads=c()
pvaluesHumanReads=c()
for(x in which(numericColumns)){
  print(x)
  thisTable=as.data.frame(cbind(subsetData[,x],myData$humanReads))
  colnames(thisTable)=c(colnames(subsetData)[x],"readCounts")
  #anova_result <- cor.test(thisTable[,1],thisTable[,2])
  P_Value=cor.test(thisTable[,1],thisTable[,2])[3]$p.value
  namesOfFactorsHumanReads=append(namesOfFactorsHumanReads,colnames(subsetData)[x])
  pvaluesHumanReads=append(pvaluesHumanReads,P_Value)
}

corFrame=as.data.frame(rbind(pvalues,pvaluesHumanReads))
colnames(corFrame)=namesOfFactors

#pdf(file="NumericalPValues.pdf",width=7.2,height=3)
Heatmap(corFrame,cluster_rows = F,cluster_columns = F,
        height = nrow(corFrame)*unit(5, "mm"),
        column_names_gp = gpar(fontsize = 8),
        col = colorRamp2(c(0,0.05,1),c("darkblue", "white", "red")))
#dev.off()



