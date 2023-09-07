#Note: importing BioC pkgs after dplyr requires explicitly using dplyr::select()
library(ggplot2)
library(stringr)
library(dplyr)
library(vroom)
library(Maaslin2)
library(deepath)
### Gather all counts files
setwd("C:/Users/tyson/OneDrive/Desktop/Coronavirus Proteins/UpdatedHackensack/GMGC")
myTable=vroom("CombinedTable.txt")
saveRowNames=myTable$`...1`
myTable=myTable[,c(2:28)]
row.names(myTable)=saveRowNames
counts.data <- mutate_all(myTable, function(x) as.numeric(as.character(x)))
#counts.data=counts.data[,-as.data.frame(which(colSums(counts.data)==0))[,1]]
#colnames(counts.data)=str_replace_all(colnames(counts.data),"\\.","-")
counts.data[is.na(counts.data)]=0
colnames(counts.data)=substr(colnames(counts.data),1,11)
metadata=vroom("C:/Users/tyson/Box/HMH-CDI COVID Data/data/metadata/cleanMetadataNonDiscretized.tsv")
idx=match(metadata$ID,colnames(counts.data))
metadata=metadata[which(!is.na(idx)),]
idx2=match(colnames(counts.data),metadata$ID)
counts.data=counts.data[,which(!is.na(idx2))]
identical(colnames(counts.data),metadata$ID)
metadata$severe="N"
metadata$severe[which(metadata$ICU_Duration>0)]="Y"
SEVERE=c(metadata$severe)
#counts.data[is.na(counts.data)] <- 0
#write.csv(counts.data,file="CountsWithComparisons.csv")
#my_data <- list()
#counter=1
#for (i in colnames(metadata)) {
#  my_data[[i]] <- as.data.frame(metadata[,which(colnames(metadata)==i)])
#  counter=counter+1
#}
#for (x in colnames(metadata)[29:length(colnames(metadata))]){
metadata=as.data.frame(metadata)
row.names(counts.data)=saveRowNames
for(x in colnames(metadata)[34]){
  if(length(unique(pull(metadata,x)))==2 & x!="Test"){
    print(paste("There are two levels for",x))
    if(length(which(is.na(metadata[,which(colnames(metadata)==x)])))>0){
      input=as.data.frame(counts.data[,-which(is.na(metadata[,which(colnames(metadata)==x)]))])
      thisMeta=as.data.frame(metadata[-which(is.na(metadata[,which(colnames(metadata)==x)])),
                                      (which(colnames(metadata)==x))])
      row.names(thisMeta)=metadata$ID[-which(is.na(metadata[,which(colnames(metadata)==x)]))]
      
    }else{
      input=as.data.frame(counts.data)
      thisMeta=as.data.frame(metadata[,(which(colnames(metadata)==x))])
      row.names(thisMeta)=metadata$ID
    }
    colnames(thisMeta)=x
    #row.names(input)=row.names(counts.data)
    Maaslin2(input_data=input,
            input_metadata=thisMeta,
            output=paste("Tweedieverse",x,"GMGC",sep="_"),
            max_significance = 0.05,
            plot_heatmap = T,
            plot_scatter = T,
            standardize = F) 
  }else{
    if(is.numeric(metadata[,(which(colnames(metadata)==x))])){
      myMin=min(na.omit(metadata[,(which(colnames(metadata)==x))]))
      myMedian=median(na.omit(metadata[,(which(colnames(metadata)==x))]))
      myMax=max(na.omit(metadata[,(which(colnames(metadata)==x))]))
      print(paste("Numeric",x))
      if(myMin==myMedian){
        print(paste("Error with",x))
        myMedian=myMin+1
      }
      q=cut.default(metadata[,(which(colnames(metadata)==x))], breaks=c(myMin,
                                                                        myMedian,
                                                                        myMax), include.lowest = TRUE)
      print(paste("I fixed",x))
      # write.csv(filtered,"filteredTest.csv")
      # write.csv(metadata,"metadataTest.csv")
      metadata[,(which(colnames(metadata)==x))]=q
      if(length(which(is.na(metadata[,which(colnames(metadata)==x)])))>0){
        input=as.data.frame(counts.data[,-which(is.na(metadata[,which(colnames(metadata)==x)]))])
        thisMeta=as.data.frame(metadata[-which(is.na(metadata[,which(colnames(metadata)==x)])),(which(colnames(metadata)==x))])
        colnames(thisMeta)=x
        row.names(thisMeta)=metadata$ID[-which(is.na(metadata[,which(colnames(metadata)==x)]))]
      }else{
        input=as.data.frame(counts.data)
        thisMeta=as.data.frame(metadata[,(which(colnames(metadata)==x))])
        colnames(thisMeta)=x
        row.names(thisMeta)=metadata$ID
      }
      
      Maaslin2(input_data=input,
                               input_metadata=thisMeta,
                               output=paste("Tweedieverse",x,"GMGC",sep="_"),
                               max_significance = 0.05,
                               plot_heatmap = T,
                               plot_scatter = T,
                               standardize = F) 
    }
    else{
      print(paste("error with", x))
    }
  }
}

db=vroom("uniqMultiSampleMapping.txt")
#db=db[-which(!is.na(db$`...4`)),]
listFolders=list.files(pattern="Tweedieverse.*GMGC")
for(y in listFolders[1]){
  fullTable=vroom(paste(y,"/all_results.tsv",sep=""))
  input=fullTable
  dbSub=db
  #dbSub=db[which(db$GMGC %in% input$feature),]
  #dbSub=dbSub[-which(duplicated(dbSub)),]
  row.names(input)=input$feature
  input=input[which(input$feature %in% db$GMGC),]
  row.names(input)=input$feature[which(input$feature %in% db$GMGC)]
  colnames(dbSub)=c("GMGC","GO","NAMES")
  x=paste(str_remove_all(y,"Tweedieverse_"),"_Deepath",sep="")
  dbSub=as.data.frame(dbSub[,c(1:3)])
  deepath_results <- deepath(input_data=as.data.frame(input),
                             output = paste("Deepath_BP_",x,sep=""),
                             #mapper_file = as.data.frame(dbSub), 
                             mapper_file = as.data.frame(dbSub), 
                             #pathway_col = "Pathway",
                             pathway_col = "NAMES",
                             #feature_col = "symbol",
                             feature_col = "GMGC",
                             score_col = 'coef',
                             pval_threshold = 0.05,
                             method = 'ks',
                             min_member = 2,
                             do_plot = TRUE)
  myRM=as.character(paste("Deepath_BP_",x,"/figures/gg_enrichment_rank.RDS",sep=""))
  unlink(myRM)
  myRM2=as.character(paste("Deepath_BP_",x,"/figures/gg_enrichment_score.RDS",sep=""))
  unlink(myRM2)
}