### set the working directory to the folder containing files to be normalized
setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/liver/")

### open metadata storing info about condition and cel file name
metadata <- read.delim("../liver.txt")

### Read necessary libraries
library(gcrma)
library(arrayQualityMetrics)

### read in cel files
Data <- ReadAffy()

### perform normalization of the data
eset <- gcrma(Data)

### perform Quality control
QC <- arrayQualityMetrics(expressionset=eset,outdir="normalized_l",force=TRUE,do.logtransform=TRUE)

### save outliers file names nto variable
outliers <- QC[["modules"]][["heatmap"]]@outliers@which

### if there are outliers
if(length(outliers) != 0){
  ### make new df from metadata containing only outliers 
  df1 <- metadata[which(metadata$Array.Data.File %in% names(outliers)),]
  ### save names of compounds present in outliers into variable
  o_compounds <- unique(df1$Parameter.Value.Compound.)
  ### save names of all files for compounds above
  o_files <- metadata[which(metadata$Parameter.Value.Compound. %in% o_compounds),]$Array.Data.File
  ### new dataframe storing only compounds which were not excluded
  df2 <- metadata[which(!(metadata$Array.Data.File %in% o_files)),]
  ### remove these files
  ### this is done since if we remove 1 cel file, we are left with just two and this is
  ### too litle to get significant results, so we can remove all files for this compound
  unlink(o_files)
}else{
  ### for constant variable names -> this can be simplified, test later on
  df2 <- metadata
}
### same steps as above for data without outliers, normalization and QC
Data <- ReadAffy()
eset <- gcrma(Data)
QC <- arrayQualityMetrics(expressionset=eset,outdir="normalized_after_outliers_l",force=TRUE,do.logtransform=TRUE)

### write new metadata, after outliers detection
write.table(df2, file = "./liver_after_o.txt",sep = "\t", row.names = F)


expres=as.data.frame(exprs(eset))
expres$PID=row.names(expres)
row.names(expres)=NULL
expres=expres[,c("PID",colnames(expres)[colnames(expres)!="PID"])]
p_to_g<- read.delim("../Array_resolved.txt")
#p_to_g$Gene.Symbol
expres <- merge(expres, p_to_g, by.x="PID", by.y="ID")
expres$mad <- apply(as.matrix(expres[,2:dim(df2)[1]+1]), 1, mad)
df<-expres[!(expres$Gene.Symbol == ""),]
unique = unique(df$Gene.Symbol)
finaldf <- data.frame()
for(gene in unique){
  temp <- df[df$Gene.Symbol == gene,]
  if (dim(temp)[1] > 1){
    temp <- temp[order(temp$mad, decreasing = TRUE),]
    temp <- cbind(temp$Gene.Symbol, temp[,1:dim(df2)[1]+1])
    finaldf <- rbind(temp[1,],finaldf)
  }
  else{
    temp <- cbind(temp$Gene.Symbol, temp[,1:dim(df2)[1]+1])
    finaldf <- rbind(temp,finaldf)
  }
}
finaldf1 <- finaldf[,-1]
rownames(finaldf1)<- finaldf[,1]
write.table(finaldf1, file = "liver_expres.txt",sep = "\t", row.names = F)
library(limma)
compounds <- unique(df2$Parameter.Value.Compound.)
for(compound in compounds){
  meta <- df2[df2$Parameter.Value.Compound. == compound,]
  cel_files <- meta[which(meta$Parameter.Value.Compound. %in% compound),]$Array.Data.File
  idx <- match(cel_files, names(finaldf1))
  intensity <- finaldf1[, idx]
  Group <- factor(meta$Parameter.Value.DoseLevel., levels=c("Control","High"))
  design <- model.matrix(~Group)
  colnames(design) <- c("Control", compound)
  fit <- lmFit(intensity, design)
  fit <- eBayes(fit)
  tab <- topTable(fit, n=Inf, coef=2)
  name <- as.character(compound)
  write.table(tab, file = paste(name, ".txt", sep = ""), sep = "\t",col.names = NA)
  #write.table(fit,results=NULL, paste(name, ".txt", sep = ""), sep="\t")
}

