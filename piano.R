setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/liver/")
library("piano")
DE <- read.delim("./doxorubicin.txt",row.names = 1)
head(as.matrix(DE[,1]))

DE=DE[ ,c('logFC','adj.P.Val')]
pval= as.matrix(DE[, 2]) #extract P as a matrix
fc= as.matrix(DE[, 1])  #extract fold changes as a matrix
row.names(pval)=row.names(DE)
row.names(fc)=row.names(DE)

gset=loadGSC("../Rattus_norvegicus_GSEA_GO_sets_all_symbols_April_2015.gmt")
gsaRes <- runGSA(pval,fc,gsc=gset, nPerm = 1000, adjMethod = "fdr")
GSAsummaryTable(gsaRes)

setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/liver/")
library("piano")
kidney <- list.files(path = "C:/Users/xNesTea/python/project_sem_3/new_dir/data/liver/")
liver <- list.files(path = "C:/Users/xNesTea/python/project_sem_3/new_dir/data/kidney/")
fnames <- intersect(kidney, liver)
fnames <- gsub(".txt", "", fnames)
gset=loadGSC("../Rattus_norvegicus_GSEA_GO_sets_all_symbols_April_2015.gmt")
for(name in fnames){
  setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/kidney/")
  DE <- read.delim(paste(name,".txt",sep=""),row.names = 1)
  DE=DE[ ,c('logFC','adj.P.Val')]
  pval= as.matrix(DE[, 2]) #extract P as a matrix
  fc= as.matrix(DE[, 1])  #extract fold changes as a matrix
  row.names(pval)=row.names(DE)
  row.names(fc)=row.names(DE)
  gsaRes <- runGSA(pval,fc,gsc=gset, nPerm = 1000, adjMethod = "fdr")
  GSAsummaryTable(gsaRes, save = TRUE, paste(name,"_piano.txt",sep=""))
}