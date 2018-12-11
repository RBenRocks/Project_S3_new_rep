setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/liver/DE/")
library("piano")
kidney <- list.files(path = "C:/Users/xNesTea/python/project_sem_3/new_dir/data/liver/DE/")
liver <- list.files(path = "C:/Users/xNesTea/python/project_sem_3/new_dir/data/kidney/DE/")
fnames <- intersect(kidney, liver)
fnames <- gsub(".txt", "", fnames)
gset=loadGSC("../../Rattus_norvegicus_GSEA_GO_sets_all_symbols_April_2015.gmt")
for(name in fnames){
  setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/liver/DE/")
  DE <- read.delim(paste(name,".txt",sep=""),row.names = 1)
  DE=DE[ ,c('logFC','P.Value')]
  pval= as.matrix(DE[, 2]) #extract P as a matrix
  fc= as.matrix(DE[, 1])  #extract fold changes as a matrix
  row.names(pval)=row.names(DE)
  row.names(fc)=row.names(DE)
  gsaRes <- runGSA(pval,fc,gsc=gset, nPerm = 1000, adjMethod = "fdr")
  name <- paste0("../piano/", name)
  GSAsummaryTable(gsaRes, save = TRUE, paste(name,"_piano.txt",sep=""))
}