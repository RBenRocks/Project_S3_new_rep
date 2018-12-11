kidney <- list.files(path = "C:/Users/xNesTea/python/project_sem_3/new_dir/data/liver/DE/")
liver <- list.files(path = "C:/Users/xNesTea/python/project_sem_3/new_dir/data/kidney/DE/")
fnames <- intersect(kidney, liver)
#fnames <- gsub("_piano.txt", "", files)
for(name in fnames){
  setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/kidney/DE/")
  drug1 <- read.delim(paste(name))
  setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/liver/DE/")
  drug2 <- read.delim(paste(name))
  drug1$adj.P.Val[drug1$adj.P.Val>0.1] <- NaN
  drug2$adj.P.Val[drug2$adj.P.Val>0.1] <- NaN
  df <- merge(drug1, drug2, by = "X")
  df <- df[complete.cases(df),]
  df <- df[sign(df$logFC.x)==sign(df$logFC.y),]
  setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/common_de")
  name <- gsub(".txt", "", name)
  write.table(df, file = paste(name, "_intersect.txt", sep = ""), sep = "\t",col.names = NA)
  print(dim(df))
}