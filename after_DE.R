kidney <- list.files(path = "C:/Users/xNesTea/python/project_sem_3/new_dir/data/liver/piano/")
liver <- list.files(path = "C:/Users/xNesTea/python/project_sem_3/new_dir/data/kidney/piano/")
fnames <- intersect(kidney, liver)
#fnames <- gsub("_piano.txt", "", files)
for(name in fnames){
  setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/kidney//piano/")
  drug1 <- read.delim(paste(name))
  break}
  setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/liver/piano/")
  drug2 <- read.delim(paste(name,".txt",sep=""))
  drug1$adj.P.Val[drug1$p.adj..dist.dir.up.>0.05] <- NaN
  drug2$adj.P.Val[drug2$p.adj..dist.dir.up.>0.05] <- NaN
  df <- merge(drug1, drug2, by = "X")
  df <- df[complete.cases(df),]
  df <- df[sign(df$logFC.x)==sign(df$logFC.y),]
  setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/common_de")
  write.table(df, file = paste(name, ".txt", sep = ""), sep = "\t",col.names = NA)
}
kidney <- list.files(path = "C:/Users/xNesTea/python/project_sem_3/For_Uppmax/liver/")
metadata <- read.delim("")