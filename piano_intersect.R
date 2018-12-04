kidney <- list.files(path = "C:/Users/xNesTea/python/project_sem_3/new_dir/data/liver/piano/")
liver <- list.files(path = "C:/Users/xNesTea/python/project_sem_3/new_dir/data/kidney/piano/")
fnames <- intersect(kidney, liver)
#fnames <- gsub("_piano.txt", "", fnames)
for(name in fnames){
  setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/kidney/piano/")
  drug1 <- read.delim(paste(name))
  drug1[is.na(drug1)] <- 1
  setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/liver/piano/")
  drug2 <- read.delim(paste(name))
  drug2[is.na(drug2)] <- 1
  drug1$p.adj..dist.dir.up.[drug1$p.adj..dist.dir.up.>0.1] <- NaN
  drug2$p.adj..dist.dir.up.[drug2$p.adj..dist.dir.up.>0.1] <- NaN
  df <- merge(drug1, drug2, by = "Name")
  df <- df[complete.cases(df),]
  n <- gsub("_piano.txt", "", name)
  #df <- df[sign(df$logFC.x)==sign(df$logFC.y),]
  setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/10padj/UP_gse/")
  write.table(df, file = paste(n, ".txt", sep = ""), sep = "\t",col.names = NA)
  print(dim(df))
}
kidney <- list.files(path = "C:/Users/xNesTea/python/project_sem_3/For_Uppmax/liver/")
metadata <- read.delim("")

kidney <- list.files(path = "C:/Users/xNesTea/python/project_sem_3/new_dir/data/liver/piano/")
liver <- list.files(path = "C:/Users/xNesTea/python/project_sem_3/new_dir/data/kidney/piano/")
fnames <- intersect(kidney, liver)
#fnames <- gsub("_piano.txt", "", fnames)
for(name in fnames){
  setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/kidney/piano/")
  drug1 <- read.delim(paste(name))
  drug1[is.na(drug1)] <- 1
  setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/liver/piano/")
  drug2 <- read.delim(paste(name))
  drug2[is.na(drug2)] <- 1
  drug1$p.adj..non.dir..[drug1$p.adj..non.dir..>0.1] <- NaN
  drug2$p.adj..non.dir..[drug2$p.adj..non.dir..>0.1] <- NaN
  df <- merge(drug1, drug2, by = "Name")
  df <- df[complete.cases(df),]
  n <- gsub("_piano.txt", "", name)
  #df <- df[sign(df$logFC.x)==sign(df$logFC.y),]
  setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/10padj/NDIR_gse/")
  write.table(df, file = paste(n, ".txt", sep = ""), sep = "\t",col.names = NA)
  print(dim(df))
}

kidney <- list.files(path = "C:/Users/xNesTea/python/project_sem_3/new_dir/data/liver/piano/")
liver <- list.files(path = "C:/Users/xNesTea/python/project_sem_3/new_dir/data/kidney/piano/")
fnames <- intersect(kidney, liver)
#fnames <- gsub("_piano.txt", "", fnames)
for(name in fnames){
  setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/kidney/piano/")
  drug1 <- read.delim(paste(name))
  drug1[is.na(drug1)] <- 1
  setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/liver/piano/")
  drug2 <- read.delim(paste(name))
  drug2[is.na(drug2)] <- 1
  drug1$p.adj..dist.dir.dn.[drug1$p.adj..dist.dir.dn.>0.1] <- NaN
  drug2$p.adj..dist.dir.dn.[drug2$p.adj..dist.dir.dn.>0.1] <- NaN
  df <- merge(drug1, drug2, by = "Name")
  df <- df[complete.cases(df),]
  n <- gsub("_piano.txt", "", name)
  #df <- df[sign(df$logFC.x)==sign(df$logFC.y),]
  setwd("C:/Users/xNesTea/python/project_sem_3/new_dir/data/10padj/DOWN_gse/")
  write.table(df, file = paste(n, ".txt", sep = ""), sep = "\t",col.names = NA)
  print(dim(df))
}