install.packages("data.table")
library(data.table)

getwd()
setwd("/mnt/sdb1/hanwoo_mRNA_DEG/umd3.1.1_star2pass_htseq_count/gene_count")

tsv_check <- read.table(file.choose(), sep="\t")
tsv_check
rename_list <- list.files(pattern="*.tsv")

for (file in rename_list) {
  data <- read.table(file, header=FALSE, sep="\t")
  df <- data[-c(1,2,3,4),]
  colnames(df)[colnames(df)=="V1"] <- "gene_name"
  colnames(df)[colnames(df)=="V2"] <- "gene_id"
  colnames(df)[colnames(df)=="V3"] <- "count"
  write.table(df, file, sep="\t", row.names=FALSE, quote=FALSE)
  }

del_list <- list.files(pattern="*.tsv")

tsv_check <- read.table(file.choose(), sep="\t")
tsv_check

library(data.table)
library(dplyr)

list.files()

for (file in del_list) {
  df2 <- fread(file)
  cleaned_data <- df2 %>% filter(gene_id != "")
  cleaned_data <- cleaned_data[,-2]
  new_file_name <- paste0("new_", basename(file))
  fwrite(cleaned_data, new_file_name, sep="\t")
  }

getwd()
list.files()

path <- "/mnt/sdb1/hanwoo_mRNA_DEG/umd3.1.1_star2pass_htseq_count/gene_count"
myoutname <- "merged_counts"

df_files <- list.files(path=path, pattern="^new_.*\\.tsv$")
df_files

labs <- gsub("new_", "", df_files)
labs <- gsub("\\.tsv", "", labs, perl=TRUE)

labs

cov <- list()

for (i in labs) {
  filepath <- file.path(path, paste("new_", i, ".tsv", sep=""))
  data3 <- read.table(filepath, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  data4 <- aggregate(. ~ gene_name, data3, sum)
  colnames(data4) <- c("gene_name", i)
  cov[[i]] <- data4
  }

data4

df <- Reduce(function(x, y) merge(x = x, y = y, by ="gene_name", all = TRUE), cov)

write.csv(df, file = paste(path, myoutname, ".csv", sep=""), row.names = FALSE)
list.files()
