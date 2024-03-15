##############################################
#                                            #
#  Rscript to combine output of htseq-count  #
#                       -Arindam Ghosh       #
#                        (7 June 2018)       #
#                                            #
##############################################
###I. emrge 전 dataframe 정리
#1. 필요한 라이브러리 설치
install.packages("data.table")
library(data.table)

#2. tsv 파일 vector 명 수정
getwd()
setwd("D:/이지연 DEG 분석/")

# 파일 리스트 생성
tsv_check <- read.table(file.choose(), sep="\t")
tsv_check
rename_list <- list.files(pattern="*.tsv")

# 각 파일을 불러와서 열 이름 변경 후 다시 저장
for (file in rename_list) {
  # 파일 불러오기
  data <- read.table(file, header=FALSE, sep="\t")
  
  # 열 이름 변경
  colnames(data)[colnames(data)=="V1"] <- "gene_name"
  colnames(data)[colnames(data)=="V2"] <- "gene_id"
  colnames(data)[colnames(data)=="V3"] <- "count"
  
  # 변경된 데이터 다시 저장
  write.table(data, file, sep="\t", row.names=FALSE, quote=FALSE)
}

#3. gene_id 열 삭제
# 파일 리스트 생성
del_list <- list.files(pattern="*.tsv")

# 각 파일을 불러와서 두 번째 열 삭제 후 다시 저장
# 필요한 패키지 로드
library(data.table)
library(dplyr)

# 'del_list'에 있는 각 파일에 대해 작업 수행
for (file in del_list) {
  
  # 파일 불러오기
  data <- fread(file)
  
  # 'gene_id' 열에서 값이 빈 문자열인 행 삭제
  cleaned_data <- data %>% filter(gene_id != "")
  
  # 두 번째 열 삭제
  cleaned_data <- cleaned_data[,-2]
  
  # 변경된 데이터 다시 저장
  new_file_name <- paste0("new_", basename(file))
  fwrite(cleaned_data, new_file_name, sep="\t")
}


-----------------------------------------------------------------------------------------------
  ###II. R script for looping through the htseq outputs and merge them into one matrix
setwd("[directory_path]")
getwd("[directory_path]")
path <- "D:/이지연 DEG 분석/HTseq_gene_name"
myoutname <- "merged_counts"

# List all files in the specified path that start with 'new_' and end with '.tsv'
files <- list.files(path=path, pattern="^new_.*\\.tsv$")

# Manipulate file names by removing 'new_' and '.tsv' from the name
labs <- gsub("new_", "", files)
labs <- gsub("\\.tsv", "", labs, perl=TRUE)

# Initialize an empty list to store the data frames
cov <- list()

# Loop through each file
for (i in labs) {
  # Construct the filepath
  filepath <- file.path(path, paste("new_", i, ".tsv", sep=""))
  
  # Read the data file into a data frame
  data <- read.table(filepath, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  
  # Aggregate data by 'gene_name', summing the values
  data <- aggregate(. ~ gene_name, data, sum)
  
  # Rename the columns of the data frame
  colnames(data) <- c("gene_name", i)
  
  # Add the data frame to the list
  cov[[i]] <- data
}

# Merge all data frames in the list into one data frame by 'gene_name'
df <- Reduce(function(x, y) merge(x = x, y = y, by ="gene_name", all = TRUE), cov)

# Write the merged data frame to a CSV file
write.csv(df, file = paste(path, myoutname, ".csv", sep=""), row.names = FALSE)
