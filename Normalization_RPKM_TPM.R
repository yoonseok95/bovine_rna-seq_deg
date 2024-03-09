#I. gene length 구하기
getwd()
setwd("[directory_path]")

#paskage download
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicFeatures")
update.packages()

#first get the gtf file from, lets say, 
library(GenomicFeatures)

ncbidb <- makeTxDbFromGFF("[directory_path]/genomic.gtf",format="gtf")
head(ncbidb)
View(ncbidb)

# then collect the exons per gene id
exons.list.per.gene <- exonsBy(ncbidb,by="gene")
head(exons.list.per.gene)
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
head(exonic.gene.sizes)
#to see them in a table format, you should unlist them
unlist_geneLength<-unlist(exonic.gene.sizes)
write.table(unlist_geneLength,"23.11.19_genelength calculation.txt")

--------------------------------------------------------------------------------------------------------
  
#II. Normalization (TPM)
# 데이터 로드
gene_lengths <- read.table("23.11.19_genelength calculation.txt",  header=T, quote = "")
gene_counts <- read.csv("merged 24.csv", row.names=1, header=T)

# gene_lengths의 유전자 이름에서 따옴표 제거
row.names(gene_lengths) <- gsub("\"", "", row.names(gene_lengths))

# RPKM 계산 함수
rpkm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(counts) * 1e6
}

# TPM 계산 함수
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

# RPKM/TPM 계산
rpkms <- apply(gene_counts, 2, function(x) rpkm(as.numeric(x), gene_lengths))
tpms <- apply(gene_counts, 2, function(x) tpm(as.numeric(x), gene_lengths))

tpms_df <- do.call(cbind, tpms)
rpkms_df <- do.call(cbind, rpkms)

colnames(tpms_df) <- colnames(gene_counts)
colnames(rpkms_df) <- colnames(gene_counts)

write.csv(rpkms_df, "rpkm_results.csv")
write.csv(tpms_df, "tpm_results.csv")


