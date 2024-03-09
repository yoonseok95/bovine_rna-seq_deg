################################################################################
#
# Course 2. DEG analysis & visualization
#
################################################################################
# Step 2.1 Data load from RINE2020 results.
################################################################################


getwd()
setwd("E:/이지연 DEG 분석/")
data = read.csv("merged 24.csv", row.names=1)
colnames(data) # To check the column names of the data



################################################################################
# Step 2.2 Extracting gene expression information from HTseq count results.
################################################################################

# Create a subset of the 'data' dataframe with only the columns corresponding to ReadCount
rc_table = data


dim(rc_table) # Print the dimensions of the resulting 'rc_table'
head(rc_table) # Print the first few rows of the 'rc_table'
summary(rc_table) # Generate summary statistics for each column of the 'rc_table'

# Change the column names of the 'rc_table' to more user-friendly names
colnames(rc_table)
colnames(rc_table) <- c("H_2nd_12_9", "H_2nd_16_9", "H_2nd_2_9",  "H_2nd_8_9", 
                        "H_3rd_14_9", "H_3rd_16_8", "H_3rd_17_9", "H_4th_14_9",
                        "H_5th_10_9", "H_5th_9_9",  "L_2nd_18_4", "L_2nd_19_4",
                        "L_3rd_20_1", "L_3rd_4_2",  "L_3rd_6_3",  "L_3rd_7_3" ,
                        "L_3rd_8_2",  "L_4th_15_3", "L_4th_17_1", "L_5th_19_1",
                        "H_RNAlater_6_7",  "H_RNAlater_7_7",  "L_RNAlater_8_4",  "L_RNAlater_9_5")
colnames(rc_table)

# Open a new tabular view of the 'rc_table'
View(rc_table)

################################################################################
# Step 2.3 DEG analysis using TCC
################################################################################

# install TCC 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCC")
library(TCC)

# initiation & filter low count genes
group <- c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,1,1,2,2)
tcc <- new("TCC", rc_table, group)
dim(tcc$count)
tcc
tcc <- filterLowCountGenes(tcc, low.count = 10)
dim(tcc$count)
tcc

# normalization & DEG estimation 
type <- "deseq2"
if (type == "edger"){
  tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger",
                         iteration=3, FDR=0.1, floorPDEG=0.05)
  tcc <- estimateDE(tcc, test.method="edger", FDR=0.1)
}else if (type == "deseq2"){
  tcc <- calcNormFactors(tcc, norm.method="deseq2", test.method="deseq2",
                         iteration=3, FDR=0.1, floorPDEG=0.05)
  tcc <- estimateDE(tcc, test.method="deseq2", FDR=0.1)
}
table(tcc$estimatedDEG)



#######2_RNAlater 샘플 제거
rc_table = data
dim(rc_table)
rc_table2 = subset(data, select=c("hw_H_2nd_12_htseq.count", "hw_H_2nd_16_htseq.count", "hw_H_2nd_2_htseq.count",  "hw_H_2nd_8_htseq.count", 
                                    "hw_H_3rd_14_htseq.count", "hw_H_3rd_16_htseq.count", "hw_H_3rd_17_htseq.count", "hw_H_4th_14_htseq.count",
                                    "hw_H_5th_10_htseq.count", "hw_H_5th_9_htseq.count",  "hw_L_2nd_18_htseq.count", "hw_L_2nd_19_htseq.count",
                                    "hw_L_3rd_20_htseq.count", "hw_L_3rd_4_htseq.count",  "hw_L_3rd_6_htseq.count",  "hw_L_3rd_7_htseq.count", 
                                    "hw_L_3rd_8_htseq.count",  "hw_L_4th_15_htseq.count", "hw_L_4th_17_htseq.count", "hw_L_5th_19_htseq.count"))
colnames(rc_table2) <- c("H_2nd_12_9", "H_2nd_16_9", "H_2nd_2_9",  "H_2nd_8_9", 
                           "H_3rd_14_9", "H_3rd_16_8", "H_3rd_17_9", "H_4th_14_9",
                           "H_5th_10_9", "H_5th_9_9",  "L_2nd_18_4", "L_2nd_19_4",
                           "L_3rd_20_1", "L_3rd_4_2",  "L_3rd_6_3",  "L_3rd_7_3" ,
                           "L_3rd_8_2",  "L_4th_15_3", "L_4th_17_1", "L_5th_19_1")
colnames(rc_table2)  


# initiation & filter low count genes
group <- c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2)
tcc <- new("TCC", rc_table2, group)
dim(tcc$count)
tcc
tcc <- filterLowCountGenes(tcc, low.count = 10)
dim(tcc$count)
tcc

# normalization & DEG estimation 
type <- "deseq2"
if (type == "edger"){
  tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger",
                         iteration=3, FDR=0.1, floorPDEG=0.05)
  tcc <- estimateDE(tcc, test.method="edger", FDR=0.1)
}else if (type == "deseq2"){
  tcc <- calcNormFactors(tcc, norm.method="deseq2", test.method="deseq2",
                         iteration=3, FDR=0.1, floorPDEG=0.05)
  tcc <- estimateDE(tcc, test.method="deseq2", FDR=0.1)
}
table(tcc$estimatedDEG)


# The getNormalizedData function can be applied to the TCC class object
# after the normalization factors have been calculated
eff_count <- getNormalizedData(tcc)
View(eff_count)

write.csv(round(eff_count, 3), "tmm(영주제외20_low.count20).csv")

# The summary statistics for top-ranked genes
result <- getResult(tcc, sort=TRUE)
dim(result)

# DEG filtering (log2fc가 2이상, -2이하인 유전자만 DEG로 선정하도록 데이터 수정)
result$estimatedDEG[abs(result$m.value) < 2] <- 0
View(result)
table(result$estimatedDEG)

#tsv로 저장
write.table(result, "deg.tsv", sep="\t", quote=F, col.names=T, row.names=F)

#csv로 저장
write.csv(result, "deg(영주제외20_lowcount40).csv")




################################################################################
# Step 2.4 DEG visualization
################################################################################

# MA plot v1
png("deg.ma.png", width=8, height=6, units="in", res=600)
plot(tcc, median.lines = TRUE, cex=0.4)
dev.off()

# MA plot v2
result_deg_up <- result[result$m.value > 2 & result$estimatedDEG > 0,]
dim(result_deg_up)
result_deg_dw <- result[result$m.value < -2 & result$estimatedDEG > 0,]
dim(result_deg_dw)

png("deg.ma_v2.png", width=8, height=6, units="in", res=600)
plot(result$a.value, result$m.value, cex=0.3, ylim=c(-10,10), xlim=c(-3, 15),
     main = "MA plot", xlab = "A value", ylab = "M value")
points(result_deg_up$a.value, result_deg_up$m.value, cex=0.3, col='red')
points(result_deg_dw$a.value, result_deg_dw$m.value, cex=0.3, col='blue')
dev.off()

# Volcano plot
png("deg.volcano.png", width=8, height=6, units="in", res=600)
plot(result$m.value, -log(result$p.value,10), cex=0.3, ylim=c(0,8), xlim=c(-10,10))
points(result_deg_up$m.value, -log(result_deg_up$p.value,10), cex=0.3, col='red')
points(result_deg_dw$m.value, -log(result_deg_dw$p.value,10), cex=0.3, col='blue')
dev.off()

png("deg.volcano2.png", width=8, height=6, units="in", res=600)
plot(result$m.value, -log(result$p.value,10), cex=0.3, ylim=c(0,8), xlim=c(-10,10))
points(result_deg_up$m.value, -log(result_deg_up$p.value,10), cex=0.3, col='red')
points(result_deg_dw$m.value, -log(result_deg_dw$p.value,10), cex=0.3, col='blue')

data$gene_id <- rownames(data)
data$GeneSymbol <- paste(data$gene_id, "(", data$Symbol, ")")
result_anno <- merge(result, data, by = "gene_id", all.x = TRUE)
result_rank <- result_anno[result_anno$rank <= 10,]
text(result_rank$m.value, -log(result_rank$p.value,10), result_rank$GeneSymbol, col = 'mediumpurple4')
dev.off()

result_anno[result_anno$rank <= 10,]

################################################################################

